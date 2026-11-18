# -*- coding: utf-8 -*-
"""verify_jcp.py -- the checks of the manuscript prepared for The Journal of Chemical Physics (deposit 3.5.0) that
src/verify.py had no check for, or whose first-order form the JCP text made wrong.

Origin.  The JCP package checked its sources with a guard of the same name (2026-11-08 to 2026-11-18): it copied
deposit 3.4.1 into a temporary folder, laid the JCP sources over its paper/, ran src/verify.py there, accepted only a
registered list of failures (each the first-order form of a number the new text re-scores, or a passage that moved),
and ran a replacement check for each.  In deposit 3.5.0 paper/ IS the JCP text, so the copy, the overlay and the
registry are gone: src/verify.py was brought to the JCP text (CHECK_CHANGES of 3.5.0, docs/KNOWN_DISCREPANCIES.md
section 36) and runs every check below as its section (12), each counted with its assertions.  The checks are the
package's, logic unchanged; only the paths are the repository's:
  data/                              the deposited data, read only;
  data/c3_frontier/second_order/     the re-scored outputs of src/frontier/rescore_second_order.py and
                                     src/frontier/null_ws7_second_order.py;
  paper/                             the sources of the manuscript and of its Supplement;
  src/jcp/number_sources.json        the registry of 706 printed numbers with their stored sources;
  src/jcp/fracgrid_numbers.py        the numbers of Table tab:fracgrid, from data/eta_sweep.json.
Dropped from the package's guard: check_dequant_moved (the dequantization paragraph moved to sec_9_sm.tex, and
src/verify.py now checks it there against data/gflow.json).

Usage (from the repository root):  python src/verify_jcp.py      (these checks alone; python src/verify.py runs
them too, as its section 12).  Exit 0 iff every check passes.  Nothing is written.

The docstrings below are the package's, dated; "the deposit" in them is this repository, and "release/" or
"verify.py section 9.9c / 9.6c" refer to deposit 3.4.1, whose checks these replace.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))                      # <repo>/src
RELEASE = os.path.dirname(HERE)                                         # <repo> (the name is the package's)
JCP_PAPER = os.path.join(RELEASE, 'paper')
DATA = os.path.join(RELEASE, 'data')

# 2026-11-16 (rules D1-D23 of the change, docs/KNOWN_DISCREPANCIES.md section 36): Theorem 1 gains the second-order
# branch.  For a real probe the certificate is min{1+w, (sqrt(1-w) + Lambda^/eta)^2} (D6); for the complex momentum
# probes it stays the first-order branch (D8); a leakage value quoted above the trivial bound is the smaller of the two
# branches (D17).  The checks below re-derive every number of that change from STORED primitives (w_S, Lambda^/eta and
# the true errors of data/) and compare them with the re-scored outputs, data/c3_frontier/second_order/<name>.
SECOND_ORDER = 'c3_frontier/second_order'

# assertions of the last check run, read by src/verify.py section 12 (one check, n assertions)
LAST_N = [0]


def _count(n):
    LAST_N[0] = int(n)


def _data_file(rel):
    """A deposited data file by its path relative to data/."""
    return os.path.join(DATA, *rel.split('/'))


def _so(name):
    """One re-scored output (cells.json, v2.json, ...), from data/c3_frontier/second_order/."""
    import io
    import json
    return json.load(io.open(_data_file(SECOND_ORDER + '/' + name), encoding='utf-8'))


def _B1(w, x):
    """First-order branch relative to ||phi||^2: (1+sqrt w)(sqrt(1-w)+x), x = Lambda^/eta."""
    import math
    return (1.0 + math.sqrt(w)) * (math.sqrt(max(1.0 - w, 0.0)) + x)


def _B2(w, x):
    """Second-order branch relative to ||phi||^2 for a real probe (Lambda~ = Lambda): (sqrt(1-w)+x)^2."""
    import math
    return (math.sqrt(max(1.0 - w, 0.0)) + x) ** 2


def _tau2(w):
    """The certificate's threshold for a real probe: B2 < 1+w exactly when x < sqrt(1+w) - sqrt(1-w)."""
    import math
    return math.sqrt(1.0 + w) - math.sqrt(max(1.0 - w, 0.0))


def _rnd(v, n, mode='near', kind='sf'):
    """v as a string, rounded to n significant figures (kind 'sf') or n decimals ('dp'): 'up' for a ceiling, 'down'
    for a floor, 'near' for a plain value or a range of stored values (working rule 7b); decimal arithmetic on repr(v)."""
    from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_UP
    d = Decimal(repr(float(v)))
    if d == 0:
        return '0'
    e = (d.adjusted() - (n - 1)) if kind == 'sf' else -n
    return format(d.quantize(Decimal(1).scaleb(e), rounding={'up': ROUND_CEILING, 'down': ROUND_FLOOR,
                                                               'near': ROUND_HALF_UP}[mode]), 'f')


def _tsci(v, n, mode='near'):
    """The LaTeX '$m\\times10^{e}$' of v with an n-figure mantissa, rounded in the given direction."""
    from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_UP
    r = {'up': ROUND_CEILING, 'down': ROUND_FLOOR, 'near': ROUND_HALF_UP}[mode]
    d = Decimal(repr(float(v)))
    e = d.adjusted()
    q = Decimal(1).scaleb(-(n - 1))
    m = d.scaleb(-e).quantize(q, rounding=r)
    if m >= 10:
        m, e = (m / 10).quantize(q, rounding=r), e + 1
    return '$%s\\times10^{%d}$' % (format(m, 'f'), e)


def _interp(x, xs, ys):
    """Linear interpolation of ys(xs) at x, xs increasing (numpy.interp inside the range); None outside it."""
    if not xs or x < xs[0] or x > xs[-1]:
        return None
    for a, b, ya, yb in zip(xs, xs[1:], ys, ys[1:]):
        if a <= x <= b:
            return ya if b == a else ya + (yb - ya) * (x - a) / (b - a)
    return ys[-1]


def _grid_rows():
    """C3_grid.json and its rows with both branches recomputed from the stored primitives (w_S, Lam_over_eta); the
    stored first-order bound_leak is reproduced first (verify.py section 10b checks the same identity)."""
    import io
    import json
    g = json.load(io.open(os.path.join(DATA, 'c3_frontier', 'C3_grid.json'), encoding='utf-8'))
    rows = []
    for i, r in enumerate(g['grid']):
        b1, b2 = _B1(r['w_S'], r['Lam_over_eta']), _B2(r['w_S'], r['Lam_over_eta'])
        if abs(b1 - r['bound_leak']) > 1e-12 * max(1.0, r['bound_leak']):
            raise AssertionError('C3_grid.json /grid[%d]: bound_leak is not B1(w_S, Lam_over_eta)' % i)
        rows.append(dict(r, idx=i, B1=b1, B2=b2, T=r['bound_trivial'], c2=min(r['bound_trivial'], b2)))
    return g, rows


def _fsum_frontier(Bs, FRs, th):
    """src/frontier/fsum.py frontier(): the first grid fraction (FR order) at which B < th, log-linear in FR
    of log(B/th) between the bracketing fractions.  Returns (FR*, lo, hi); (nan, last, nan) if never below."""
    import math
    s = sorted(zip(FRs, Bs, th))
    for i, (f, b, t) in enumerate(s):
        if b < t:
            if i == 0:
                return f, float('nan'), f
            f0, b0, t0 = s[i - 1]
            y0, y1 = math.log(b0 / t0), math.log(b / t)
            return f0 + (f - f0) * (0.0 - y0) / (y1 - y0), f0, f
    return float('nan'), s[-1][0], float('nan')


def _second_order_cells():
    """The sixteen (L, eta) cells of C3_grid.json 'frontier' with the certificate of D6.  Per cell: FR* (fsum.py
    convention, which first reproduces the stored first-order FR_star), its bracket and the convergence of its ends;
    the true error at FR* (log-linear in FR_real over the rows of the cell, the convention of check_threshold_error
    since 2026-11-08) and the same at the stored first-order FR_star; bound over error at the crossing ((1+w) there,
    w linear in FR over the CONVERGED rows, as calib_oos.py, over that error); 1-w at FR* (log-linear in FR over the
    CONVERGED rows); the least loose certified CONVERGED grid point; the margin over the reference drift of the
    published-subspace file of the same L (pt_L*_PUB.json meta.ref_drift).  Every value is compared with the
    re-scored second_order/cells.json (deposit file SECOND_ORDER/cells.json).  Returns (cells, data failures)."""
    import glob
    import io
    import json
    import math
    g, rows = _grid_rows()
    drift = {}
    for fn in glob.glob(os.path.join(DATA, 'c3_frontier', 'published_fraction', 'pt_L*_PUB.json')):
        x = json.load(io.open(fn, encoding='utf-8'))
        x = x[0] if isinstance(x, list) else x
        drift[x['meta']['L']] = x['meta']['ref_drift']
    ref = {(c['L'], c['eta']): c for c in _so('cells.json')['cells']}
    cells, bad = [], []

    def err_at(cl, fs):
        s = sorted(cl, key=lambda r: r['FR_real'])
        for a, b in zip(s, s[1:]):
            if a['FR_real'] <= fs <= b['FR_real']:
                t = (fs - a['FR_real']) / (b['FR_real'] - a['FR_real'])
                return math.exp((1 - t) * math.log(a['relL1']) + t * math.log(b['relL1']))
        return None
    for f in g['frontier']:
        L, eta = f['L'], f['eta']
        cl = [r for r in rows if r['L'] == L and abs(r['eta'] - eta) < 1e-12]
        one = _fsum_frontier([r['B1'] for r in cl], [r['FR'] for r in cl], [r['T'] for r in cl])
        if abs(one[0] - f['FR_star']) > 1e-12:
            bad.append('L=%d eta=%.2f: the fsum.py convention does not reproduce the stored FR_star' % (L, eta))
        fs, lo, hi = _fsum_frontier([r['B2'] for r in cl], [r['FR'] for r in cl], [r['T'] for r in cl])
        rlo = [r for r in cl if r['FR'] == lo][0]
        rhi = [r for r in cl if r['FR'] == hi][0]
        cv = sorted((r for r in cl if r['conv'] == 'CONVERGED' and r['relL1']), key=lambda r: r['FR'])
        wC = _interp(fs, [r['FR'] for r in cv], [r['w_S'] for r in cv])
        lomw = _interp(fs, [r['FR'] for r in cv], [math.log(max(1.0 - r['w_S'], 1e-300)) for r in cv])
        err = err_at(cl, fs)
        nt = [r['c2'] / r['relL1'] for r in cl if r['B2'] < r['T'] and r['relL1'] and r['conv'] == 'CONVERGED']
        c = dict(L=L, eta=eta, FRstar=fs, lo=lo, hi=hi, conv=(rlo['conv'], rhi['conv']), err=err,
                 FRstar_old=f['FR_star'], err_old=err_at(cl, f['FR_star']), ratio=(1.0 + wC) / err,
                 omw=math.exp(lomw), least=min(nt), margin=err / drift[L]['%g' % eta],
                 FR05=f['FR_true']['0.05']['FR'])
        cells.append(c)
        s = ref.get((L, eta))
        if s is None:
            bad.append('second_order/cells.json has no cell L=%d eta=%.2f' % (L, eta))
            continue
        for k, k2 in (('FRstar', 'FRstar_new'), ('err', 'err_new'), ('ratio', 'ratio_new'),
                      ('least', 'least_loose_new_conv'), ('margin', 'drift_margin_new_PUB'),
                      ('omw', 'omw_new_logFRconv')):
            if abs(c[k] - s[k2]) > 1e-9 * abs(s[k2]):
                bad.append('L=%d eta=%.2f: %s recomputed %.12g, second_order/cells.json %s %.12g' % (L, eta, k, c[k], k2, s[k2]))
        if [lo, hi] != s['br_new']:
            bad.append('L=%d eta=%.2f: bracket [%g, %g], second_order/cells.json %r' % (L, eta, lo, hi, s['br_new']))
    return cells, bad


def _verdict(exp, label):
    """Compare one expectation (a dict of lists) with the printed text and print the result.  Keys:
      data_bad  failures found in the stored data themselves;
      want      (file, literal) that must be printed;
      gone      (file, literal) that must not be printed (the superseded form);
      regex     (file, pattern, value, what): at least one hit, and every hit equal to value;
      near      (file, anchor, words, before, after, what): at least one occurrence of anchor has every word within
                `before` characters before it and `after` after it ('re:' marks a regular expression);
      nsp_want, nsp_gone, nsp_any   the same for formulas, spaces and spacing macros ignored."""
    bad = list(exp.get('data_bad', []))
    cache = {}

    def T(f):
        if f not in cache:
            cache[f] = _printed(f)
        return cache[f]
    for f, lit in exp.get('want', []):
        if lit not in T(f):
            bad.append('%s does not print %r' % (f, lit))
    for f, lit in exp.get('gone', []):
        if lit in T(f):
            bad.append('%s still prints the superseded %r' % (f, lit))
    for f, pat, val, what in exp.get('regex', []):
        hits = re.findall(pat, T(f))
        if not hits:
            bad.append('%s does not print %s' % (f, what))
        for h in hits:
            if h != val:
                bad.append('%s prints %r for %s; the stored data give %r' % (f, h, what, val))
    for f, anchor, words, before, after, what in exp.get('near', []):
        t, ok, seen = T(f), False, False
        for m in re.finditer(re.escape(anchor), t):
            seen = True
            win = t[max(0, m.start() - before):m.end() + after]
            if all((re.search(w[3:], win) if w.startswith('re:') else w.lower() in win.lower()) for w in words):
                ok = True
                break
        if not seen:
            bad.append('%s does not print %r' % (f, anchor))
        elif not ok:
            bad.append('%s: %s (%r not within %d characters before / %d after %r)'
                       % (f, what, words, before, after, anchor))
    # formulas, compared with whitespace and the spacing macros \, \; \: \! removed (nsp_want / nsp_gone: a literal;
    # nsp_any: one of several literals)
    NS = lambda f: re.sub(r'\s+|\\[,;:!]', '', T(f))
    for f, lit in exp.get('nsp_want', []):
        if lit not in NS(f):
            bad.append('%s does not print %r (spaces ignored)' % (f, lit))
    for f, lit in exp.get('nsp_gone', []):
        if lit in NS(f):
            bad.append('%s still prints the superseded %r (spaces ignored)' % (f, lit))
    for f, alts, what in exp.get('nsp_any', []):
        if not any(a in NS(f) for a in alts):
            bad.append('%s does not print %s (one of %r, spaces ignored)' % (f, what, alts))
    for b in bad:
        print('  FAIL  %s' % b)
    print('  %s  %s' % ('ok  ' if not bad else 'FAIL', label))
    _count(1 + sum(len(exp.get(k, [])) for k in ('want', 'gone', 'regex', 'near', 'nsp_want', 'nsp_gone', 'nsp_any')))
    return not bad


def _expect_five_percent():
    """Data half of check_fraction_five_percent() (see there)."""
    import io
    g, rows = _grid_rows()
    exp = dict(data_bad=[], want=[], gone=[], regex=[], near=[])

    def first_below(L, eta, lev=0.05):
        cl = sorted((r for r in rows if r['L'] == L and abs(r['eta'] - eta) < 1e-12), key=lambda r: r['FR'])
        for i, r in enumerate(cl):
            if r['c2'] < lev:
                return (cl[i - 1] if i else None), r
        return cl[-1], None
    br = {}
    for L in (8, 10, 12):
        a, b = first_below(L, 0.18)
        if a is None or b is None:
            exp['data_bad'].append('L=%d eta=0.18: the certificate does not cross 0.05 inside the grid' % L)
            continue
        for r in (a, b):
            if abs(r['FR_real'] - r['FR']) >= 5e-3:
                exp['data_bad'].append('L=%d FR=%.2f: FR_real %.5f is not FR at two decimals' % (L, r['FR'], r['FR_real']))
        br[L] = (a, b)
    a8, b8 = br[8]
    if (a8['conv'], b8['conv']) != ('CONVERGED', 'CONVERGED'):
        exp['data_bad'].append('L=8 eta=0.18: the 0.05 bracket of the certificate is not two CONVERGED rows')
    for eta in (0.10, 0.15, 0.18):
        if first_below(6, eta)[1] is not None:
            exp['data_bad'].append('L=6 eta=%.2f: the certificate now reaches 0.05 on the grid (D19)' % eta)
    if first_below(6, 0.25)[1] is None:
        exp['data_bad'].append('L=6 eta=0.25: the certificate no longer reaches 0.05 on the grid (D19)')
    d19 = _so('v2.json')['D19']['cells']
    for L in (8, 10, 12):
        if L in br and d19['L%d_eta0.18' % L]['new'].get('bracket') != [br[L][0]['FR'], br[L][1]['FR']]:
            exp['data_bad'].append('L=%d eta=0.18: the bracket disagrees with second_order/v2.json D19' % L)
    lo, hi = '%.2f' % a8['FR'], '%.2f' % b8['FR']
    sec = r"between \$([\d.]+)\$ and \$([\d.]+)\$ of the sector at \$L=8\$"
    for rel in ('sec_1.tex', 'sec_5_body.tex', 'sec_10.tex', 'sm_frontier.tex'):
        exp['regex'].append((rel, sec, (lo, hi), '"between $lo$ and $hi$ of the sector at $L=8$"'))
        exp['gone'] += [(rel, 'about $0.98$ of the sector'), (rel, 'at least $0.98$ of the sector'),
                        (rel, '$[0.977,0.983]$'), (rel, 'not reached at all at $L=6$, $10$ or $12$'),
                        (rel, 'does not reach that criterion at $L=6$, $10$ or $12$')]
    for rel in ('sec_5_body.tex', 'sec_10.tex', 'sm_frontier.tex'):
        for L in (10, 12):
            if L in br:
                exp['want'].append((rel, 'between $%.2f$ and $%.2f$' % (br[L][0]['FR'], br[L][1]['FR'])))
        exp['regex'].append((rel, r"(\\eta)\s*\\le\s*0\.18", '\\eta', 'the qualifier of D19, \\eta\\le0.18'))
    ca, cb = _rnd(a8['c2'], 2, 'up'), _rnd(b8['c2'], 2, 'up')
    exp['want'] += [('sm_frontier.tex', '$%s$' % ca), ('sm_frontier.tex', '$%d$' % a8['S']),
                    ('sm_frontier.tex', '$%s$' % cb), ('sm_frontier.tex', '$%d$' % b8['S'])]
    exp['gone'] += [('sm_frontier.tex', '$0.084$ at $3802$'), ('sm_frontier.tex', '$0.060$ at $3830$'),
                    ('sm_frontier.tex', '$0.039$ at $3855$'), ('sm_frontier.tex', 'until $\\mathrm{FR}\\approx0.98$')]
    # the first-order fact behind the superseded wording, kept as a data check (labelled): the fine scan brackets the
    # FIRST-ORDER crossing of 0.05 at integer sizes
    zf = io.open(os.path.join(DATA, 'c3_frontier', 'adversarial', 'z_fine_L8.txt'), encoding='utf-8').read()
    fine = [(int(m.group(1)), float(m.group(2))) for m in
            re.finditer(r"\|S\|=(\d+) FR=[\d.]+ eta=0\.18 \| .*?relR5=([\d.e+-]+)", zf)]
    above = [s for s, v in fine if v >= 0.05]
    below = [s for s, v in fine if v < 0.05]
    if not (above and below and max(above) < min(below)):
        exp['data_bad'].append('z_fine_L8.txt no longer brackets the first-order crossing of 0.05')
    exp['info'] = dict(L8=(lo, hi, ca, a8['S'], cb, b8['S'], a8['conv'], b8['conv']),
                       L10=('%.2f' % br[10][0]['FR'], '%.2f' % br[10][1]['FR'], br[10][0]['conv'], br[10][1]['conv']),
                       L12=('%.2f' % br[12][0]['FR'], '%.2f' % br[12][1]['FR'], br[12][0]['conv'], br[12][1]['conv']),
                       first_order_fine=(max(above or [0]), min(below or [0])))
    exp['label'] = ('where the certificate falls below 0.05 at eta=0.18t: L=8 [%s, %s] (%s at %d, %s at %d of 3920, both '
                    'converged), L=10 [%s, %s], L=12 [%s, %s] (upper end %s); L=6 not reached at eta<=0.18t; first-order '
                    'fine bracket (%d, %d] kept as data'
                    % ((lo, hi, ca, a8['S'], cb, b8['S']) + exp['info']['L10'][:2] + exp['info']['L12'][:2]
                       + (exp['info']['L12'][3],) + exp['info']['first_order_fine']))
    return exp


def check_fraction_five_percent():
    """2026-11-08 (final check, repair B; until 2026-11-16 check_fraction_about_098), rewritten 2026-11-16 for the
    second-order branch (rule D6, D19; the number sheet rows I2, V22, V23, C6, C6b, F-2, F-31, F-32;
    the check of the number sheet section 5, first caution): where the certificate itself, min{1+w, B2} from the stored
    w_S and Lam_over_eta of C3_grid.json, first falls below the 0.05 acceptance criterion, in grid order.
      L=8, eta=0.18t: bracketed by two CONVERGED grid rows, printed 'between $lo$ and $hi$ of the sector at $L=8$' in
        sec_1, sec_5_body, sec_10 and sm_frontier (grid fractions at two decimals; |FR_real - FR| < 5e-3 is asserted,
        so the outward and the nearest rounding agree), and in sm_frontier the two rows as their certified values,
        rounded UP at two figures, and their sizes |S| of 3920.
      L=10 and L=12 at eta=0.18t: 'between $lo$ and $hi$' in sec_5_body, sec_10 and sm_frontier (the L=12 upper end is
        MARGINAL in depth, and both brackets hold at eta=0.18t only: the text must say so; not checked here).
      L=6: not reached at eta = 0.10, 0.15, 0.18t, reached at 0.25t (D19): '\\eta\\le0.18' must be printed in
        sec_5_body, sec_10 and sm_frontier.
    The superseded wording must be gone: 'about $0.98$' / 'at least $0.98$' of the sector, the bracket [0.977,0.983],
    the fine-scan values 0.084 / 0.060 / 0.039 at 3802 / 3830 / 3855, 'not reached at all at $L=6$, $10$ or $12$'.  The
    first-order fact behind the old wording is kept as a data check, labelled: z_fine_L8.txt still brackets the
    FIRST-ORDER crossing of 0.05 in (3830, 3855]/3920."""
    exp = _expect_five_percent()
    return _verdict(exp, exp['label'])


def _expect_fig_results():
    """Data half of check_fig_results() (see there): the plotted rows, the marks and the caption literals."""
    import io
    import json
    import math
    exp = dict(data_bad=[], want=[], gone=[], regex=[], near=[])
    g, grows = _grid_rows()
    rows = sorted((r for r in grows if r['L'] == 8 and abs(r['eta'] - 0.18) < 1e-12), key=lambda r: r['FR_real'])
    src = [(r['FR_real'], r['relL1'], min(r['B1'], r['B2']), r['bound_trivial'], 1.0 - r['w_S']) for r in rows]
    zf = io.open(os.path.join(DATA, 'c3_frontier', 'adversarial', 'z_fine_L8.txt'), encoding='utf-8').read()
    for m in re.finditer(r"\|S\|=(\d+) FR=[\d.]+ eta=0\.18 \| w_S=([\d.e+-]+) relR0=([\d.e+-]+) "
                         r"relR5=([\d.e+-]+) Lamhat/eta=([\d.e+-]+)", zf):
        w, b1, x = float(m.group(2)), float(m.group(4)), float(m.group(5))
        if abs(_B1(w, x) - b1) > 1e-6 * b1:
            exp['data_bad'].append('z_fine_L8.txt |S|=%s: relR5 is not the first-order branch of its w_S and Lamhat/eta'
                                   % m.group(1))
        src.append((int(m.group(1)) / 3920.0, float(m.group(3)), min(b1, _B2(w, x)), 1.0 + w, 1.0 - w))
    ref = _so('v2.json')['D17']['fig1']['rows']
    if len(ref) != len(src) or any(abs(r['leak_D17'] - s[2]) > 1e-9 * s[2] for r, s in zip(ref, src)):
        exp['data_bad'].append('the smaller-branch column does not agree with second_order/v2.json D17.fig1')
    sw = [r for r in json.load(io.open(os.path.join(DATA, 'support_witness.json'), encoding='utf-8'))['results']
          if r['L'] == 8][0]
    st = [t for t in sw['table'] if abs(t['eta'] - 0.18) < 1e-12][0]
    fr = [f for f in g['frontier'] if f['L'] == 8 and abs(f['eta'] - 0.18) < 1e-12][0]
    cells, cbad = _second_order_cells()
    exp['data_bad'] += cbad
    fs = [c for c in cells if c['L'] == 8 and abs(c['eta'] - 0.18) < 1e-12][0]['FRstar']
    cert = [min(s[2], s[3]) for s in src]
    i05 = [i for i in range(len(src)) if cert[i] < 0.05][0]
    lo05 = '%.2f' % (math.floor(src[i05 - 1][0] * 100 + 1e-9) / 100)
    hi05 = '%.2f' % (math.ceil(src[i05][0] * 100 - 1e-9) / 100)
    ranked56 = [s for s in src if abs(s[0] - 0.5599489796) < 1e-6][0]
    printed = (  # (literal in the caption, deposited value, printed value, tolerance on the printed rounding)
        ('falls from $0.40$ to', src[0][1], 0.40, 0.005),
        ('to $1.2\\times10^{-6}$', src[-1][1] * 1e6, 1.2, 0.05),
        ('falls from $20$ to', src[0][2], 20.0, 0.5),
        ('to $3.8\\times10^{-4}$', src[-1][2] * 1e4, 3.8, 0.05),
        ('$\\mathrm{FR}^{\\ast}=0.74$', fs, 0.74, 0.005),
        ('$\\lvert S\\rvert/D\\approx0.525$', fr['FR_true']['0.05']['FR'], 0.525, 0.0005),
        ('($\\lvert S\\rvert/D=0.571$)', sw['frac'], 0.571, 0.0005),
        ('$0.35$, is ten times', st['relL1'], 0.35, 0.005),
        ('nearly the same size ($0.560$)', ranked56[0], 0.560, 0.0005),
        ('is lower ($0.9988$)', 1.0 - ranked56[4], 0.9988, 0.00005),
    )
    F = 'figs/fig_results_native.tex'
    for lit, want, val, tol in printed:
        if abs(val - want) > tol:
            exp['data_bad'].append('caption literal %r: the deposit gives %r' % (lit, want))
        exp['want'].append((F, lit))
    exp['want'] += [(F, 'between $%s$ and $%s$' % (lo05, hi05)), (F, 'smaller leakage branch')]
    exp['gone'] += [(F, 'between $0.977$'), (F, 'and $0.983$'), (F, '$\\mathrm{FR}^{\\ast}=0.79$'), (F, 'to $0.039$'),
                    (F, 'and locate the crossing')]
    ratio = st['relL1'] / ranked56[1]
    if not 9.5 <= ratio <= 10.5:
        exp['data_bad'].append('"ten times": the deposit gives %.2f' % ratio)
    exp.update(src=src, star=(sw['frac'], st['relL1'], st['bound_theoremB'], st['bound_trivial']), FRstar=fs,
               FRtrue=fr['FR_true']['0.05']['FR'], bracket05=(lo05, hi05), n_caption=len(printed) + 1)
    return exp


def check_fig_results():
    """Every value of fig_results_L8.dat / fig_results_star.dat / fig_results_marks.tex against the
    deposit, and every number printed in the caption of fig:results against the same values.
    2026-11-16 (second-order branch; rule D14, D17, D22; the number sheet rows F1-F12, FS2-FS4): the plotted
    'leak' column is, row by row, the SMALLER of the two leakage branches -- the first-order one (stored bound_leak;
    relR5 of z_fine_L8.txt, checked against its w_S and Lamhat/eta) and the second-order one (sqrt(1-w)+Lambda^/eta)^2
    from the stored w_S and Lam_over_eta (Lamhat/eta) -- and agrees with second_order/v2.json D17.fig1; FigResFRstar is the
    certificate's crossing (_second_order_cells(), second_order/cells.json); the caption prints 20 (first row, first-order),
    3.8e-4 (last row, second-order, two figures), FR*=0.74 and the plotted rows between which the bound itself first
    falls below 0.05 (0.90 and 0.92, outward), names the curve 'smaller leakage branch', and no longer says that the
    last three points locate that crossing (F8).  The axis range must still contain every plotted value (F2), and
    the two labels at y=70 must straddle FR* within 0.02 of it (F3, F3b)."""
    import io
    exp = _expect_fig_results()
    src, bad = exp['src'], exp['data_bad']
    dat = [l.split() for l in io.open(os.path.join(JCP_PAPER, 'figs', 'fig_results_L8.dat'),
                                      encoding='utf-8').read().splitlines() if l and not l.startswith('#')][1:]
    if len(dat) != len(src) or len(src) != 18:
        bad.append('fig_results_L8.dat has %d rows, the deposit %d' % (len(dat), len(src)))
    for d, s in zip(dat, src):
        for x, y in zip(d[:5], s):
            if abs(float(x) - y) > 1e-8 * max(1.0, abs(y)):
                bad.append('fig_results_L8.dat value %s is not the deposited %r' % (x, y))
    star = [l.split() for l in io.open(os.path.join(JCP_PAPER, 'figs', 'fig_results_star.dat'),
                                       encoding='utf-8').read().splitlines() if l and not l.startswith('#')][1]
    for x, y in zip(star, exp['star']):
        if abs(float(x) - y) > 1e-8 * max(1.0, abs(y)):
            bad.append('fig_results_star.dat value %s is not the deposited %r' % (x, y))
    mk = io.open(os.path.join(JCP_PAPER, 'figs', 'fig_results_marks.tex'), encoding='utf-8').read()
    for name, want in (('FigResFRstar', exp['FRstar']), ('FigResFRtrue', exp['FRtrue'])):
        m = re.search(r'\\def\\%s\{([\d.]+)\}' % name, mk)
        if not m or abs(float(m.group(1)) - want) > 5e-7:
            bad.append('fig_results_marks.tex %s is not %r' % (name, want))
    raw = io.open(os.path.join(JCP_PAPER, 'figs', 'fig_results_native.tex'), encoding='utf-8').read()
    m = re.search(r'ymin=([\d.e+-]+),\s*ymax=([\d.e+-]+)', raw)
    ys = [v for s in src for v in s[1:5]] + list(exp['star'][1:])
    if not m or not (float(m.group(1)) < min(ys) and max(ys) < float(m.group(2))):
        bad.append('the y axis of fig:results does not contain every plotted value (%.3g to %.3g)' % (min(ys), max(ys)))
    xs = [float(x) for x in re.findall(r'at \(axis cs:([\d.]+),70\)', raw)]
    if len(xs) != 2 or not (min(xs) < exp['FRstar'] < max(xs)) or max(abs(x - exp['FRstar']) for x in xs) > 0.02:
        bad.append('the two labels at y=70 do not straddle FR*=%.4f (x = %r)' % (exp['FRstar'], xs))
    ok = _verdict(exp, 'fig:results: %d plotted rows (smaller leakage branch), the S* point, 2 marks (FR*=%.6f) and %d '
                       'caption numbers against the deposit' % (len(dat), exp['FRstar'], exp['n_caption']))
    _count(LAST_N[0] + 5 * len(dat) + len(star) + 2 + 2)
    return ok
BUILD_PRODUCTS = ('.aux', '.log', '.out', '.pdf', '.synctex.gz')


def _expect_looseness():
    """Data half of check_fixedP_looseness() (see there)."""
    import math
    g, rows = _grid_rows()
    by, cell = {}, {}
    for r in rows:
        if r['conv'] == 'CONVERGED' and r['relL1']:
            q = r['B2'] / r['relL1']
            by.setdefault((r['L'], round(r['FR'], 6)), []).append(q)
            cell.setdefault((r['L'], r['eta']), []).append(q)
    q = {}
    for (L, fr), v in by.items():
        if len(v) >= 2:
            q.setdefault(L, []).append(max(v) / min(v))
    m68 = max(max(q[6]), max(q[8]))
    mall = max(max(v) for v in q.values())
    cc = [max(v) / min(v) for v in cell.values()]
    up1 = lambda x: math.ceil(10 * x - 1e-9) / 10
    exp = dict(data_bad=[], want=[], gone=[], regex=[], near=[])
    ref = _so('v2.json')['D22']['looseness_second_order_branch']
    for k, v in (('fixedP_over_eta_max_L68', m68), ('fixedP_over_eta_max_all', mall)):
        if abs(ref[k] - v) > 1e-9 * v:
            exp['data_bad'].append('%s recomputed %.12g, second_order/v2.json %.12g' % (k, v, ref[k]))
    if any(abs(a - b) > 1e-9 * b for a, b in zip(ref['across_subspaces_of_cell'], (min(cc), max(cc)))):
        exp['data_bad'].append('across-subspace range does not agree with second_order/v2.json D22')
    fac, fac68 = '%.1f' % up1(mall), '%.1f' % up1(m68)
    if fac68 != fac:
        exp['data_bad'].append('the L=6,8 maximum (%s) differs from the maximum over every size (%s): the sentence must '
                               'print both' % (fac68, fac))
    lo_c, hi_c = _rnd(min(cc), 1, 'near', 'dp'), _rnd(max(cc), 1, 'near', 'dp')
    for rel in ('sec_3_body.tex', 'sm_frontier.tex'):
        exp['regex'].append((rel, r'by at most a factor \$([\d.]+)\$', fac,
                             'the fixed-subspace looseness of the second-order branch, "by at most a factor $X$"'))
        exp['near'].append((rel, 'by at most a factor $%s$' % fac, ['second-order'], 300, 0,
                            'the factor must name the second-order branch (D22, T6)'))
        exp['want'].append((rel, 'varies by a factor $%s$ to $%s$' % (lo_c, hi_c)))
        exp['gone'] += [(rel, 'by at most a factor $1.8$'), (rel, '$2.0$ up to $L=12$'), (rel, 'by a factor $23$ to $148$'),
                        (rel, '\\pm17\\%'), (rel, '\\pm25\\%')]
    exp['gone'].append(('sec_3_body.tex', 'changes little with the resolution'))
    exp['label'] = ('fixed-subspace looseness of the second-order branch over eta: max %.4f (L=6,8), %.4f (L=6..12), '
                    'printed rounded up as %s; across the subspaces of a cell %.4f to %.4f, printed %s to %s'
                    % (m68, mall, fac, min(cc), max(cc), lo_c, hi_c))
    return exp


def check_fixedP_looseness():
    """2026-11-08 (final check R1), rewritten 2026-11-16 for the second-order branch (rule D22; the number sheet rows
    T6, T7, F-40, F-41 and convention 1: the second-order branch, named).  On the CONVERGED rows of C3_grid.json the
    ratio B2 / true error of the second-order branch (B2 recomputed from w_S and Lam_over_eta):
      (a) at a fixed subspace (L, FR), its max/min over the eta's present; the maximum over L=6,8 and over L=6..12,
          printed rounded UP at one decimal ('by at most a factor $X$'; the two agree, so it is printed once), in a
          sentence that names the second-order branch (within 300 characters before it);
      (b) across the subspaces of one (L, eta) cell, its max/min per cell; the range over the sixteen cells printed to
          nearest at one decimal ('varies by a factor $a$ to $b$').
    Both in sec_3_body.tex and sm_frontier.tex.  The first-order values of the same two statements (1.8 / 2.0 and 23 to
    148, slack_leakbranch), 'changes little with the resolution' (D22) and the unsourced '+-17%' / '+-25%' clause must
    be gone."""
    exp = _expect_looseness()
    return _verdict(exp, exp['label'])

def check_drift_margin():
    """2026-11-08: Sec. S4 says all twenty published-subspace rows clear the reference-drift rule, the
    weakest (L=14, eta=0.10t, reference drift 1.9e-3 on halving its depth) by a factor 51.  Re-derived:
    margin = relL1 / ref_drift at the same eta (pt_L*_PUB.json meta.ref_drift; L14_FR008.json
    ref_drift, deepest rows).  The weakest row, its drift (2 s.f.) and its margin (nearest integer)
    must be the printed ones."""
    import glob
    import io
    import json
    pub = os.path.join(DATA, 'c3_frontier', 'published_fraction')
    m = []
    for f in sorted(glob.glob(os.path.join(pub, 'pt_L*_PUB.json'))):
        x = json.load(io.open(f, encoding='utf-8'))
        x = x[0] if isinstance(x, list) else x
        seen = set()
        for r in x['rows']:
            if 'relL1' in r and r['eta'] not in seen:
                seen.add(r['eta'])
                d = x['meta']['ref_drift']['%g' % r['eta']]
                m.append((r['relL1'] / d, x['meta']['L'], r['eta'], d))
    b = json.load(io.open(os.path.join(pub, 'L14_FR008.json'), encoding='utf-8'))[0]
    deep = max(r['n'] for r in b['rows'])
    for r in b['rows']:
        if r['n'] == deep and 'relL1' in r:
            d = b['ref_drift']['%g' % r['eta']]
            m.append((r['relL1'] / d, b['L'], r['eta'], d))
    q, L, eta, d = min(m)
    txt = ' '.join(io.open(os.path.join(JCP_PAPER, 'sm_frontier.tex'), encoding='utf-8').read().split())
    mant, expo = ('%.1e' % d).split('e')
    want = ('$L=%d$ at $\\eta=%.2f\\,t$, where the reference itself drifts by $%s\\times10^{%d}$ on '
            'halving its depth---by a factor $%d$' % (L, eta, mant, int(expo), round(q)))
    ok = len(m) == 20 and min(x[0] for x in m) > 20 and want in txt
    if not ok:
        print('  FAIL  Sec. S4 does not print %r (or a row fails the rule)' % want)
    print('  %s  reference-drift rule: 20 rows, weakest L=%d eta=%.2f, drift %.2e, margin %.1f'
          % ('ok  ' if ok else 'FAIL', L, eta, d, q))
    _count(len(m) + 1)
    return ok

def check_withdrawn_unsourced():
    """2026-11-08 (rule: no printed number without a source): the SI numbers withdrawn as unsourced must
    stay out of the printed text, and the two sourced statements that replaced them must hold:
    Lambda_K and the dense Lambda_S agree to within 1% where both were computed
    (data/c3_frontier/validation.json, ratio_K_S), and the orthogonality defect of the published-subspace
    rows is at most 4.9e-15 (published_fraction/pt_L*_PUB.json orth_defect_sampled, L14_FR008.json
    orth_defect)."""
    import glob
    import io
    import json
    ok = True
    v = json.load(io.open(os.path.join(DATA, 'c3_frontier', 'validation.json'), encoding='utf-8'))
    r = [x['ratio_K_S'] for blk in v.values() if isinstance(blk, dict) for x in blk.get('rows', [])
         if 'ratio_K_S' in x]
    if not r or max(abs(q - 1.0) for q in r) >= 0.01:
        print('  FAIL  validation.json ratio_K_S is not within 1%% of one: %r' % ((min(r), max(r)) if r else None))
        ok = False
    pub = os.path.join(DATA, 'c3_frontier', 'published_fraction')
    od = []
    for fn in glob.glob(os.path.join(pub, 'pt_L*_PUB.json')):
        x = json.load(io.open(fn, encoding='utf-8'))
        x = x[0] if isinstance(x, list) else x
        od.append(x['orth_defect_sampled'] if not isinstance(x['orth_defect_sampled'], dict)
                  else max(x['orth_defect_sampled'].values()))
    od.append(json.load(io.open(os.path.join(pub, 'L14_FR008.json'), encoding='utf-8'))[0]['orth_defect'])
    if max(od) > 4.95e-15:
        print('  FAIL  orthogonality defect of the published rows reaches %.3e > 4.9e-15' % max(od))
        ok = False
    gone = ('$1169$', '0.5258', '1.4868', '[2.34,5.61]', 'exhaustive search over', '1.000000000000',
            'descends to $1.61$', 'inflated thirtyfold')
    for rel in ('sm_leakage.tex', 'sm_frontier.tex', 'sm_nogo.tex', 'sec_6_body.tex'):
        lines = io.open(os.path.join(JCP_PAPER, rel), encoding='utf-8').read().splitlines()
        txt = ' '.join(l for l in lines if not l.lstrip().startswith('%'))
        for g in gone:
            if g in txt:
                print('  FAIL  %s prints the withdrawn %r' % (rel, g))
                ok = False
    for rel, w in (('sm_leakage.tex', 'they agree to within $1\\%$'),
                   ('sm_leakage.tex', 'is at most $4.9\\times10^{-15}$')):
        if w not in ' '.join(io.open(os.path.join(JCP_PAPER, rel), encoding='utf-8').read().split()):
            print('  FAIL  %s no longer prints %r' % (rel, w))
            ok = False
    print('  %s  withdrawn unsourced numbers absent; ratio_K_S %.4f..%.4f; orthogonality defect max %.2e'
          % ('ok  ' if ok else 'FAIL', min(r), max(r), max(od)))
    _count(len(r) + len(od) + 4 * len(gone) + 2)
    return ok

def _expect_table_III():
    """Data half of the Table III part of check_main_numbers() (see there)."""
    import io
    import json
    import math
    J = lambda *p: json.load(io.open(os.path.join(DATA, *p), encoding='utf-8'))
    pub = J('c3_frontier', 'published_fraction', 'PUB_rows.json')
    l14 = J('c3_frontier', 'published_fraction', 'L14_FR008.json')[0]
    deep = max(r['n'] for r in l14['rows'])
    exp = dict(data_bad=[], want=[], gone=[], regex=[], near=[])
    vac1, vac2, rel = {}, {}, {}
    for r in pub:
        w = r['w_S']
        x = r['B'] / (1.0 + math.sqrt(w)) - math.sqrt(1.0 - w)
        if abs(r['vac'] - r['B'] / r['triv']) > 1e-12 * r['vac']:
            exp['data_bad'].append('PUB_rows.json L=%d eta=%g: vac is not B/triv' % (r['L'], r['eta']))
        vac1[(r['L'], r['eta'])] = r['vac']
        vac2[(r['L'], r['eta'])] = _B2(w, x) / r['triv']
        if abs(r['eta'] - 0.15) < 1e-12:
            rel[r['L']] = r['relL1']
    for r in l14['rows']:
        if r['n'] == deep:
            if abs(_B1(l14['w_S'], r['LamHat_over_eta']) - r['bound_leak']) > 1e-9 * r['bound_leak']:
                exp['data_bad'].append('L14_FR008.json eta=%g: bound_leak is not B1(w_S, LamHat_over_eta)' % r['eta'])
            vac1[(14, r['eta'])] = r['vacuity_factor']
            vac2[(14, r['eta'])] = _B2(l14['w_S'], r['LamHat_over_eta']) / r['bound_trivial']
            if abs(r['eta'] - 0.15) < 1e-12:
                rel[14] = r['relL1']
    vac = {k: min(vac1[k], vac2[k]) for k in vac1}
    for t in _so('v2.json')['D17']['table_III']:
        k = (t['L'], t['eta'])
        if k not in vac or abs(vac[k] - t['printed_full']) > 1e-9 * t['printed_full']:
            exp['data_bad'].append('Table III L=%d eta=%g: smaller branch disagrees with second_order/v2.json D17' % k)
    if len(vac) != 20 or min(vac.values()) <= 1.0:
        exp['data_bad'].append('Table III: a published subspace is no longer vacuous (D13) or a cell is missing')
    s5 = 'sec_5_body.tex'
    for L in (6, 8, 10, 12, 14):
        exp['want'].append((s5, ' & '.join('$%.2f$' % vac[(L, e)] for e in (0.10, 0.15, 0.18, 0.25))
                            + ' & $%.4f$' % rel[L]))
    lo, hi = min(vac.values()), max(vac.values())
    e18 = [v for (L, e), v in vac.items() if abs(e - 0.18) < 1e-12]
    for rel_f in ('sec_1.tex', 'sec_10.tex'):
        exp['want'].append((rel_f, 'factors of $%.2f$ to $%.2f$' % (lo, hi)))
        exp['gone'].append((rel_f, 'factors of $1.66$ to $8.61$'))
    exp['want'] += [(s5, '$%.2f$ to $%.2f$ at the $\\eta=0.18\\,t$' % (min(e18), max(e18))),
                    (s5, 'smaller of the two leakage branches'), ('sec_1.tex', 'smaller of the two leakage branches')]
    exp['gone'] += [(s5, 'from $1.66$ ($L=6$')]
    exp['want'].append((s5, '$%.4f$--$%.4f$ against the threshold $0.05$' % (min(rel.values()), max(rel.values()))))
    nsec = sum(1 for k in vac if vac2[k] < vac1[k])
    exp['label'] = ('Table III: 20 cells, the smaller of the two leakage branches over the trivial one (second-order in '
                    '%d cell(s)), %.2f to %.2f; %.2f to %.2f at eta=0.18t' % (nsec, lo, hi, min(e18), max(e18)))
    return exp


def check_main_numbers():
    """2026-11-08 (systematic sweep, the number sweep of 2026-11-08): main-text numbers that no
    check of verify.py prints, re-derived here from the deposit and required verbatim in the source.
      Table III (sec_5_body.tex): the twenty vacuity factors and the five relative L1 errors at
        eta=0.15t, from published_fraction/PUB_rows.json and L14_FR008.json (deepest rows);
      Fig. 3 caption (carried/fig_method.tex): 0.59 -> 0.0185 (factor 32), 75 -> 242 determinants,
        2.4e-4 and 2.7%, from method_max.json;
      Sec. III D: 378 deposited subspaces, 357 with K_S=-1 (cert_stress/akw/teqsci.json);
      Sec. IV: 296 of 300 (z_suppK_sym_L6.json 'allowed'); Table II from moments_table.json.
    2026-11-16 (second-order branch; rule D7, D13, D17; the number sheet rows I1, V5, V6, V12, V14, V15, C5): each
    vacuity factor of Table III is the SMALLER of the two leakage branches over the trivial branch -- the first-order
    one is the stored vac (PUB_rows.json) / vacuity_factor (L14_FR008.json, deepest rows), as before; the second-order
    one is (sqrt(1-w)+x)^2/(1+w), x recovered from the stored first-order B (PUB_rows.json) or read as LamHat_over_eta
    (L14_FR008.json) -- and agrees with second_order/v2.json D17.table_III; all twenty stay above one (D13).  The range is
    printed in Secs. I and X, the eta=0.18t range in Sec. V, and Secs. I and V name the rule ('smaller of the two
    leakage branches'); the first-order range 1.66--8.61 must be gone."""
    import io
    import json
    bad = []
    J = lambda *p: json.load(io.open(os.path.join(DATA, *p), encoding='utf-8'))
    txt = lambda rel: ' '.join(io.open(os.path.join(JCP_PAPER, rel), encoding='utf-8').read().split())
    # Table III and the sentences that quote it are read as printed (comments removed, 2026-11-16): a commented-out
    # row satisfied the earlier raw-text test, and a commented-out superseded value must not fail it
    t3 = _expect_table_III()
    bad += t3['data_bad']
    for f, lit in t3['want']:
        if lit not in _printed(f):
            bad.append('%s does not print %r' % (f, lit))
    for f, lit in t3['gone']:
        if lit in _printed(f):
            bad.append('%s still prints the superseded %r' % (f, lit))
    m = J('method_max.json')
    s0 = [x for x in m['sampled'] if x['K'] == 0][0]
    s12 = [x for x in m['sampled'] if x['K'] == 12][0]
    fm = txt('carried/fig_method.tex')
    for lit in ('from $%.2f$ at $K=0$ to $%.4f$ at $K=12$' % (s0['l1'], s12['l1']),
                'by a factor $%d$' % round(s0['l1'] / s12['l1']),
                'grows from $%d$ to $%d$ determinants' % (round(s0['S']), round(s12['S'])),
                '$%.1f\\times10^{-4}$ relative' % (m['sum_rule']['reconstruction_deviation_relative'] * 1e4),
                '$%.1f\\%%$ short of $0.5$' % (100 * m['sum_rule']['window_loss_relative'])):
        if lit not in fm:
            bad.append('Fig. 3 caption does not print %r' % lit)
    n = k = 0
    for fn in ('cert_stress.json', 'cert_akw.json', 'cert_teqsci.json'):
        rows = [r for r in J(fn)['rows'] if 'K_S' in r]
        n += len(rows)
        k += sum(1 for r in rows if r['K_S'] == -1)
    if 'across the $%d$ deposited subspaces, of which $%d$ have $K_S=-1$' % (n, k) not in txt('sec_3_1.tex'):
        bad.append('Sec. III D counts are not %d / %d' % (n, k))
    a6 = J('c3_frontier', 'adversarial', 'z_suppK_sym_L6.json')
    if 'stalls at $%d$ of $%d$ determinants' % (a6['allowed'], a6['sector_dimension']) not in txt('sec_4_body.tex'):
        bad.append('Sec. IV does not print %d of %d' % (a6['allowed'], a6['sector_dimension']))
    mt = J('moments_table.json')['table']['rows']
    t2 = txt('sec_4_body.tex')
    want = ['$0\\ \\longrightarrow\\ %.1f\\times10^{-4}$' % (mt[0]['rel_eps_2K2'] * 1e4),
            '$%.1f\\times10^{-16}\\ \\longrightarrow\\ %.1f\\times10^{-8}$' % (mt[1]['max_rel_eps_upto_2K1'] * 1e16,
                                                                           mt[1]['rel_eps_2K2'] * 1e8),
            '$(%.2f\\times10^{-2})$' % (mt[0]['abs_eps_2K2'] * 1e2), '$(%.2f\\times10^{-4})$' % (mt[1]['abs_eps_2K2'] * 1e4)]
    if mt[0]['max_rel_eps_upto_2K1'] != 0.0 or mt[2]['max_rel_eps_upto_2K1'] != 0.0:
        bad.append('Table II prints 0 where the deposit is not zero')
    for w in want:
        if w not in t2:
            bad.append('Table II does not print %r' % w)
    for b in bad:
        print('  FAIL  %s' % b)
    print('  %s  main-text numbers: Table III (25 cells), vacuity range, Fig. 3 caption, 378/357, 296/300, Table II'
          % ('ok  ' if not bad else 'FAIL'))
    _count(len(t3['want']) + len(t3['gone']) + 5 + 1 + 1 + len(want) + 1)
    return not bad

SOURCES = os.path.join(HERE, 'jcp', 'number_sources.json')
SCRIPTS = os.path.join(HERE, 'jcp')                                     # fracgrid_numbers.py


def _printed(rel):
    """The printed text of one JCP source: comment lines and trailing comments removed, whitespace
    collapsed."""
    import io
    lines = io.open(os.path.join(JCP_PAPER, rel), encoding='utf-8').read().replace('\r\n', '\n').split('\n')
    return ' '.join(' '.join(re.sub(r'(?<!\\)%.*', '', l) for l in lines).split())


def _walk(obj, key):
    """Value at a key path written as /a/b[3]/c (the format of the builder of the registry (2026-11-08))."""
    for part in re.findall(r'/([^/\[]+)|\[(\d+)\]', key):
        obj = obj[part[0]] if part[0] else obj[int(part[1])]
    return obj


def _source_value(e, cache):
    """The stored value a record of number_sources.json points to (key path, or expression on the deposited files
    'file' -> d and 'files' -> D); files are resolved by _data_file(), so a re-scored output can be named by its deposit
    path SECOND_ORDER/<name>.  Raises KeyError / IndexError / TypeError / ValueError when it cannot be read."""
    import io
    import json
    for f in [e['file']] + list(e.get('files') or []):
        if f not in cache:
            fp = _data_file(f)
            cache[f] = (io.open(fp, encoding='utf-8').read() if fp.endswith('.txt')
                        else json.load(io.open(fp, encoding='utf-8')))
    if e.get('expr'):
        # a quantity derived from the stored values (a ratio, a fit); the expression is part of
        # this guard's input and reads only deposited files: d is 'file', D the list 'files'
        D = [cache[f2] for f2 in (e.get('files') or [])]
        return float(eval(e['expr'], {'d': cache[e['file']], 'D': D, 're': re}))
    return float(_walk(cache[e['file']], e['key']))


def _literal_matches(e, x):
    """True if the printed literal of record e is the stored value x at its printed number of significant digits,
    or, for a record with 'round' up / down, x rounded in that direction at the printed decimals."""
    s = e['lit'].replace('\\,', '').replace(',', '')
    nsig = max(1, len(s.replace('.', '').lstrip('0')))
    v = float(s) * 10 ** e['exp']
    targets = [v] + ([v / 100.0] if e['pct'] else [])
    same = lambda a, b: (a == b) if (a == 0 or b == 0) else ('%.*e' % (nsig - 1, abs(a))) == ('%.*e' % (nsig - 1, abs(b)))
    ndec = len(s.split('.')[1]) if '.' in s else 0
    # printed precision read as decimal places too ('0.010' is 0.0097 rounded to three decimals)
    dec_ok = (not e['exp'] and not e['pct'] and ndec > 0 and ('%.*f' % (ndec, abs(x))) == s)
    if e.get('round') in ('up', 'down'):
        # 2026-11-10 (7b audit): a number printed as a bound ("at most", "within", "<=" -> 'up'; "at
        # least", a floor -> 'down'; the ends of a confidence interval outward) must be the stored value
        # rounded in that direction at its printed decimals, in the printed unit (x 10^-exp, or % if
        # pct); 'sign': -1 marks a literal printed after a minus sign (the lower end of an interval).
        import math
        X = (100.0 * x if e['pct'] else x) / 10 ** e['exp']
        q = 10 ** ndec
        want = (math.ceil(X * q - 1e-9) if e['round'] == 'up' else math.floor(X * q + 1e-9)) / q
        dec_ok = abs(e.get('sign', 1) * float(s) - want) < 0.5 / q
        targets = []
    return dec_ok or any(same(x, t) for t in targets)


def check_number_sources():
    """2026-11-08 (pass 7, systematic sweep): every number listed in src/jcp/number_sources.json is
    re-read from the deposited file and key path recorded there, scaled by the power of ten or the
    percent sign printed with it, and must equal the printed literal at its printed number of
    significant digits; the literal must still be printed in that source file.
    2026-11-16 (second-order branch, the number sheet rows GD1, GD13, GD14): the records of the numbers the change
    re-scores point to the re-scored outputs by their deposit path SECOND_ORDER/<name> (data/c3_frontier/
    second_order/, deposited since 3.5.0; see _data_file);
    reading and comparing are factored into _source_value() and _literal_matches(), unchanged."""
    import io
    import json
    if not os.path.isfile(SOURCES):
        print('  FAIL  %s is missing' % SOURCES)
        return False
    entries = json.load(io.open(SOURCES, encoding='utf-8'))
    cache, bad = {}, []
    for e in entries:
        try:
            x = _source_value(e, cache)
        except (KeyError, IndexError, TypeError, ValueError):
            bad.append('%s %s: key %s not found in %s' % (e['tex'], e['lit'], e['key'], e['file']))
            continue
        if not _literal_matches(e, x):
            bad.append('%s prints %s but %s %s = %r%s' % (e['tex'], e['lit'], e['file'], e['key'], x,
                                                         (' (rounded %s)' % e['round']) if e.get('round') else ''))
        if not re.search(r'(?<![\d.])' + re.escape(e['lit']) + r'(?![\d])', _printed(e['tex'])):
            bad.append('%s no longer prints %s' % (e['tex'], e['lit']))
    for b in bad:
        print('  FAIL  %s' % b)
    print('  %s  number sources: %d printed numbers re-read from %d deposited files'
          % ('ok  ' if not bad else 'FAIL', len(entries), len(cache)))
    _count(2 * len(entries))
    return not bad


def check_gapscaling_caption():
    """2026-11-08 (pass 7): the caption of the gap-scaling figure is generated by the deposited
    src/recovered/make_fig_gapscaling.py from data/gap_scaling.json, spinqw_L12.json and
    akw_lanczos_L12.json, and was then edited by hand (wording only).  The generator is run on a
    scratch copy of those four files, and every number printed in the JCP caption must appear in the
    generated caption (at the JCP rounding)."""
    import io
    tmp = tempfile.mkdtemp(prefix='p2jcp_gap_')
    try:
        for sub in ('src/recovered', 'data', 'paper/figs', 'build'):
            os.makedirs(os.path.join(tmp, *sub.split('/')))
        shutil.copy2(os.path.join(RELEASE, 'src', 'recovered', 'make_fig_gapscaling.py'),
                     os.path.join(tmp, 'src', 'recovered'))
        for f in ('gap_scaling.json', 'spinqw_L12.json', 'akw_lanczos_L12.json'):
            shutil.copy2(os.path.join(DATA, f), os.path.join(tmp, 'data'))
        p = subprocess.run([sys.executable, '-I', '-B', os.path.join(tmp, 'src', 'recovered', 'make_fig_gapscaling.py')],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=tmp)
        gen = io.open(os.path.join(tmp, 'build', 'fig_gapscaling_caption.tex'), encoding='utf-8').read() \
            if p.returncode == 0 else ''
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    if not gen:
        print('  FAIL  the deposited gap-scaling generator did not run')
        return False
    num = lambda t: re.findall(r'(?<![A-Za-z\d.])\d+(?:\.\d+)?', ' '.join(re.sub(r'(?<!\\)%.*', '', l)
                                                                     for l in t.split('\n')))
    have = set()
    for g in num(gen):
        have.add(g)
        if '.' in g:
            for d in range(0, len(g.split('.')[1])):
                have.add(('%.*f' % (d, float(g))))
    mine = num(_printed('figs/fig_gapscaling_caption.tex'))
    # 2026-11-10 (7b audit): 'agree to within $8.4\%$' is the bound rounded UP (max |1 - ratio| = 8.32 %);
    # the generator prints the nearest rounding.  The literal is re-derived by check_directed_rounding().
    mine = [x for x in mine if x not in DIRECTED_GAP]
    miss = [x for x in mine if x not in have]
    if miss:
        print('  FAIL  gap-scaling caption numbers not produced by the deposited generator: %r' % miss)
    print('  %s  gap-scaling caption: %d printed numbers, all produced by src/recovered/make_fig_gapscaling.py'
          % ('ok  ' if not miss else 'FAIL', len(mine)))
    _count(len(mine))
    return not miss


def check_fracgrid():
    """2026-11-08 (pass 7): Table tab:fracgrid and the sentences that quote it were rewritten from
    data/eta_sweep.json by src/jcp/fracgrid_numbers.py; every rewritten passage is recomputed
    from that file here and must be printed verbatim."""
    sys.path.insert(0, SCRIPTS)
    try:
        import fracgrid_numbers as rf
    finally:
        sys.path.pop(0)
    n = rf.numbers()
    bad = []
    for fn, old, new in rf.edits(n):
        if ' '.join(new.split()) not in _printed(fn):
            bad.append('%s does not print the recomputed %r' % (fn, ' '.join(new.split())[:70]))
    for b in bad:
        print('  FAIL  %s' % b)
    print('  %s  fraction grid: table, fits and %d quoting passages recomputed from eta_sweep.json'
          % ('ok  ' if not bad else 'FAIL', len(rf.edits(n))))
    _count(len(rf.edits(n)))
    return not bad


def check_frontier_anatomy():
    """2026-11-08 (pass 7): the rung-by-rung displacement of the threshold (Sec. V F, SM) and the
    Cauchy-Schwarz margin, recomputed from data/c3_frontier/verdict_frontier.json and ladder_L8.json.
    Shares: median over the uncensored (eta, target) cells of (FR*_k - FR*_{k-1}) / (FR*_5 - FR*_0),
    k = 1 pointwise, 5 operator norm, 4 Feshbach/Minkowski, 2 Cauchy-Schwarz; L=8, and L=6 with L=8
    pooled.  Margin: R2/R1 over the FR rows of ladder_L8.json, median and maximum."""
    import io
    import json
    import math
    import statistics as st
    vf = json.load(io.open(os.path.join(DATA, 'c3_frontier', 'verdict_frontier.json'), encoding='utf-8'))

    def shares(Ls):
        rows = [r for r in vf if r['L'] in Ls and all(math.isfinite(r['FRstar%d' % k]) for k in range(6))]
        return len(rows), {k: st.median([(r['FRstar%d' % k] - r['FRstar%d' % (k - 1)]) / r['d_total'] for r in rows])
                           for k in range(1, 6)}
    n8, s8 = shares((8,))
    n68, s68 = shares((6, 8))
    lad = json.load(io.open(os.path.join(DATA, 'c3_frontier', 'ladder_L8.json'), encoding='utf-8'))
    q = [r['R2'] / r['R1'] for c in lad['cases'] if c['tag'].startswith('FR=') for r in c['rows'] if r['R1'] > 1e-13]
    pct = lambda x: ('%.1f' % (100 * x)) if 100 * x < 10 else ('%d' % round(100 * x))
    a8 = [pct(s8[k]) for k in (1, 5, 4, 2)]
    a68 = [pct(s68[k]) for k in (1, 5, 4, 2)]
    want = [
        ('sec_5_body.tex', 'the displacement of the threshold is $%s\\%%$ the pointwise step' % a8[0]),
        ('sec_5_body.tex', '$%s\\%%$ the operator-norm step' % a8[1]),
        ('sec_5_body.tex', '$%s\\%%$ the Feshbach/Minkowski split and $%s\\%%$ Cauchy--Schwarz' % (a8[2], a8[3])),
        ('sec_5_body.tex', 'a factor $%.2f$ at the median and $%.2f$ at worst' % (st.median(q), max(q))),
        ('sm_frontier.tex', 'accounts for $%s\\%%$ of it' % a8[0]),
        ('sm_frontier.tex', 'for $%s\\%%$, the Feshbach/Minkowski split for $%s\\%%$, and the Cauchy--Schwarz step in'
         ' $\\omega$ for $%s\\%%$ (medians over $%d$ of $18$ cells' % (a8[1], a8[2], a8[3], n8)),
        ('sm_frontier.tex', 'a factor $%.2f$ at the median and $%.2f$ at the worst point measured, over the $%d$'
         ' ladder rows at $L=8$' % (st.median(q), max(q), len(q))),
        ('sec_5_app.tex', 'medians over $%d$ of the $18$' % n8),
        ('sec_5_app.tex', 'shifts them to $%s\\%%$, $%s\\%%$, $%s\\%%$ and $%s\\%%$' % tuple(a68)),
    ]
    bad = [(f, w) for f, w in want if w not in _printed(f)]
    # 2026-11-16 (second-order branch; rule D10; the number sheet rows V33, H3, F-54b; the check of the number sheet
    # discrepancy 2): these shares and margins are those of the chain of inequalities of the FIRST-ORDER branch; they are
    # kept, and each file that prints them must name the first-order chain within 400 characters before its first
    # share (and Sec. V within 150 characters of 'rung by rung'); Sec. V or Sec. S4 must say that the second-order proof
    # has no pointwise step.
    for f in ('sec_5_body.tex', 'sm_frontier.tex', 'sec_5_app.tex'):
        t = _printed(f)
        pos = [t.find(w) for g_, w in want if g_ == f and t.find(w) >= 0]
        if pos and 'first-order' not in t[max(0, min(pos) - 400):min(pos) + 200].lower():
            bad.append((f, 'the label "first-order" within 400 characters before its first share (D10)'))
    t5 = _printed('sec_5_body.tex')
    i = t5.find('rung by rung')
    if i >= 0 and 'first-order' not in t5[max(0, i - 150):i + 150].lower():
        bad.append(('sec_5_body.tex', 'the label "first-order" within 150 characters of "rung by rung" (D10)'))
    if not any('no pointwise step' in _printed(f) for f in ('sec_5_body.tex', 'sm_frontier.tex')):
        bad.append(('sec_5_body.tex', '"no pointwise step" (D10: the second-order proof has none), here or in sm_frontier.tex'))
    for f, w in bad:
        print('  FAIL  %s does not print %r' % (f, w))
    print('  %s  threshold anatomy (first-order chain, labelled): shares %s%% (L=8, %d cells), %s%% (L=6+8); Cauchy-Schwarz margin %.2f / %.2f'
          % ('ok  ' if not bad else 'FAIL', '/'.join(a8), n8, '/'.join(a68), st.median(q), max(q)))
    _count(len(want) + 3 + 1 + 1)
    return not bad


def _expect_threshold():
    """Data half of check_threshold_error() (see there)."""
    import io
    import json
    import math
    import statistics
    cells, bad = _second_order_cells()
    exp = dict(data_bad=list(bad), want=[], gone=[], regex=[], near=[])
    err = [c['err'] for c in cells]
    if len(err) != 16 or any(e is None for e in err):
        exp['data_bad'].append('C3_grid.json: %d cells with a crossing error' % len(err))
        return exp
    lo, hi, med = min(err), max(err), statistics.median(err)
    sp, mid, half = hi / lo, math.sqrt(lo * hi), math.sqrt(hi / lo)
    by = {(c['L'], c['eta']): c['err'] for c in cells}
    E = (0.1, 0.15, 0.18, 0.25)
    if not all(by[(a, x)] > by[(b, x)] for x in E for a, b in ((6, 8), (8, 10), (10, 12))) or \
            not all(by[(L, a)] < by[(L, b)] for L in (6, 8, 10, 12) for a, b in ((0.1, 0.15), (0.15, 0.18), (0.18, 0.25))):
        exp['data_bad'].append('the crossing error no longer falls with L at every eta and rises with eta at every L')
    if not all(c['FRstar'] > c['FR05'] for c in cells):
        exp['data_bad'].append('at the fraction where the true error first falls below 0.05 the certificate is no longer '
                               'trivial in every cell (D13)')
    if half > 3.0:
        exp['data_bad'].append('half-spread %.3f > 3: "within a factor of three" fails' % half)
    ref = _so('cells.json')['crossing']['new']
    for k, v in (('min', lo), ('max', hi), ('median', med), ('spread', sp)):
        if abs(ref[k] - v) > 1e-9 * v:
            exp['data_bad'].append('crossing %s recomputed %.12g, second_order/cells.json %.12g' % (k, v, ref[k]))
    margin = min(c['margin'] for c in cells)
    e18 = [by[(L, 0.18)] for L in (6, 8, 10, 12)]
    f612 = [by[(6, x)] / by[(12, x)] for x in E]
    fe = [by[(L, 0.25)] / by[(L, 0.1)] for L in (6, 8, 10, 12)]
    ratio = [c['ratio'] for c in cells]
    least = [c['least'] for c in cells]
    omw = sorted(c['omw'] for c in cells)
    es = sorted(err)
    fr8 = [c['FR05'] for c in cells if c['L'] == 8 and abs(c['eta'] - 0.18) < 1e-12][0]
    lo1, hi1 = min(c['err_old'] for c in cells), max(c['err_old'] for c in cells)
    L_, H_, M_ = _tsci(lo, 2), _tsci(hi, 2), _tsci(med, 2)
    S_ = _rnd(sp, 2)
    R_ = '%s to %s' % (_tsci(min(ratio), 2), _tsci(max(ratio), 2))
    Q_ = '%s to %s' % (_tsci(min(least), 2), _tsci(max(least), 2))
    F1 = '$%.1f\\times10^{-3}$ to $%.1f\\times10^{-3}$' % (lo1 * 1e3, hi1 * 1e3)
    exp['want'] += [
        ('main.tex', '%s--%s over sixteen $(L,\\eta)$ cells' % (L_, H_)),
        ('sec_1.tex', '%s and %s: a spread of $%s$' % (L_, H_, S_)),
        ('sec_5_body.tex', '%s to %s in every one of sixteen' % (L_, H_)),
        ('sec_5_body.tex', '%s to %s, median %s: a spread of $%s$' % (L_, H_, M_, S_)),
        ('sm_frontier.tex', '%s to %s, median %s: a spread of $%s$' % (L_, H_, M_, S_)),
        ('sec_10.tex', 'between %s and %s' % (L_, H_)),
        ('sec_5_body.tex', 'at least $%d$ times the drift of our own reference' % math.floor(margin)),
        # 2026-11-18 (round 2, defect 6): the floor reads "at least" here too (correction round 2 R2-22)
        ('sec_5_app.tex', 'standing at a factor of at least $%d$' % math.floor(margin)),
        ('sm_frontier.tex', 'at least $%d$ times the drift of our own' % math.floor(margin)),
        ('sm_frontier.tex', '$%s$, $%s$, $%s$ and $%s\\times10^{-3}$ at $\\eta=0.18\\,t$ for $L=6$--$12$, a factor '
         '$%.1f$--$%.1f$' % tuple([_rnd(e * 1e3, 3) for e in e18] + [min(f612), max(f612)])),
        ('sm_frontier.tex', '$%.1f$--$%.1f$, the product of the two being the observed $%s$' % (min(fe), max(fe), S_)),
        ('sm_frontier.tex', 'over a factor $%s$ against the true error\'s $%s$, and still over $%s$ against $%s$'
         % (_rnd(omw[-1] / omw[0], 2), S_, _rnd(omw[-2] / omw[1], 2), _rnd(es[-2] / es[1], 2))),
        ('sec_5_body.tex', '$\\mathrm{FR}\\approx%.3f$, the certificate requires' % fr8),
        ('sm_frontier.tex', 'is reached at $\\mathrm{FR}\\approx%.3f$, while' % fr8),
        ('sec_9.tex', '${\\approx}%d\\times10^{-3}$ to within a factor of three' % round(mid * 1e3)),
        ('sm_frontier.tex', '${\\approx}%d\\times10^{-3}$' % round(mid * 1e3)),
        ('sec_1.tex', R_), ('sec_3_body.tex', R_), ('sm_frontier.tex', R_),
        ('sec_3_body.tex', Q_), ('sm_frontier.tex', Q_)]
    # D11: the first-order crossing error is kept ONCE in the body, in the section on the crossing, labelled
    exp['near'].append(('sec_5_body.tex', F1, ['first-order'], 200, 0,
                        'the first-order crossing error (D11) must be labelled as the first-order branch alone'))
    exp['count'] = [('sec_5_body.tex', F1, 1)]
    exp['gone'] += [
        ('main.tex', '$%.1f$--$%.1f\\times10^{-3}$' % (lo1 * 1e3, hi1 * 1e3)),
        ('sec_1.tex', '$%.1f\\times10^{-3}$ and $%.1f\\times10^{-3}$' % (lo1 * 1e3, hi1 * 1e3)),
        ('sec_10.tex', '$%.1f\\times10^{-3}$ and $%.1f\\times10^{-3}$' % (lo1 * 1e3, hi1 * 1e3)),
        ('sm_frontier.tex', F1), ('sec_5_body.tex', F1 + ', median'), ('sec_5_body.tex', F1 + ' in every one'),
        ('sec_1.tex', 'a spread of $5.9$'), ('sec_5_body.tex', 'a spread of $5.9$'), ('sm_frontier.tex', 'a spread of $5.9$'),
        ('sec_5_body.tex', '$46$ times the drift'), ('sm_frontier.tex', '$46$ times the drift'),
        ('sec_5_app.tex', 'a factor $46$'), ('sec_9.tex', '${\\approx}3\\times10^{-3}$'),
        ('sm_frontier.tex', '${\\approx}3\\times10^{-3}$'),
        ('sec_1.tex', '$3\\times10^{2}$ to $1.6\\times10^{3}$'), ('sec_3_body.tex', '$3\\times10^{2}$ to $1.6\\times10^{3}$'),
        ('sm_frontier.tex', '$3\\times10^{2}$ to $1.6\\times10^{3}$'),
        ('sec_3_body.tex', '$3\\times10^{2}$ to $2.4\\times10^{3}$'), ('sm_frontier.tex', '$3\\times10^{2}$ to $2.4\\times10^{3}$'),
        ('sm_frontier.tex', '$5.65$, $3.50$, $2.52$')]
    exp['label'] = ('true error at the certificate\'s crossing: %d cells, %.3e to %.3e (median %.3e, spread %.3f, middle '
                    '%.1e x/: %.2f, at least %.1f x the drift); bound/error there %.0f..%.0f, least loose %.0f..%.0f; '
                    'first-order %.2e..%.2e kept once (D11)'
                    % (len(err), lo, hi, med, sp, mid, half, margin, min(ratio), max(ratio), min(least), max(least),
                       lo1, hi1))
    return exp


def check_threshold_error():
    """2026-11-08 (pass 7): the true relative L1 error at the non-vacuity threshold, for each of the
    sixteen (L, eta) cells of C3_grid.json 'frontier': relL1 of the grid rows of that cell interpolated
    log-linearly in the fraction at FR*.
    2026-11-16, rewritten for the second-order branch (rule D6, D11, D13; the number sheet rows M1, I3, I4, V2,
    V25, V25b, V26, V27, H1, T8, T9, C1, N1, F-34 to F-39, F-42, F-48, F-52; convention 4 of the sheet: a range of
    stored values rounds to nearest): FR* is now where the CERTIFICATE (min{1+w, B2}, recomputed from w_S and
    Lam_over_eta) first improves on the trivial bound, by the fsum.py convention that reproduces the stored FR_star
    (_second_order_cells(), second_order/cells.json).  Printed: the range of the error there at two figures, its median and
    spread (abstract, Secs. I, V, X, S4); the margin over the reference drift, rounded DOWN (V, S15, S4); the four
    eta=0.18t values at three figures and the L=6/L=12 factor (S4); the eta factor and the 1-w comparison (S4); the
    geometric middle (IX, S4); bound/error at the crossing and at the least loose certified CONVERGED point (I, III,
    S4).  The first-order range is printed ONCE, in Sec. V, labelled (D11); its other printings, and every superseded
    first-order value of these sentences, must be gone.  Data: the crossing error falls with L and rises with eta in
    every cell, and the certificate is trivial where the true error first falls below 0.05 (D13), in all sixteen."""
    exp = _expect_threshold()
    for f, lit, n in exp.get('count', []):
        k = _printed(f).count(lit)
        if k != n:
            exp['data_bad'].append('%s prints %r %d times; D11 wants it exactly %d time(s)' % (f, lit, k, n))
    return _verdict(exp, exp.get('label', 'true error at the certificate\'s crossing'))


# A main-text or SM number that restates one which verify.py checks in another source file (shown by
# the mutation probe of 2026-11-08): the restatement must print the same literal as the anchor.
RESTATED = [
    ('sec_3_1.tex', '606', 'figs/caption_thm1iii_violation.tex', 'violated on none of the $606$'),
    ('sec_5_body.tex', '606', 'figs/caption_thm1iii_violation.tex', 'violated on none of the $606$'),
    ('sm_frontier.tex', '606', 'figs/caption_thm1iii_violation.tex', 'violated on none of the $606$'),
    # 2026-11-16 (second-order branch, rule D18; the number sheet row GD8): the two restatements of the pooled count
    # 61 (first-order branch) are replaced by check_pooled_split(), which re-derives 85 = 57 + 28 from the stored pool
    # and requires it, with its split, in all four places (the caption, Secs. S9, V and S4).
]


def check_printed_arithmetic():
    """2026-11-08 (pass 8): numbers that are arithmetic on other printed numbers of the same sentence --
    the two-sided p-values of the paired t statistics of App. dequant (9 dof), the difference of the two
    companion-study values, and the exact floors 2/2^5 and 2/2^10 of the sign and signed-rank tests."""
    from scipy import stats
    txt = _printed('app_dequant.tex')
    bad = []
    for t9, p in (('1.19', 0.26), ('3.27', 0.010)):
        got = 2 * stats.t.sf(float(t9), 9)
        if not (re.search(re.escape(t9), txt) and abs(round(got, 3 if p < 0.1 else 2) - p) < 1e-12):
            bad.append('t_9=%s gives p=%.4f, printed %s' % (t9, got, p))
    for w in ('($p=0.26$)', '$p=0.010$', '$p=2/2^{5}=0.0625$', 'floor is $0.002$',
              '$-1.62\pm3.46$', '$-1.83$', 'second to $0.21$'):
        if w not in txt:
            bad.append('App. dequant no longer prints %r' % w)
    if abs(1.83 - 1.62 - 0.21) > 1e-9 or abs(2 / 2 ** 5 - 0.0625) > 1e-12 or round(2 / 2 ** 10, 3) != 0.002:
        bad.append('arithmetic')
    for b in bad:
        print('  FAIL  %s' % b)
    print('  %s  printed arithmetic of App. dequant: two p-values, one difference, two test floors'
          % ('ok  ' if not bad else 'FAIL'))
    _count(2 + 7 + 1)
    return not bad


DIRECTED_GAP = {'8.4', '0.9168'}   # gap-scaling caption literals the generator does not print at this precision:
# the bound 8.4 rounded UP and (2026-11-14, closing reading A2) the ratio that sets it, printed with four decimals
# so that 1 - 0.9168 = 8.32 % visibly lies within 8.4 %; both are re-derived in check_directed_rounding()


def check_directed_rounding():
    """2026-11-10 (7b audit, the rounding audit of 2026-11-10): a number printed as a bound is the stored
    value rounded in the direction of the bound -- "at most", "up to", "within", "<=" round UP; "at least" and
    a floor round DOWN -- and not to nearest.  The entries that src/jcp/number_sources.json cannot express
    (a percentage computed from a formula, a floor that is a product, a spread, a word) are re-derived here
    from the deposit; the old nearest-rounded literals must be gone."""
    import glob
    import io
    import json
    import math
    import statistics
    J = lambda *p: json.load(io.open(os.path.join(DATA, *p), encoding='utf-8'))
    up = lambda x, n: math.ceil(x * 10 ** n - 1e-9) / 10 ** n
    dn = lambda x, n: math.floor(x * 10 ** n + 1e-9) / 10 ** n
    want, got = [], {}
    # (1) the assumed floor w_S >= 0.995: tau >= 0.9280, so the threshold moves by at most 1 - tau(0.995)
    #     2026-11-14 (closing reading A2): printed at one decimal, rounded UP (7.196 % -> 7.2 %); the integer
    #     ceiling 8 % could not be recomputed by a reader from the printed 0.9280
    #     2026-11-16 (second-order branch; rule D6, D12; the number sheet rows V28, V29, F-44, F-45, F-46): the
    #     certificate's threshold for a real probe is tau(w) = sqrt(1+w) - sqrt(1-w), tau(1) = sqrt 2, so the floor gives
    #     tau >= tau(0.995) = 1.3417 (printed 1.341, rounded DOWN, the same literal in all three places) and moves the
    #     threshold by at most 1 - tau(0.995)/tau(1) = 5.125 % (printed 5.2 %, rounded UP); 0.9280 / 0.928 / 7.2 % go
    w = 0.995
    tau = _tau2(w)
    mv = up(100 * (1 - tau / _tau2(1.0)), 1)
    fl = dn(tau, 3)
    got['threshold move'] = '%.3f%% -> at most %.1f%%; tau(0.995) %.5f -> %.3f' % (100 * (1 - tau / _tau2(1.0)), mv,
                                                                                  tau, fl)
    want += [('sec_5_body.tex', 'moving the threshold by at most $%.1f\\%%$' % mv),
             ('sm_frontier.tex', 'at the cost of at most $%.1f\\%%$ in the threshold' % mv),
             ('sec_5_body.tex', '$\\tau\\ge%.3f$' % fl), ('sm_frontier.tex', '$\\tau\\ge%.3f$' % fl),
             ('sm_frontier.tex', '<%.3f$' % fl)]
    # (2) the L=14 shot budget: the L=12 overhead T/|S| (run B) times |S|_eta(14).  2026-11-14 (closing reading A1):
    #     an EXTRAPOLATION, not a floor -- the overhead exponent is not resolved (honest_sampling.json
    #     fits/honest_overhead_TperS_vs_L: ci95 contains zero, significant_positive_95 false), so the L=12 overhead
    #     is not a lower bound on the L=14 one.  A plain value rounds to NEAREST (2 s.f.: 1.268e8 -> 1.3e8), the
    #     floor wording must be gone, and the comparison with the largest budget drawn (16 x T of the Fig. 3 run,
    #     8.3e7) must still hold.
    hs = J('honest_sampling.json')
    pub = J('c3_frontier', 'published_fraction', 'PUB_rows.json')
    S12 = [r['S'] for r in pub if r['L'] == 12 and abs(r['eta'] - 0.15) < 1e-12][0]
    S14 = round(J('eta_sweep.json')['published_sequence_L4_to_L14']['values'][5] * 10306296)
    fl = hs['per_L']['12']['honest_T_mean'] / S12 * S14
    m = round(fl / 1e8, 1)
    if round(154 * S14 / 1e8, 1) != m:
        want.append(('sec_7.tex', 'the printed overhead 154 gives another extrapolation'))
    if hs['fits']['honest_overhead_TperS_vs_L']['significant_positive_95']:
        want.append(('sec_7.tex', 'the overhead fit now resolves a positive exponent: re-judge the wording'))
    big = 16 * J('fig3_finite_shot_T5.2e6.json')['config']['T']
    if not fl > big:
        want.append(('sec_7.tex', 'the extrapolation no longer exceeds the largest budget drawn (%.3g)' % big))
    got['L=14 extrapolation'] = '%.4e -> %.1f x 1e8 (> %.2e drawn)' % (fl, m, big)
    # 2026-11-14 (closing round 2, D7): the overhead HELD FIXED is the one against the PUBLISHED support (153.8 at
    #     L=12), so the fit that decides "extrapolated, not a floor" is its own: OLS of log(T/|S|_published) on
    #     L=6..12, t interval with n-2 = 2 degrees of freedom (printed in Sec. 7.3 as 0.08, [-0.03, 0.19]).  If its
    #     95% interval excludes zero, the wording must be re-judged.  The stored fit read above is the one against
    #     the populated support; both are now checked.
    from scipy import stats
    Ls = [6, 8, 10, 12]
    Sp = {L: [r['S'] for r in pub if r['L'] == L and abs(r['eta'] - 0.15) < 1e-12][0] for L in Ls}
    ofit = stats.linregress(Ls, [math.log(hs['per_L'][str(L)]['honest_T_mean'] / Sp[L]) for L in Ls])
    tq = stats.t.ppf(0.975, len(Ls) - 2)
    olo, ohi = ofit.slope - tq * ofit.stderr, ofit.slope + tq * ofit.stderr
    if not olo < 0 < ohi:
        want.append(('sec_7.tex', 'the published-support overhead fit excludes zero, [%.4f, %.4f]: re-judge the wording'
                     % (olo, ohi)))
    got['overhead fit, published support'] = '%.4f [%.4f, %.4f]' % (ofit.slope, olo, ohi)
    # 2026-11-14 (closing round 2, D6): the deposit's own L=14 value is that of the log-linear T-vs-L fit
    #     (extrapolation_L14/honest_T14_from_fit, re-derived here from fits/honest_T_vs_L); it is printed beside the
    #     overhead extrapolation as a plain value (2 s.f., nearest; since round 3 once, with its interval, see
    #     below) and must exceed the budget drawn too.
    ft = hs['fits']['honest_T_vs_L']
    t14 = math.exp(ft['intercept'] + 14 * ft['slope'])
    if abs(t14 / hs['extrapolation_L14']['honest_T14_from_fit'] - 1) > 1e-9:
        want.append(('sec_7.tex', 'the stored L=14 fit value is not the T-vs-L fit at L=14'))
    if not t14 > big:
        want.append(('sec_7.tex', 'the fitted L=14 value no longer exceeds the largest budget drawn (%.3g)' % big))
    f14 = round(t14 / 1e8, 1)
    got['L=14 from the T-vs-L fit'] = '%.4e -> %.1f x 1e8' % (t14, f14)
    # 2026-11-14 (closing round 3, D2): the fitted value is printed ONCE (S7 opening), with the 95% confidence
    #     interval of the fitted value at L=14: t(0.975, n-2) x se x sqrt(Sxx/n + (14 - mean L)^2) on log T, from the
    #     stored fit (n = 4, L = 6..12; that its stored slope interval is t(0.975, 2) x se is checked here).  The ends
    #     round OUTWARD (lower DOWN, upper UP, two figures).  The lower end lies below the largest budget drawn, so
    #     the comparison of Sec. 7.3 names only the extrapolation at the L=12 overhead; the round-2 clause that set
    #     the fit beside it is forbidden below.  If the lower end ever exceeds the budget drawn, re-judge the wording.
    Lf = [6, 8, 10, 12]
    tq2 = stats.t.ppf(0.975, len(Lf) - 2)
    if ft['n'] != len(Lf) or abs(ft['slope'] - tq2 * ft['se'] - ft['ci95'][0]) > 1e-9:
        want.append(('sec_7.tex', 'the stored T-vs-L fit is not the n=4, t(0.975,2) fit the interval assumes'))
    xb = sum(Lf) / len(Lf)
    hw = tq2 * ft['se'] * math.sqrt(sum((L - xb) ** 2 for L in Lf) / len(Lf) + (14 - xb) ** 2)
    clo = math.exp(ft['intercept'] + 14 * ft['slope'] - hw)
    chi = math.exp(ft['intercept'] + 14 * ft['slope'] + hw)
    if not clo < big:
        want.append(('sec_7.tex', 'the lower end of the fit interval (%.3g) no longer lies below the budget drawn: '
                     're-judge the wording of Sec. 7.3' % clo))
    c_lo, c_hi = dn(clo / 1e7, 1), up(chi / 1e8, 1)
    got['L=14 fit, 95% CI of the fitted value'] = '[%.4e, %.4e] -> [%.1fe7, %.1fe8]' % (clo, chi, c_lo, c_hi)
    want += [('sec_7.tex', 'for a single seed and branch (the log-linear fit itself gives ${\\approx}%.1f\\times10^{8}$, '
              '$95\\%%$ confidence interval $[%.1f\\times10^{7},%.1f\\times10^{8}]$).' % (f14, c_lo, c_hi)),
             ('sec_7.tex', '(Table~\\ref{tab:shots}), exceeds the largest shot budget drawn anywhere in this work, '
              'the $8.3\\times10^{7}$ shots')]
    want += [('sec_7.tex', 'its $L=14$ value, extrapolated at the $L=12$ overhead, being ${\\approx}%.1f\\times10^{8}$ shots'
              % m),
             ('sec_7.tex', 'Only $L=14$ is extrapolated, not measured: the overhead is held at its $L=12$ value.'),
             ('sec_7.tex', '${\\approx}%.1f\\times10^{8}$ & $154$ & extrapolated, not measured' % m),
             ('sec_7.tex', 'The $L=14$ extrapolation at the $L=12$ overhead, ${\\approx}%.1f\\times10^{8}$ shots for a' % m)]
    # (3) the true error at the crossing: "track each other to within a factor X", X = spread rounded UP
    #     2026-11-16 (second-order branch): the crossing is the certificate's (_second_order_cells()); its spread is
    #     printed once, by check_threshold_error()
    err = [c['err'] for c in _second_order_cells()[0]]
    sp = max(err) / min(err)
    # 2026-11-14 (closing reading A2): the sentence no longer prints the spread a second time (5.9 and, rounded up,
    #     6.0 were one stored value printed two ways); it refers to the spread printed above it, which
    #     check_threshold_error() pins as 'a spread of $5.9$' (since 2026-11-16 the certificate's crossing, $5.8$)
    got['crossing spread'] = '%.4f -> "within the observed spread"' % sp
    want.append(('sm_frontier.tex', 'That the two track each other to within the observed spread across five decades of error'))
    # (4) Sec. S7 opening: the fraction exponent moves by the factor of Table tab:fracgrid (eta_sweep.json)
    sys.path.insert(0, SCRIPTS)
    try:
        import fracgrid_numbers as rf
        fac = rf.numbers()['factor']
    finally:
        sys.path.pop(0)
    got['fraction-exponent factor'] = '%.4f' % fac
    want.append(('sec_7.tex', 'its exponent moving by a factor $%.1f$ over a decade of the resolution' % fac))
    # (5) gap-scaling caption: the two spin channels agree to within max |1 - ratio| (L >= 6), rounded UP;
    #     2026-11-14 (closing reading A2): the ratio that sets it is printed with four decimals, the others with
    #     three, so that the printed list shows the bound (1 - 0.917 read as 8.3 % against 'within 8.4 %')
    gs = J('gap_scaling.json')
    sizes = gs['provenance']['parameters']['sizes']
    rat = [r for L, r in zip(sizes, gs['spin_channel_control']['ratio_hubbard_over_heisenberg']) if L >= 6]
    gp = up(100 * max(abs(1 - r) for r in rat), 1)
    worst = max(rat, key=lambda r: abs(1 - r))
    lst = ', '.join(('%.4f' if r == worst else '%.3f') % r for r in rat)
    got['spin channels'] = '%.4f%% -> within %.1f%%; ratios %s' % (100 * max(abs(1 - r) for r in rat), gp, lst)
    want.append(('figs/fig_gapscaling_caption.tex', 'For $L\\ge6$ the two channels agree to within $%.1f\\%%$' % gp))
    want.append(('figs/fig_gapscaling_caption.tex', '(ratios $%s$)' % lst))
    if DIRECTED_GAP != {'%.1f' % gp, '%.4f' % worst}:
        want.append(('figs/fig_gapscaling_caption.tex', 'DIRECTED_GAP is not {%r, %r}' % ('%.1f' % gp, '%.4f' % worst)))
    # (6) Lambda along the deposited depth series at L=12: rises by up to X % between consecutive depths
    #     below depth 100 (pt_L12_*_PUB.json, Lambda_K), X rounded UP
    rise = []
    for fn in glob.glob(os.path.join(DATA, 'c3_frontier', 'published_fraction', 'pt_L12_*_PUB.json')):
        x = json.load(io.open(fn, encoding='utf-8'))
        x = x[0] if isinstance(x, list) else x
        for eta in set(r['eta'] for r in x['rows']):
            s = sorted((r for r in x['rows'] if r['eta'] == eta), key=lambda r: r['n'])
            rise += [b['Lambda_K'] / a['Lambda_K'] - 1 for a, b in zip(s, s[1:]) if b['n'] <= 100]
    rp = up(100 * max(rise), 1)
    got['depth rise L=12'] = '%.4f%% -> up to %.1f%%' % (100 * max(rise), rp)
    want.append(('sm_frontier.tex', 'rises by up to $%.1f\\%%$ from one depth to the next below depth $100$ at $L=12$' % rp))
    # (7) Sec. S8: this work's eta L against the dynamical-DMRG constants it quotes (printed numbers)
    ours = [12 * e for e in (0.10, 0.15, 0.18)]
    consts = (5.0, 12.0, 6.4, 12.8)
    words = {1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five', 6: 'six', 7: 'seven', 8: 'eight', 9: 'nine',
             10: 'ten', 11: 'eleven', 12: 'twelve'}
    lo, hi = math.floor(min(consts) / max(ours)), math.ceil(max(consts) / min(ours))
    got['eta L factor'] = '%.2f..%.2f -> %s to %s' % (min(consts) / max(ours), max(consts) / min(ours),
                                                      words[lo], words[hi])
    want += [('sm_physics.tex', '$\\eta L=1.2$--$2.2\\,t$, a factor of %s to %s below every constant' % (words[lo], words[hi])),
             ('sm_physics.tex', 'use $\\eta=5t/L$ in one place'), ('sm_physics.tex', 'declare $\\eta L=12\\,t$'),
             ('sm_physics.tex', 'with $c=6.4$ and $12.8\\,t$')]
    # (8) 2026-11-14 (closing reading A3): two ceilings the 7b audit saw and left -- the ladder/chain support ratio
    #     (resource_master.json, max 2.2374 -> "up to 2.24x") and the CASCI residual (n19_suite.json, max
    #     1.137e-13 -> "converged to within 1.2e-13 Ha"); both rounded UP at the printed digits
    rm = J('resource_master.json')
    lad = [x for x in rm['points'] if x['family'] == 'ladder']
    ch2 = [x for x in rm['points'] if x['family'] == 'chain2']
    lr = up(max(lad[k]['S'] / ch2[k]['S'] for k in range(3)), 2)
    n19 = J('n19_suite.json')
    res = 1e13 * max(max(abs(m_.get('dFCI_eq') or 0), abs(m_.get('dFCI_diss') or 0)) for m_ in n19['mols'])
    # 2026-11-14 (closing round 2, D3): the residual is printed one way in all four places, at three figures
    #     (1.137 -> 1.14), where the ceiling ("converged to within", two captions) and the nearest value (the
    #     appendix, twice) agree; if they ever differ, the four printings must be re-judged
    rs = up(res, 2)
    if rs != round(res, 2):
        want.append(('app_carried_repro.tex', 'the CASCI residual rounds differently up and to nearest: re-judge'))
    got['ladder/chain, CASCI residual'] = 'up to %.2f; within %.2fe-13 Ha (nearest %.2f)' % (lr, rs, round(res, 2))
    want += [('sm_nogo.tex', "support exceeds the chain's by up to $%.2f\\times$" % lr),
             ('sm_nogo.tex', 'differing by up to $%.2f\\times$' % lr),
             ('table_molecules_v3.tex', 'converged to within $%.2f\\times10^{-13}$~Ha' % rs),
             ('carried/fig_hero.tex', 'converged to within $%.2f\\times10^{-13}$~Ha' % rs),
             ('app_carried_repro.tex', 'the largest residual across all $38$ points is $%.2f\\times10^{-13}$~Ha'
              % round(res, 2)),
             ('app_carried_repro.tex', 'the residual the appendix records is $%.2f\\times10^{-13}$~Ha' % round(res, 2))]
    bad = [(f, w_) for f, w_ in want if w_ not in _printed(f)]
    # the nearest-rounded literals these replace must not come back
    gone = [('sec_5_body.tex', 'at most $7\\%$'), ('sm_frontier.tex', 'at most $7\\%$'),
            ('sm_frontier.tex', 'within a factor $5.9$'),
            ('sec_7.tex', 'by a factor $2.2$ over a decade'), ('figs/fig_gapscaling_caption.tex', 'within $8.3\\%$'),
            ('sm_frontier.tex', 'up to $2.2\\%$'), ('sm_physics.tex', 'two to ten below'),
            ('sm_physics.tex', 'at least $1.08$'), ('sec_5_body.tex', '0.717'), ('sm_frontier.tex', '0.717'),
            ('figs/fig5_caption.tex', 'at most $8.1\\%$'), ('figs/captions.tex', 'within $0.9\\%$'),
            ('sec_7.tex', '[-0.02,0.19]'), ('figs/caption_fig12.tex', '[-0.02,0.19]'),
            ('sec_7.tex', '0.172]'), ('figs/caption_fig12.tex', '0.172]'),
            ('sec_3_body.tex', 'by at most a factor $1.7$'), ('sm_frontier.tex', 'by at most a factor $1.7$'),
            ('sec_2.tex', 'up to $2.60\\times10^{6}$'),
            # 2026-11-14 (closing reading A1-A3): the L=14 floor wording, the second printing of the crossing spread,
            # the integer ceiling of the threshold move, the 3-decimal worst ratio, 'up to 2.61e6', '~2.2x', 1.1e-13
            ('sec_7.tex', 'floor of $'), ('sec_7.tex', 'floor, not measured'), ('sec_7.tex', '\\gtrsim1.2'),
            ('sec_7.tex', 'which is a lower bound and not an estimate'), ('sm_frontier.tex', 'within a factor $6.0$'),
            ('sec_5_body.tex', 'at most $8\\%$'), ('sm_frontier.tex', 'at most $8\\%$'),
            ('figs/fig_gapscaling_caption.tex', '0.917,'), ('sec_2.tex', 'up to $2.61\\times10^{6}$'),
            ('sm_nogo.tex', '{\\sim}2.2\\times'), ('table_molecules_v3.tex', 'converged to $1.1\\times10^{-13}$'),
            ('carried/fig_hero.tex', 'converged to $1.1\\times10^{-13}$'),
            # 2026-11-14 (closing round 2, D3): the second printing of the residual must not come back
            ('table_molecules_v3.tex', 'within $1.2\\times10^{-13}$'), ('carried/fig_hero.tex', 'within $1.2\\times10^{-13}$'),
            ('app_carried_repro.tex', '$1.1\\times10^{-13}$'),
            # 2026-11-14 (closing round 3, D2): the fit's L=14 value must not come back into the comparison of
            # Sec. 7.3 (the lower end of its interval lies below the budget drawn), nor without its interval
            ('sec_7.tex', 'the log-linear fit gives'),
            ('sec_7.tex', 'the log-linear fit itself gives ${\\approx}2.2\\times10^{8}$).'),
            # 2026-11-16 (second-order branch): the first-order threshold move and floor must not come back
            ('sec_5_body.tex', 'at most $7.2\\%$'), ('sm_frontier.tex', 'at most $7.2\\%$'),
            ('sec_5_body.tex', '0.9280'), ('sm_frontier.tex', '0.9280'), ('sm_frontier.tex', '<0.928$')]
    bad += [(f, 'still prints ' + w_) for f, w_ in gone if w_ in _printed(f)]
    for f, w_ in bad:
        print('  FAIL  %s does not print %r' % (f, w_) if not w_.startswith('still') else '  FAIL  %s %s' % (f, w_))
    print('  %s  directed rounding of bounds: %s' % ('ok  ' if not bad else 'FAIL',
                                                     '; '.join('%s %s' % kv for kv in got.items())))
    _count(len(want) + len(gone))
    return not bad


def check_restatements():
    bad = []
    for tex, lit, anchor, phrase in RESTATED:
        if phrase not in _printed(anchor):
            bad.append('%s no longer prints %r' % (anchor, phrase))
        if not re.search(r'(?<![\d.])' + re.escape('$' + lit + '$'), _printed(tex)):
            bad.append('%s does not print $%s$' % (tex, lit))
    for b in bad:
        print('  FAIL  %s' % b)
    print('  %s  restatements: %d numbers equal to the verify.py-checked figure they restate'
          % ('ok  ' if not bad else 'FAIL', len(RESTATED)))
    _count(2 * len(RESTATED))
    return not bad


# =====================================================================================================================
# 2026-11-16: checks of the NEW statements of the second-order branch (rule D1-D23,
# the number sheet).  Each has a data half (_expect_*: stored primitives and re-scored outputs only, run by
# the unit checks of the package's guards without the manuscript) and a text half (_verdict).
# =====================================================================================================================
def _expect_tau():
    """Data half of check_second_order_tau() (see there)."""
    import io
    import json
    import math
    exp = dict(data_bad=[], want=[], gone=[], regex=[], near=[], nsp_want=[], nsp_gone=[], nsp_any=[])
    n, prev = 200000, 0.0
    for i in range(1, n + 1):
        t = _tau2(i / float(n))
        if not t > prev:
            exp['data_bad'].append('tau(w) = sqrt(1+w)-sqrt(1-w) is not positive and increasing at w = %.6f' % (i / float(n)))
            break
        prev = t
    if abs(_tau2(1.0) - math.sqrt(2.0)) > 1e-15:
        exp['data_bad'].append('tau(1) is not sqrt 2')
    g, rows = _grid_rows()
    mism = [r['idx'] for r in rows if (r['B2'] < r['T']) != (r['Lam_over_eta'] < _tau2(r['w_S']))]
    if mism:
        exp['data_bad'].append('B2 < 1+w and x < tau(w) disagree on grid rows %r' % mism[:5])
    pub = json.load(io.open(os.path.join(DATA, 'c3_frontier', 'published_fraction', 'PUB_rows.json'), encoding='utf-8'))
    w12 = [r['w_S'] for r in pub if r['L'] == 12][0]
    ref = _so('supplement.json')['tau']
    if abs(ref['tau_new_L12pub'] - _tau2(w12)) > 1e-12 or abs(ref['tau_new_0995'] - _tau2(0.995)) > 1e-12:
        exp['data_bad'].append('tau values disagree with second_order/supplement.json')
    t12 = _rnd(_tau2(w12), 4, 'near', 'dp')
    for rel in ('sec_5_body.tex', 'sm_frontier.tex'):
        exp['want'].append((rel, '=%s' % t12))
        exp['gone'].append((rel, '=0.9625'))
    exp['nsp_want'] += [('sec_5_body.tex', '\\tau(w_S):=\\sqrt{1+w_S}-\\sqrt{1-w_S}'),
                        ('sec_5_body.tex', '\\tau(1)=\\sqrt'), ('sm_frontier.tex', '\\tau(1)=\\sqrt')]
    exp['nsp_gone'] += [('sec_5_body.tex', '\\frac{1+w_S}{1+\\sqrt{w_S}}-\\sqrt{1-w_S}'),
                        ('sec_5_body.tex', '\\tau(1)=1,'), ('sm_frontier.tex', '\\tau(1)=1$')]
    exp['gone'] += [('sec_5_body.tex', 'drops below unity'), ('sm_frontier.tex', 'for $w_S\\ge0.067$'),
                    ('sm_frontier.tex', 'positive for $w_S>0.296$')]
    exp['nsp_any'].append(('sec_6_body.tex', ['W/(2\\sqrt2\\eta)', 'W/(2\\sqrt{2}\\eta)'],
                           'the threshold of Prop. 2 for real problems, W/(2 sqrt2 eta) (P1)'))
    exp['label'] = ('threshold of the certificate tau(w) = sqrt(1+w)-sqrt(1-w): positive and increasing on (0,1], '
                    'tau(1) = sqrt 2, equivalent to B2 < 1+w on the 204 grid rows; tau(w_S) = %s at the L=12 published '
                    'subspace' % t12)
    return exp


def check_second_order_tau():
    """2026-11-16 (second-order branch; rule D1, D6, D12, D22; the number sheet rows V3, V4, V8, F-4, F-43, F-44,
    P1): the threshold of the certificate for a real probe, B2 < 1+w exactly when Lambda^/eta < tau(w) :=
    sqrt(1+w) - sqrt(1-w).  Data: tau is positive and increasing on (0, 1] (200000 points), tau(1) = sqrt 2, and the
    equivalence holds on every grid row of C3_grid.json (both sides recomputed).  Text: Eq. (tau) of Sec. V defines the
    new tau with tau(1) = sqrt2 and the first-order definition and tau(1)=1 are gone ('drops below unity' too, V4); the
    L=12 published subspace asks for tau(w_S)=1.3766 (stored w_S, four decimals, nearest) in Secs. V and S4, where
    0.9625 is gone; Sec. S4 no longer says tau is decreasing for w_S >= 0.067 or positive for w_S > 0.296 (F-43); Prop. 2
    states W/(2 sqrt2 eta) for real problems (P1; W/2eta may stay for the first-order branch).  The floor 1.341 and the
    5.2 % move are checked by check_directed_rounding()."""
    exp = _expect_tau()
    return _verdict(exp, exp['label'])


def _expect_table_S1():
    """Data half of check_table_S1() (see there)."""
    exp = dict(data_bad=[], want=[], gone=[], regex=[], near=[])
    g, rows = _grid_rows()
    ref = {(t['L'], t['eta']): t for t in _so('supplement.json')['table_S1']}

    def mark(a, b):
        c = (a['conv'], b['conv'])
        return 'n.c.' if 'NOT-CONVERGED' in c else ('' if c == ('CONVERGED', 'CONVERGED') else '*')

    def cell(a, b):
        mk = mark(a, b)
        if mk == 'n.c.':
            return '\\emph{n.c.}', 'n.c.'
        return ('$[%.2f,%.2f]%s$' % (a['FR'], b['FR'], '^{\\ast}' if mk else ''),
                '[%.2f,%.2f]%s' % (a['FR'], b['FR'], mk))
    printed = []
    for eta in (0.18, 0.15):
        for L in (6, 8, 10, 12):
            cl = sorted((r for r in rows if r['L'] == L and abs(r['eta'] - eta) < 1e-12), key=lambda r: r['FR'])
            fs, lo, hi = _fsum_frontier([r['B2'] for r in cl], [r['FR'] for r in cl], [r['T'] for r in cl])
            cols = [cell([r for r in cl if r['FR'] == lo][0], [r for r in cl if r['FR'] == hi][0])]
            for lev in (0.5, 0.2):
                i = [k for k, r in enumerate(cl) if r['c2'] < lev]
                if not i:
                    cols.append(('$>%.2f$' % cl[-1]['FR'], '>%.2f' % cl[-1]['FR']))
                elif i[0] == 0:
                    cols.append(('$\\le%.2f$' % cl[0]['FR'], '<=%.2f' % cl[0]['FR']))
                else:
                    cols.append(cell(cl[i[0] - 1], cl[i[0]]))
            t = ref.get((L, eta))
            if t is None or [c[1] for c in cols] != [t['FRstar_new'], t['cert_lt_050_new'], t['cert_lt_020_new']]:
                exp['data_bad'].append('Table S1 L=%d eta=%.2f: %r disagrees with second_order/supplement.json' % (L, eta, cols))
            row = ' & '.join(c[0] for c in cols)
            exp['want'].append(('sm_frontier.tex', row))
            printed.append(row)
    exp['printed_rows'] = printed
    exp['label'] = 'Table S1, columns 5-7 (FR* and the certificate below 0.50 and 0.20), 8 rows from the stored grid'
    return exp


def check_table_S1():
    """2026-11-16 (second-order branch; rule D6, D19; the number sheet rows F-13a-d, F-14a-d, F-20; check
    discrepancy 4): Table S1 (tab:thresholds), columns 5-7, re-derived from the stored grid rows with the certificate of
    D6: column 5 the bracket of FR* (fsum.py convention), columns 6-7 where min{1+w, B2} first falls below 0.50 and
    0.20 in grid order; '$[lo,hi]$', '^{\\ast}' if an end is MARGINAL, '\\emph{n.c.}' if an end is NOT-CONVERGED,
    '$>FR$' if not reached on the grid.  Each of the eight rows (eta = 0.18 and 0.15 t, L = 6..12) must be printed
    verbatim and agree with second_order/supplement.json table_S1.  (The number_sources records of this table were keyed to the
    stored first-order FR_star by a coincidence of the printed digits; they are re-keyed to the printed column.)"""
    exp = _expect_table_S1()
    return _verdict(exp, exp['label'])


def _expect_akw():
    """Data half of check_akw_channels() (see there)."""
    import io
    import json
    exp = dict(data_bad=[], want=[], gone=[], regex=[], near=[])
    ca = [r for r in json.load(io.open(os.path.join(DATA, 'cert_akw.json'), encoding='utf-8'))['rows']
          if r['L'] == 8 and abs(r['FR'] - 0.85) < 1e-9]
    real = [r for r in ca if abs(r['k_over_pi']) < 1e-9 or abs(r['k_over_pi'] - 1.0) < 1e-9]
    cplx = [r for r in ca if not (abs(r['k_over_pi']) < 1e-9 or abs(r['k_over_pi'] - 1.0) < 1e-9)]
    rc, cc = [], []
    for r in real:
        w, x = r['w_S'], r['Lambda_hat_over_eta']
        if abs(_B1(w, x) - r['rel_upper']) > 1e-6 * r['rel_upper']:
            exp['data_bad'].append('cert_akw.json k/pi=%g: rel_upper is not B1(w_S, Lambda_hat_over_eta)' % r['k_over_pi'])
        c = min(1.0 + w, _B2(w, x))
        rc.append((c, c / r['relL1_true'], c < 1.0 + w))
    for r in cplx:
        cc.append((r['rel_upper'], r['slack'], bool(r['bound_nontrivial'])))
    if (len(ca), len(real), len(cplx)) != (16, 4, 12) or not all(x[2] for x in rc + cc):
        exp['data_bad'].append('cert_akw.json L=8 FR=0.85: not 4 real + 12 complex non-vacuous channels')
    nd = lambda rs: len(set((round(r['w_S'], 8), _rnd(r['Lambda_hat_over_eta'], 4)) for r in rs))
    if (nd(real), nd(cplx)) != (2, 3):
        exp['data_bad'].append('symmetry-distinct channels %d + %d, not 2 + 3 ("five of them")' % (nd(real), nd(cplx)))
    d18 = _so('v2.json')['D18']
    for got, want in ((min(x[0] for x in rc), d18['real']['cert'][0]), (max(x[0] for x in rc), d18['real']['cert'][1]),
                      (min(x[1] for x in rc), d18['real']['cert_over_true'][0]),
                      (max(x[1] for x in rc), d18['real']['cert_over_true'][1]),
                      (min(x[0] for x in cc), d18['complex']['first_order'][0]),
                      (max(x[0] for x in cc), d18['complex']['first_order'][1]),
                      (min(x[1] for x in cc), d18['complex']['first_order_over_true'][0]),
                      (max(x[1] for x in cc), d18['complex']['first_order_over_true'][1])):
        if abs(got - want) > 1e-9 * want:
            exp['data_bad'].append('cert_akw.json ranges disagree with second_order/v2.json D18 (%.12g vs %.12g)' % (got, want))
    pr = lambda v, f: '$%s$--$%s$' % (f(min(v)), f(max(v)))
    R = pr([x[0] for x in rc], lambda y: _rnd(y, 2, 'near', 'dp'))
    RS = pr([x[1] for x in rc], lambda y: _rnd(y, 2))
    C = pr([x[0] for x in cc], lambda y: _rnd(y, 2, 'near', 'dp'))
    CS = pr([x[1] for x in cc], lambda y: _rnd(y, 2))
    for rel in ('sec_1.tex', 'sec_5_body.tex', 'sec_10.tex', 'sm_frontier.tex'):
        exp['want'] += [(rel, R), (rel, RS), (rel, C), (rel, CS)]
        exp['near'].append((rel, C, ['first-order'], 300, 300,
                            'the complex channels carry the first-order branch only and must say so (D8, D18)'))
        exp['gone'] += [(rel, '$0.77$--$1.75$'), (rel, '$350$--$900$')]
    exp['label'] = ('momentum-resolved certificate at L=8, FR=0.85: real channels k=0, pi %s (%s x the true error), '
                    'complex channels first-order %s (%s x), 16 of 16 non-vacuous, 2 + 3 symmetry-distinct'
                    % (R, RS, C, CS))
    return exp


def check_akw_channels():
    """2026-11-16 (second-order branch; rule D8, D18; the number sheet rows I6a-c, V20a-c, C4a-b, F-10a-b; the
    replacement of verify.py 9.9c 'relative bound, low/high' and 'slack, low/high (2 s.f.)' in sec_1, sec_5_body,
    sec_10 and sm_frontier, registered at the top of this file): the momentum-resolved certificate of Fig. 3, L=8,
    FR=0.85, cert_akw.json, sixteen channels.  The real channels k=0, pi carry the certificate of D6, min{1+w, B2}
    from the stored w_S and Lambda_hat_over_eta (their stored first-order rel_upper is reproduced first); the complex
    channels carry the first-order branch only (stored rel_upper and slack; Lambda~ was not stored, D8).  All sixteen are
    non-vacuous, two + three symmetry-distinct.  Printed in Secs. I, V, X and S4: the real range at two decimals and
    its ratio to the true error at two figures, the complex first-order range and ratio the same way, the latter
    labelled first-order within 300 characters; the all-channel first-order 0.77--1.75 and 350--900 must be gone."""
    exp = _expect_akw()
    return _verdict(exp, exp['label'])


def _pool():
    """The pooled sweep of Fig. S10 from its five deposited sources (src/recovered/build_n3.py), both leakage branches
    in absolute units: thm3_results.json (norm_phi2, wS, Lam, eta; RHS_new is the first-order branch, reproduced),
    thm3_L8b.json (LB = ||phi||^2 (1-w), Lam, Lam_over_eta = Lam/eta), cert_{stress,akw,teqsci}.json (norm_phi,
    norm_phiS, norm_phiQ, Lambda, eta; upper_unclipped is the first-order branch, reproduced).  A cert_akw.json row with
    k/pi other than 0 and 1 is a complex probe: its certificate is the first-order branch only (D8); every other row is
    real and its certificate is min{trivial, B2} (D6)."""
    import io
    import json
    import math
    J = lambda *p: json.load(io.open(os.path.join(DATA, *p), encoding='utf-8'))
    P = []
    for i, r in enumerate(J('thm1iii_violation', 'thm3_results.json')):
        nphi, w = math.sqrt(r['norm_phi2']), r['wS']
        nQ, nS = nphi * math.sqrt(max(1.0 - w, 0.0)), nphi * math.sqrt(w)
        if abs((nphi + nS) * (nQ + r['Lam'] / r['eta']) - r['RHS_new']) > 1e-8 * r['RHS_new']:
            raise AssertionError('thm3_results.json [%d]: RHS_new is not the first-order branch' % i)
        P.append(dict(src='thm1iii_violation/thm3_results.json [%d]' % i, err=r['LHS'], rel=r['relL1'],
                      triv=r['RHS_triv'], b1=r['RHS_new'], b2=(nQ + r['Lam'] / r['eta']) ** 2, frac=r['frac'],
                      eta=r['eta'], real=True, dep=False))
    for i, r in enumerate(J('thm1iii_violation', 'thm3_L8b.json')):
        w = r['wS']
        n2 = r['LB'] / (1.0 - w)
        nQ = math.sqrt(n2 * (1.0 - w))
        P.append(dict(src='thm1iii_violation/thm3_L8b.json [%d]' % i, err=r['LHS'], rel=r['relL1'], triv=n2 * (1.0 + w),
                      b1=r['RHS_new'], b2=(nQ + r['Lam_over_eta']) ** 2, frac=r['FR'], eta=r['Lam'] / r['Lam_over_eta'],
                      real=True, dep=False))
    for fn in ('cert_stress.json', 'cert_akw.json', 'cert_teqsci.json'):
        for i, r in enumerate(J(fn)['rows']):
            nQ = r['norm_phiQ']
            if abs((r['norm_phi'] + r['norm_phiS']) * (nQ + r['Lambda'] / r['eta']) - r['upper_unclipped']) > \
                    1e-8 * r['upper_unclipped']:
                raise AssertionError('%s /rows[%d]: upper_unclipped is not the first-order branch' % (fn, i))
            cplx = fn == 'cert_akw.json' and not (abs(r['k_over_pi']) < 1e-9 or abs(r['k_over_pi'] - 1.0) < 1e-9)
            P.append(dict(src='%s /rows[%d]' % (fn, i), err=r['L1_true'], rel=r['relL1_true'], triv=r['upper_trivial'],
                          b1=r['upper_unclipped'], b2=(nQ + r['Lambda'] / r['eta']) ** 2, frac=r['frac'],
                          eta=r['eta'], real=not cplx, dep=True))
    for p in P:
        p['cert'] = min(p['triv'], p['b2'] if p['real'] else p['b1'])
    return P


def _expect_pooled():
    """Data half of check_pooled_split() (see there)."""
    exp = dict(data_bad=[], want=[], gone=[], regex=[], near=[])
    P = _pool()
    live = [p for p in P if p['rel'] > 1e-10]
    real = [p for p in live if p['real']]
    cplx = [p for p in live if not p['real']]
    for p in live:
        p['ratio'] = p['err'] / p['cert']
    n_old = sum(p['b1'] < p['triv'] for p in live)
    n_r = sum(p['b2'] < p['triv'] for p in real)
    n_c = sum(p['b1'] < p['triv'] for p in cplx)
    viol = sum(p['ratio'] > 1.0 + 1e-9 for p in live)
    rmin, rmax = min(p['ratio'] for p in live), max(p['ratio'] for p in live)
    sm = _so('pooled.json')['summary']
    if (len(P), len(live), len(real), len(cplx), n_old, n_r, n_c, viol) != (
            sm['n_pooled'], sm['n_live'], sm['n_real_live'], sm['n_complex_live'], sm['leak_wins_old'],
            sm['real_rows_B2_wins'], sm['complex_rows_B1_wins'], sm['violations_new']):
        exp['data_bad'].append('pool counts disagree with second_order/pooled.json summary')
    if abs(rmin - sm['ratio_new'][0]) > 1e-9 * rmin or abs(rmax - sm['ratio_new'][1]) > 1e-12:
        exp['data_bad'].append('pool ratio range disagrees with second_order/pooled.json summary')
    if viol:
        exp['data_bad'].append('THE NEW CERTIFICATE IS VIOLATED on %d pooled rows' % viol)
    tot = n_r + n_c
    places = ('sec_3_app.tex', 'sec_5_body.tex', 'sm_frontier.tex', 'figs/caption_thm1iii_violation.tex')
    for rel in places:
        exp['near'].append((rel, '$%d$' % tot, ['re:(?<![\\d.])%d(?![\\d])' % n_r, 're:(?<![\\d.])%d(?![\\d])' % n_c],
                            250, 250, 'the pooled count %d must carry its split %d real + %d complex (D18)'
                            % (tot, n_r, n_c)))
        exp['gone'].append((rel, 'on only $%d$' % n_old))
    rng = '%s to $%.10f$' % (_tsci(rmin, 3), rmax)
    exp['want'] += [('figs/caption_thm1iii_violation.tex', rng),
                    ('figs/caption_thm1iii_violation.tex', 'violated on none of the $%d$' % len(live))]
    exp['gone'].append(('figs/caption_thm1iii_violation.tex', '$1.18\\times10^{-5}$'))
    exp['dat'] = {'n3_thmB_frac.dat': sorted(('%.6g' % p['frac'], '%.6g' % p['ratio']) for p in live),
                  'n3_thmB_eta.dat': sorted(('%.6g' % p['eta'], '%.6g' % p['ratio']) for p in live)}
    exp['label'] = ('pooled sweep of Fig. S10 with the new certificate: %d live rows (%d real, %d complex), %d violations, '
                    'beats the trivial bound on %d = %d + %d (first-order branch alone: %d), ratios %s'
                    % (len(live), len(real), len(cplx), viol, tot, n_r, n_c, n_old, rng))
    return exp


def check_pooled_split():
    """2026-11-16 (second-order branch; rule D6, D8, D18; the number sheet rows X1, X3, V10, F-9, S10b, S10d, S10e,
    GD8; the replacement of verify.py 9.6c 'rows on which the leakage branch is the smaller' (Fig. S10 caption and
    Sec. S9), 'Theorem 1 ratio range' and the two n3_thmB_*.dat tables, registered at the top of this file, and of the
    two RESTATED entries of the count 61): the pooled sweep of Fig. S10 (_pool(): 621 subspaces, 606 above round-off)
    with the new certificate -- real rows min{T, B2}, complex momentum rows min{T, B1}.  It is violated on none; it
    beats the trivial bound on 85 = 57 of 526 real rows + 28 of 80 complex rows (61 with the first-order branch alone);
    its error-to-bound ratios span 3.55e-4 to 0.9999999997 (agreeing with second_order/pooled.json).  Printed: '$85$' with its
    split (57 and 28 within 250 characters) in Sec. S9, Sec. V, Sec. S4 and the Fig. S10 caption, where 'on only $61$'
    is gone; the ratio range in the caption; the plotted n3_thmB_frac.dat / n3_thmB_eta.dat hold one point per live
    row, (fraction or eta, error / certificate) at six figures."""
    import io
    exp = _expect_pooled()
    for name, want in exp['dat'].items():
        got = []
        for ln in io.open(os.path.join(JCP_PAPER, 'figs', name), encoding='utf-8').read().splitlines()[1:]:
            if ln.strip():
                a, b = ln.split()
                got.append((a, b))
        if sorted(got) != want:
            exp['data_bad'].append('figs/%s does not hold the (x, error / new certificate) of the %d live rows'
                                   % (name, len(want)))
    ok = _verdict(exp, exp['label'])
    _count(LAST_N[0] + sum(len(v) for v in exp['dat'].values()))
    return ok


def _expect_fitted():
    """Data half of check_fitted_rule() (see there)."""
    exp = dict(data_bad=[], want=[], gone=[], regex=[], near=[])
    n2, n1 = _so('null_ws7_second_order.json'), _so('null_ws7_first_order_copy.json')
    cells, cbad = _second_order_cells()
    exp['data_bad'] += cbad
    for c in cells:
        k = 'L%d_eta%g' % (c['L'], c['eta'])
        if abs(n2['summary']['FRstar'][k] - c['FRstar']) > 1e-12:
            exp['data_bad'].append('%s: the re-run fitted-rule comparison did not use the certificate\'s crossing' % k)
        if abs(n1['summary']['FRstar'][k] - c['FRstar_old']) > 1e-12:
            exp['data_bad'].append('%s: the first-order control did not use the stored FR_star' % k)

    def counts(n):
        s = n['splits']
        return (len(s), sum(x['nullB_test_spread'] < x['cert_test_spread'] for x in s),
                sum(x['nullA_test_spread'] > x['cert_test_spread'] for x in s),
                [v for x in s for v in x['fitB']['firing_level_at_eta'].values()])
    N, tight, loose, lv = counts(n2)
    N1, tight1, loose1, lv1 = counts(n1)
    sm = n2['summary']
    if (N, tight, loose) != (sm['n_splits'], sm['two_constant_rule_tighter_than_cert'],
                             sm['eta_blind_rule_looser_than_cert']) or \
            [min(lv), max(lv)] != sm['two_constant_firing_level_range']:
        exp['data_bad'].append('second_order/null_ws7_second_order.json: its summary is not the count of its 14 split records')
    if (N1, tight1, loose1) != (14, 11, 13):
        exp['data_bad'].append('the first-order control no longer reproduces the deposit\'s 11 and 13 of 14')
    fl = (5e-3 / max(lv), 5e-3 / min(lv))
    if any(abs(a - b) > 1e-9 * b for a, b in zip(fl, sm['floor_5e3_over_two_constant_levels'])):
        exp['data_bad'].append('floor 5e-3 over the firing level disagrees with second_order/null_ws7_second_order.json')
    W = {10: 'ten', 11: 'eleven', 12: 'twelve', 13: 'thirteen', 14: 'fourteen'}
    if N != 14 or tight not in W or loose not in W:
        exp['data_bad'].append('fitted-rule counts outside the printed words: %d, %d of %d' % (tight, loose, N))
        return exp
    for rel in ('sec_1.tex', 'sec_5_body.tex', 'sec_10.tex', 'sm_frontier.tex'):
        exp['regex'].append((rel, r'more tightly(?: than the certificate)? on (\w+) of (?:the )?fourteen', W[tight],
                             'the count of splits on which the two-constant rule localizes the crossing more tightly'))
        exp['gone'].append((rel, '%s of' % W[tight1]))
    for rel in ('sec_5_body.tex', 'sm_frontier.tex'):
        exp['regex'].append((rel, r'looser of the two on (\w+) of the fourteen', W[loose],
                             'the count of splits on which the eta-blind rule is the looser'))
        exp['want'] += [(rel, _tsci(min(lv), 2)), (rel, _tsci(max(lv), 2))]
        exp['gone'] += [(rel, '$10^{-5}$--$10^{-4}$'), (rel, 'one to two orders of magnitude')]
    exp['want'] += [('sm_frontier.tex', '$%s$' % _rnd(fl[0], 2)), ('sm_frontier.tex', _tsci(fl[1], 2))]
    exp['label'] = ('fitted rule on the certificate\'s crossings: tighter on %s of %d splits, eta-blind looser on %s '
                    '(first-order control %d, %d); fires at %s to %s, the floor 5e-3 a factor %s to %s above it'
                    % (W[tight], N, W[loose], tight1, loose1, _tsci(min(lv), 2), _tsci(max(lv), 2), _rnd(fl[0], 2),
                       _tsci(fl[1], 2)))
    return exp


def check_fitted_rule():
    """2026-11-16 (second-order branch; rule D20; the number sheet rows I5, V30, V31, V32, C3, F-49, F-50, F-51;
    the check of the number sheet discrepancy 6, whose two-figure firing level 2.4e-5 to 1.0e-3 is pinned here, since with
    the one-figure 2e-5 the printed floor ratio would be 2.5e2, not 2.1e2): the fitted-rule comparison of Secs. V D and
    S4, re-run on the certificate's crossings by the copy null_ws7_second_order.py of the deposit's
    src/frontier/null_ws7.py (output null_ws7_second_order.json, deposit file SECOND_ORDER/null_ws7_second_order.json);
    its control, the same copy on the stored first-order FR_star (null_ws7_first_order_copy.json), reproduces the
    deposit's own 11 and 13 of 14.  Checked: the crossings it used are the certificate's (_second_order_cells()) and the
    control's the stored ones; its summary is the count of its fourteen split records; the text prints the
    two-constant rule tighter on thirteen of fourteen splits (Secs. I, V, X, S4), the eta-blind rule looser on twelve
    (V, S4), the firing level 2.4e-5 to 1.0e-3 (V, S4) and the floor 5e-3 as a factor 5.0 to 2.1e2 above it (S4); the
    first-order 'eleven of', the eta-blind 'thirteen', '10^{-5}--10^{-4}' and 'one to two orders of magnitude' are gone."""
    exp = _expect_fitted()
    return _verdict(exp, exp.get('label', 'fitted rule'))


def _expect_never_below():
    """Data half of check_never_below() (see there)."""
    import glob
    import io
    import json
    exp = dict(data_bad=[], want=[], gone=[], regex=[], near=[])
    J = lambda *p: json.load(io.open(os.path.join(DATA, *p), encoding='utf-8'))
    cen = []

    def add(src, rel, err, b2, b1, real=True, full=False):
        cls = ('full' if full else 'noerr' if rel is None else 'roundoff' if rel <= 1e-10 else
               'complex' if not real else 'tested')
        cen.append(dict(src=src, rel=rel, err=err, b2=b2, b1=b1, cls=cls))
    g, rows = _grid_rows()
    for r in rows:
        add('C3_grid.json', r['relL1'], r['relL1'], r['B2'], r['B1'])
    for f in sorted(glob.glob(os.path.join(DATA, 'c3_frontier', 'published_fraction', 'pt_L*_PUB.json'))):
        d = json.load(io.open(f, encoding='utf-8'))
        d = d[0] if isinstance(d, list) else d
        for q in d['rows']:
            if 'relL1' in q:
                add(os.path.basename(f), q['relL1'], q['relL1'], _B2(d['w_S'], q['LamHat_over_eta']),
                    _B1(d['w_S'], q['LamHat_over_eta']))
    l14 = J('c3_frontier', 'published_fraction', 'L14_FR008.json')[0]
    for q in l14['rows']:
        if 'relL1' in q:
            add('L14_FR008.json', q['relL1'], q['relL1'], _B2(l14['w_S'], q['LamHat_over_eta']),
                _B1(l14['w_S'], q['LamHat_over_eta']))
    zf = io.open(os.path.join(DATA, 'c3_frontier', 'adversarial', 'z_fine_L8.txt'), encoding='utf-8').read()
    for m in re.finditer(r"\|S\|=(\d+) FR=[\d.]+ eta=0\.18 \| w_S=([\d.e+-]+) relR0=([\d.e+-]+) "
                         r"relR5=([\d.e+-]+) Lamhat/eta=([\d.e+-]+)", zf):
        w, e = float(m.group(2)), float(m.group(3))
        add('z_fine_L8.txt', e, e, _B2(w, float(m.group(5))), float(m.group(4)))
    for q in J('c3_frontier', 'adversarial', 'z_fine_L6.json'):
        add('z_fine_L6.json', q['relR0'], q['relR0'], (q['phiQ_rel'] + q['Lamhat_over_eta']) ** 2, q['relR5'])
    for z in J('c3_frontier', 'adversarial', 'z_rank_L6.json'):
        add('z_rank_L6.json', z['relR0'], z['relR0'], _B2(z['w_S'], z['Lamhat_over_eta']), _B1(z['w_S'], z['Lamhat_over_eta']))
    for res in J('support_witness.json')['results']:
        w = min(res['w_S'], 1.0)
        for t in res['table']:
            add('support_witness.json', t['relL1'], t['relL1'], _B2(w, t['Lambda_hat_over_eta']),
                _B1(w, t['Lambda_hat_over_eta']))
    for p in _pool():
        add('pool', p['rel'], p['err'], p['b2'] if p['real'] else None, p['b1'], real=p['real'])
    for fn in ('cert_frontier_6-8.json', 'cert_frontier_10.json'):
        for Lk, cv in J(fn)['curves'].items():
            for q in cv['points']:
                full = bool(q.get('S_is_full_sector'))
                add(fn, q.get('relL1_true'), q.get('L1_true'), None if full else (q['norm_phiQ'] + q['Lambda'] / q['eta']) ** 2,
                    None if full else q['upper_unclipped'], full=full)
    for Lf in (4, 6, 8):
        d = J('c3_frontier', 'ladder_L%d.json' % Lf)
        for c in d['cases']:
            n2 = c.get('norm_phi2', d.get('norm_phi2'))
            for q in c['rows']:
                add('ladder_L%d.json' % Lf, q['R0'] / n2, q['R0'], (c['norm_phiQ'] + q['Lam_over_eta']) ** 2, q['R5'])
    by = lambda k: [q for q in cen if q['cls'] == k]
    tested, cplx = by('tested'), by('complex')
    viol = [q for q in tested if q['err'] > q['b2']]
    viol1 = [q for q in cplx if q['err'] > q['b1']]
    got = dict(n_rows_all=len(cen), n_tested=len(tested), violations_strict=len(viol),
               n_excluded_roundoff=len(by('roundoff')), n_excluded_full_sector=len(by('full')),
               n_no_stored_error=len(by('noerr')), complex_rows=len(cplx), complex_violations_B1=len(viol1))
    d21 = _so('v2.json')['D21']
    for k, v in got.items():
        if d21[k] != v:
            exp['data_bad'].append('D21 census %s: recomputed %d, second_order/v2.json %r' % (k, v, d21[k]))
    if viol or viol1:
        exp['data_bad'].append('THE CERTIFICATE IS BELOW THE STORED TRUE ERROR on %d real / %d complex rows'
                               % (len(viol), len(viol1)))
    mr = min(q['b2'] / q['err'] for q in tested)
    # the 378 deposited rows of Sec. S9 (cert_*.json): those that define a ratio, and the certificate on each
    dep = [p for p in _pool() if p['dep']]
    dep_def = [p for p in dep if p['rel'] > 1e-10]
    held = sum(1 for p in dep_def if p['err'] <= p['cert'])
    if (len(dep), len(dep_def), held) != (378, d21['deposited_378']['n_ratio_defined'],
                                          d21['deposited_378']['each_leakage_branch_holds_on_ratio_defined']):
        exp['data_bad'].append('the 378 deposited rows: %d define a ratio, the certificate holds on %d' % (len(dep_def), held))
    # the check of the number sheet section 5: the grid-point files (pt_L*_FR*_nl*.json, not the published ones) hold more
    # rows with a stored true error; the second-order branch must hold there too (counted apart from the census)
    extra, ev, emin = 0, 0, None
    for f in sorted(glob.glob(os.path.join(DATA, 'c3_frontier', 'pt_L*_FR*_nl*.json'))):
        if f.endswith('_PUB.json'):
            continue
        d = json.load(io.open(f, encoding='utf-8'))
        d = d[0] if isinstance(d, list) else d
        for q in d['rows']:
            if q.get('relL1') and q['relL1'] > 1e-10:
                b2 = _B2(d['w_S'], q['LamHat_over_eta'])
                extra += 1
                ev += q['relL1'] > b2
                emin = b2 / q['relL1'] if emin is None else min(emin, b2 / q['relL1'])
    if ev:
        exp['data_bad'].append('the second-order branch is below the stored true error on %d grid-point-file rows' % ev)
    exp['want'] += [('sec_3_app.tex', '$%d$' % len(dep_def)), ('sec_3_app.tex', 'round-off')]
    exp['gone'].append(('sec_3_app.tex', '$378$ for each side of Eq.~\\eqref{eq:twobudget}'))
    exp['info'] = dict(got, min_ratio=mr, deposited=(len(dep), len(dep_def), held), grid_point_rows=(extra, ev, emin))
    exp['label'] = ('never below the true error: %d stored rows, %d real rows above 1e-10 tested, 0 below (smallest '
                    'B2/error %.10f); %d round-off and %d full-sector rows excluded, %d complex rows first-order (0 below), '
                    '%d without an error; deposited rows %d, %d define a ratio, certificate holds on %d; grid-point files '
                    '%d more rows, %d below (smallest %.1f)'
                    % (len(cen), len(tested), mr, got['n_excluded_roundoff'], got['n_excluded_full_sector'], len(cplx),
                       got['n_no_stored_error'], len(dep), len(dep_def), held, extra, ev, emin or 0.0))
    return exp


def check_never_below():
    """2026-11-16 (second-order branch; rule D21; the number sheet header D21 and row X4; the check of the number sheet
    section 5): 'never below the true error', with its qualifier.  Every stored row with a true error (1194 rows: the
    grid, the published-subspace files, the fine scans, the oracle ranking, the support witness, the pooled sweep, the
    certificate-frontier files and the ladders) is classified as in rescore_second_order.py section 18: rows whose true
    relative L1 error is at most 1e-10 (round-off) and the two full-sector rows are excluded and counted, the complex
    momentum rows are tested with the first-order branch only (D8), rows without a stored error are counted.  The
    second-order branch must be at or above the stored true error on every other real-probe row (1059, strict) and the
    first-order one on the 80 complex rows; the census agrees with second_order/v2.json D21.  The grid-point files hold more
    rows with a true error: the branch must hold there too.  On the 378 deposited rows of Sec. S9 the certificate holds
    on the 376 that define a ratio.  Printed (sec_3_app.tex): '$376$' and 'round-off'; '$378$ for each side of
    Eq. (twobudget)' is gone (X4)."""
    exp = _expect_never_below()
    return _verdict(exp, exp['label'])


def _expect_d17():
    """Data half of check_d17_named() (see there)."""
    import io
    import json
    import math
    exp = dict(data_bad=[], want=[], gone=[], regex=[], near=[])
    J = lambda *p: json.load(io.open(os.path.join(DATA, *p), encoding='utf-8'))
    xb = lambda B, w: B / (1.0 + math.sqrt(w)) - math.sqrt(max(1.0 - w, 0.0))
    best = 0.0
    for blk in J('c3_frontier', 'ties.json'):
        for p in blk['points']:
            for k in sorted(p):
                if k.startswith('f32_eta') and isinstance(p[k], dict):
                    k64 = k.replace('f32', 'f64')
                    w32, w64, b32, b64 = p['f32_w_S'], p['f64_w_S'], p[k]['B'], p[k64]['B']
                    s32, s64 = min(b32, _B2(w32, xb(b32, w32))), min(b64, _B2(w64, xb(b64, w64)))
                    best = max(best, abs(s32 - s64) / s64)
    ref = _so('v2.json')['D17']['ties']['max_smaller_rel']
    if abs(best - ref) > 1e-9 * ref:
        exp['data_bad'].append('ties: smaller-branch change %.12g, second_order/v2.json %.12g' % (best, ref))
    tie = _tsci(best, 2, 'up')
    for rel in ('sec_5_app.tex', 'sm_leakage.tex'):
        exp['near'].append((rel, tie, ['branch'], 250, 0,
                            'the tie effect is a change of the smaller leakage branch and must name it (D17; check '
                            'discrepancy 3)'))
        exp['gone'].append((rel, '$7.8\\times10^{-4}$'))
    p12 = [r for r in J('c3_frontier', 'published_fraction', 'PUB_rows.json') if r['L'] == 12 and abs(r['eta'] - 0.18) < 1e-12][0]
    if not p12['B'] < _B2(p12['w_S'], xb(p12['B'], p12['w_S'])):
        exp['data_bad'].append('L=12 published, eta=0.18: the first-order branch is no longer the smaller')
    exp['near'].append(('sm_frontier.tex', '$%.3f$' % p12['B'], ['first-order'], 250, 100,
                        'the L=12 published-subspace value names its branch, the first-order one (D17, F-6)'))
    l14 = J('c3_frontier', 'published_fraction', 'L14_FR008.json')[0]
    deep = max(r['n'] for r in l14['rows'])
    d14 = sorted((r for r in l14['rows'] if r['n'] == deep), key=lambda r: r['eta'])
    if not all(r['bound_leak'] < _B2(l14['w_S'], r['LamHat_over_eta']) for r in d14):
        exp['data_bad'].append('L=14 published, depth %d: the first-order branch is no longer the smaller' % deep)
    exp['near'].append(('sm_frontier.tex', '$%.2f$' % d14[0]['bound_leak'], ['first-order'], 250, 150,
                        'the L=14 published-subspace values name their branch, the first-order one (D17, F-11)'))
    exp['want'].append(('sec_3_body.tex', 'smaller of the two'))
    exp['label'] = ('values quoted above the trivial bound name their branch (D17): ties %s (smaller branch), L=12 %.3f '
                    'and L=14 %.2f.. (first-order, the smaller there); Sec. III states the rule'
                    % (tie, p12['B'], d14[0]['bound_leak']))
    return exp


def check_d17_named():
    """2026-11-16 (second-order branch; rule D17, working rule 7b; the number sheet rows T0, F-6, F-11, H2, G2;
    the check of the number sheet discrepancy 3): a leakage value quoted above the trivial bound is the SMALLER of the two
    leakage branches; Sec. III says so once ('smaller of the two'), and each such value names its branch where it is
    printed.  Data: the largest relative change that ties at the cut make in the smaller branch (c3_frontier/ties.json,
    both precisions, x recovered from the stored first-order B), 1.077e-3, printed rounded UP at two figures (1.1e-3)
    in Secs. S15 and S2 with 'branch' within 250 characters before it, where 7.8e-4 (the first-order branch alone) is
    gone; the first-order branch is the smaller at the L=12 published subspace (8.544) and at the four L=14 values
    (17.21, ...), printed in Sec. S4 with 'first-order' within 250 characters before (or 100-150 after)."""
    exp = _expect_d17()
    return _verdict(exp, exp['label'])


# The checks of the JCP text, in the order the package ran them (check_dequant_moved excepted), as (label, function).
# src/verify.py section 12 runs this list; each is one check with LAST_N assertions.
CHECKS = [
    ('fig:results (Fig. 1): plotted rows, marks and caption numbers', 'check_fig_results'),
    ('the certificate reaches 0.05 at eta=0.18t: L=8, 10, 12 brackets; not reached at L=6', 'check_fraction_five_percent'),
    ('fixed-subspace looseness of the second-order branch', 'check_fixedP_looseness'),
    ('reference-drift rule of the published subspaces (Sec. S4)', 'check_drift_margin'),
    ('withdrawn unsourced numbers stay out; their sourced replacements hold', 'check_withdrawn_unsourced'),
    ('main-text numbers: Table III (smaller leakage branch), Fig. 3 caption, Sec. III D, Sec. IV, Table II',
     'check_main_numbers'),
    ('src/jcp/number_sources.json: every registered printed number re-read from its stored source',
     'check_number_sources'),
    ('gap-scaling caption numbers produced by src/recovered/make_fig_gapscaling.py', 'check_gapscaling_caption'),
    ('Table tab:fracgrid and its quoting passages, from data/eta_sweep.json', 'check_fracgrid'),
    ('threshold anatomy (first-order chain, labelled)', 'check_frontier_anatomy'),
    ("true error at the certificate's crossing (sixteen cells)", 'check_threshold_error'),
    ('restatements of a number another check anchors', 'check_restatements'),
    ('printed arithmetic of App. dequant', 'check_printed_arithmetic'),
    ('directed rounding of bounds', 'check_directed_rounding'),
    ("threshold of the certificate tau(w) = sqrt(1+w) - sqrt(1-w)", 'check_second_order_tau'),
    ('Table S1, columns 5-7, from the stored grid', 'check_table_S1'),
    ('momentum-resolved certificate of Fig. 4 (real and complex channels)', 'check_akw_channels'),
    ('pooled sweep of Fig. S10 with the new certificate (85 = 57 + 28; the n3_thmB tables)', 'check_pooled_split'),
    ('fitted-rule comparison on the certificate\'s crossings (null_ws7_second_order.py)', 'check_fitted_rule'),
    ('never below the true error (the census of every stored row)', 'check_never_below'),
    ('values quoted above the trivial bound name their branch (D17)', 'check_d17_named'),
]


def run_all(echo=True):
    """Run every check of CHECKS; returns [(label, function name, ok, assertions, printed lines)]."""
    import contextlib
    import io
    out = []
    for label, name in CHECKS:
        buf = io.StringIO()
        LAST_N[0] = 0
        try:
            with contextlib.redirect_stdout(buf):
                ok = bool(globals()[name]())
        except Exception as exc:                                            # noqa: BLE001
            ok = False
            buf.write('  FAIL  %s raised %s: %s\n' % (name, type(exc).__name__, exc))
        lines = buf.getvalue().splitlines()
        if echo:
            for ln in lines:
                print(ln)
        out.append((label, name, ok, max(LAST_N[0], 1), lines))
    return out


def main():
    res = run_all()
    bad = [r for r in res if not r[2]]
    n = sum(r[3] for r in res)
    if bad:
        print('RESULT: FAIL -- %d of %d checks of the JCP text fail: %s' % (len(bad), len(res), [r[1] for r in bad]))
        return 1
    print('RESULT: PASS -- %d checks of the JCP text, %d assertions' % (len(res), n))
    return 0


if __name__ == '__main__':
    sys.exit(main())
