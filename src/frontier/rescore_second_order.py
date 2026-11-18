# -*- coding: utf-8 -*-
r"""rescore_second_order.py -- P2 (deposit 3.5.0, the manuscript prepared for JCP): every number that depends on the
upper branch of Theorem 1, before and after the adoption of the second-order branch (rules D1-D23 of that change,
docs/KNOWN_DISCREPANCIES.md section 36; the D17-D23 additions of 2026-11-16 are section 18 and v2.json).

NO SIMULATION.  Every quantity is arithmetic on fields already stored in data/ (and, for the literals the
guards record, src/jcp/number_sources.json).  No eigen-solver, no sampling, no exact diagonalization, no fit
or optimisation here (the fitted-rule comparison of D20 is the separate copy src/frontier/null_ws7_second_order.py
of src/frontier/null_ws7.py, run after this one: it reads FRstar_new from cells.json),
no network, no file written outside its output folder.  Standard library and numpy only; deterministic; about one
second.

Deposit 3.5.0: this is the script of the JCP package (written 2026-11-15/16), with its paths taken from its place in
the repository: it reads data/ and src/jcp/number_sources.json and writes data/c3_frontier/second_order/ (--out DIR
writes elsewhere), and it records its inputs as repository-relative paths.  Nothing else changed: the outputs it
writes are those of the package, number for number (docs/C3_FRONTIER.md section 8 says how that was checked).

Notation, relative to ||phi||^2 (Theorem 1, D1):
    w  = captured weight w_P                 a  = sqrt(1-w)          x = Lambda^_P / eta  (Lambda / (||phi|| eta))
    B1 = (1+sqrt w)(a+x)                     first-order (leakage) branch of the present text
    B2 = (a+x)^2                             second-order branch for real H, P, phi_P  (Lambda~ = Lambda)
    T  = 1+w                                 trivial branch
    old certificate = min(T, B1)             new certificate (real probes, D6) = min(T, B2)
Complex momentum-resolved probes (cert_akw.json rows with k not in {0, pi}): first-order only (D8).
The Krylov object P = P_K inside S is real whenever H_S and phi_S are (all site-seeded rows), so D6 applies to it.

Usage (from the repository root):
        python src/frontier/rescore_second_order.py [data dir] [number_sources.json] [--out DIR]
Writes: grid_rows.json, cells.json, published.json, five_percent.json, supplement.json, pooled.json, checks.json,
        summary.json, v2.json into data/c3_frontier/second_order/ (see the end of the file).  summary.json and
        v2.json record their run times (runtime_s); they are the two outputs that differ from run to run.
"""
import os
import sys
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import glob
import json
import math
import re
import statistics
import time
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_UP

import numpy as np

T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))                    # <repo>/src/frontier
REPO = os.path.normpath(os.path.join(HERE, '..', '..'))              # <repo>
_ARGS = sys.argv[1:]
OUT = os.path.join(REPO, 'data', 'c3_frontier', 'second_order')
if '--out' in _ARGS:
    _i = _ARGS.index('--out')
    OUT = os.path.abspath(_ARGS[_i + 1])
    del _ARGS[_i:_i + 2]
DATA = os.path.abspath(_ARGS[0]) if len(_ARGS) > 0 else os.path.join(REPO, 'data')
NSRC = os.path.abspath(_ARGS[1]) if len(_ARGS) > 1 else os.path.join(REPO, 'src', 'jcp', 'number_sources.json')
os.makedirs(OUT, exist_ok=True)


def _rel(p):
    """a path as it is recorded in summary.json: repository-relative when it lies in the repository"""
    try:
        r = os.path.relpath(p, REPO)
    except ValueError:                                               # another drive
        return p
    return p if r.startswith('..') else r.replace('\\', '/')


def J(*p):
    with open(os.path.join(DATA, *p), encoding='utf-8') as fh:
        return json.load(fh)


def T(*p):
    with open(os.path.join(DATA, *p), encoding='utf-8') as fh:
        return fh.read()


# ----------------------------------------------------------------------------------------------- formulas
def B1(w, x):
    """first-order branch, relative: (1+sqrt w)(sqrt(1-w)+x)"""
    return (1.0 + math.sqrt(w)) * (math.sqrt(max(1.0 - w, 0.0)) + x)


def B2(w, x):
    """second-order branch, relative, real probes: (sqrt(1-w)+x)^2"""
    return (math.sqrt(max(1.0 - w, 0.0)) + x) ** 2


def tau_old(w):
    """first-order non-vacuity threshold, present Eq. (tau): B1 < 1+w  <=>  x < tau_old(w)"""
    return (1.0 + w) / (1.0 + math.sqrt(w)) - math.sqrt(max(1.0 - w, 0.0))


def tau_new(w):
    """second-order threshold, real probes: B2 < 1+w  <=>  x < sqrt(1+w) - sqrt(1-w)"""
    return math.sqrt(1.0 + w) - math.sqrt(max(1.0 - w, 0.0))


def x_from_B1(B, w):
    """invert the stored first-order branch: x = B/(1+sqrt w) - sqrt(1-w)"""
    return B / (1.0 + math.sqrt(w)) - math.sqrt(max(1.0 - w, 0.0))


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# ----------------------------------------------------------------------------------------------- rounding
_MODE = {'up': ROUND_CEILING, 'down': ROUND_FLOOR, 'near': ROUND_HALF_UP}


def rnd(v, n, mode='near', kind='sf'):
    """v rounded to n significant figures (kind='sf') or n decimals (kind='dp'), as a string.
    mode 'up' for a ceiling ('at most', '<=', an upper end), 'down' for a floor ('at least', a lower end),
    'near' for a plain value (rounding rule 7b, docs/KNOWN_DISCREPANCIES.md section 36)."""
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return None
    d = Decimal(repr(float(v)))
    if d == 0:
        return '0'
    e = (d.adjusted() - (n - 1)) if kind == 'sf' else -n
    return format(d.quantize(Decimal(1).scaleb(e), rounding=_MODE[mode]), 'f')


def sci(v, n, mode='near'):
    """'m x10^e' with an n-figure mantissa, rounded in the given direction."""
    if v is None:
        return None
    d = Decimal(repr(float(v)))
    e = d.adjusted()
    q = Decimal(1).scaleb(-(n - 1))
    m = d.scaleb(-e).quantize(q, rounding=_MODE[mode])
    if m >= 10:
        m = (m / 10).quantize(q, rounding=_MODE[mode])
        e += 1
    return '%sx10^%d' % (format(m, 'f'), e)


# ----------------------------------------------------------------------------------------------- fsum.py logic
def frontier(Bs, FRs, th):
    """src/frontier/fsum.py frontier(), verbatim logic: first FR at which B < th, log-linear in FR of
    log(B/th) between the two bracketing grid fractions."""
    order = np.argsort(FRs)
    FRs = np.asarray(FRs, dtype=float)[order]
    Bs = np.asarray(Bs, dtype=float)[order]
    th = np.asarray(th, dtype=float)[order]
    below = Bs < th
    if below.all():
        return (float(FRs[0]), float('nan'), float(FRs[0]), 'below-grid')
    if not below.any():
        return (float('nan'), float(FRs[-1]), float('nan'), 'above-grid')
    i = int(np.argmax(below))
    if i == 0:
        return (float(FRs[0]), float('nan'), float(FRs[0]), 'below-grid')
    x0, x1 = FRs[i - 1], FRs[i]
    y0 = math.log(Bs[i - 1] / th[i - 1])
    y1 = math.log(Bs[i] / th[i])
    return (float(x0 + (x1 - x0) * (0.0 - y0) / (y1 - y0)), float(x0), float(x1), 'bracketed')


ETAS = (0.10, 0.15, 0.18, 0.25)
LS = (6, 8, 10, 12)
EK = lambda e: '%g' % e                                     # key of ref_drift dicts: '0.1', '0.15', ...

# =============================================================================================== 0. the grid
# Stored fields used (c3_frontier/C3_grid.json, /grid[i]):
#   w_S            -> w
#   Lam_over_eta   -> x = Lambda^_K/eta (Lambda_K at the depth run; NOT Lam_out_over_eta, the n->inf 'limit'
#                     variant; verified below: B1(w_S, Lam_over_eta) reproduces the stored bound_leak)
#   bound_trivial  -> T (checked against 1+w_S)
#   bound_leak, cert, slack, slack_leakbranch, nontrivial  -> the present (first-order) values, for the check
#   relL1          -> true relative L1 error ||A - A_K||_1 / ||phi||^2 (whole line)
#   conv           -> CONVERGED / MARGINAL / NOT-CONVERGED (fsum.py depth test)
#   FR (grid fraction, used by fsum.py frontier) and FR_real (= |S|/D, used by verify_jcp.py interpolation)
G = J('c3_frontier', 'C3_grid.json')
GRID, FRONT = G['grid'], G['frontier']
rows = []
dev = dict(B1_vs_bound_leak=0.0, T_vs_1pw=0.0, cert_vs_min=0.0, slack_lb=0.0, x_vs_out=0.0)
for i, r in enumerate(GRID):
    w, x, Tv = r['w_S'], r['Lam_over_eta'], r['bound_trivial']
    b1, b2 = B1(w, x), B2(w, x)
    dev['B1_vs_bound_leak'] = max(dev['B1_vs_bound_leak'], rel(b1, r['bound_leak']))
    dev['T_vs_1pw'] = max(dev['T_vs_1pw'], rel(Tv, 1.0 + w))
    dev['cert_vs_min'] = max(dev['cert_vs_min'], rel(r['cert'], min(Tv, r['bound_leak'])))
    if r['relL1']:
        dev['slack_lb'] = max(dev['slack_lb'], rel(r['slack_leakbranch'], r['bound_leak'] / r['relL1']))
    dev['x_vs_out'] = max(dev['x_vs_out'], rel(x, r['Lam_out_over_eta']))
    c1, c2 = min(Tv, b1), min(Tv, b2)
    rows.append(dict(
        idx=i, L=r['L'], eta=r['eta'], FR=r['FR'], FR_real=r['FR_real'], S=r['S'], D=r['D'], nl=r['nl'],
        w=w, one_minus_w=1.0 - w, a=math.sqrt(max(1.0 - w, 0.0)), x=x,
        B1=b1, B2=b2, trivial=Tv, cert_old=c1, cert_new=c2,
        B1_over_T=b1 / Tv, B2_over_T=b2 / Tv,
        nontrivial_old=b1 < Tv, nontrivial_new=b2 < Tv, B2_le_B1=b2 <= b1,
        relL1=r['relL1'], conv=r['conv'], in_over_out=r['in_over_out'],
        B2_over_true=(b2 / r['relL1']) if r['relL1'] else None,
        cert_new_over_true=(c2 / r['relL1']) if r['relL1'] else None,
        cert_old_over_true=(c1 / r['relL1']) if r['relL1'] else None,
        src='c3_frontier/C3_grid.json /grid[%d]' % i))
assert dev['B1_vs_bound_leak'] < 1e-12 and dev['T_vs_1pw'] < 1e-12 and dev['cert_vs_min'] < 1e-12
# Lemma (second_order_branch.tex): wherever B1 < T, B2 <= B1.
lemma_viol = [r['idx'] for r in rows if r['nontrivial_old'] and not r['B2_le_B1']]


def cell(L, eta, key='FR'):
    return sorted([r for r in rows if r['L'] == L and abs(r['eta'] - eta) < 1e-9], key=lambda r: r[key])


# grid point files: their own reference drift (an alternative source for the drift margin, see cells below)
PT = {}
for f in sorted(glob.glob(os.path.join(DATA, 'c3_frontier', 'pt_L*_nl*.json'))):
    if f.endswith('_PUB.json'):
        continue
    d = json.load(open(f, encoding='utf-8'))
    d = d[0] if isinstance(d, list) else d
    k = (d['meta']['L'], round(d['FR'], 4))
    if k in PT and PT[k]['nl'] >= d['nl']:          # fsum.py: keep the deepest file per (L, FR)
        continue
    PT[k] = dict(nl=d['nl'], drift=d['meta']['ref_drift'], file='c3_frontier/' + os.path.basename(f))
# published-subspace files: reference drift per L (the convention of src/verify_jcp.py check_threshold_error)
DRIFT_PUB, PUBFILE = {}, {}
for f in sorted(glob.glob(os.path.join(DATA, 'c3_frontier', 'published_fraction', 'pt_L*_PUB.json'))):
    d = json.load(open(f, encoding='utf-8'))
    d = d[0] if isinstance(d, list) else d
    DRIFT_PUB[d['meta']['L']] = d['meta']['ref_drift']
    PUBFILE[d['meta']['L']] = 'c3_frontier/published_fraction/' + os.path.basename(f)


def interp_V(frs, cl, ykey='relL1', log=True):
    """verify_jcp.py convention: every row of the cell sorted by FR_real, the first bracketing pair, value at
    FR* (an FR-coordinate number) interpolated in FR_real; log-linear for the error."""
    s = sorted(cl, key=lambda r: r['FR_real'])
    for a, b in zip(s, s[1:]):
        if a['FR_real'] <= frs <= b['FR_real']:
            t = (frs - a['FR_real']) / (b['FR_real'] - a['FR_real'])
            if log:
                return math.exp((1 - t) * math.log(a[ykey]) + t * math.log(b[ykey])), (a['idx'], b['idx'])
            return (1 - t) * a[ykey] + t * b[ykey], (a['idx'], b['idx'])
    return None, None


def interp_C(frs, cl, ykey='relL1', log=True):
    """calib_oos.py / null_ws7.py convention: CONVERGED rows with relL1>0, sorted by FR, np.interp in FR."""
    s = [r for r in cl if r['conv'] == 'CONVERGED' and r['relL1'] and r['relL1'] > 0]
    xs = np.array([r['FR'] for r in s])
    if frs < xs[0] or frs > xs[-1]:
        return None
    ys = np.array([r[ykey] for r in s])
    if log:
        return float(math.exp(np.interp(frs, xs, np.log(ys))))
    return float(np.interp(frs, xs, ys))


def first_below(cl, key, lev):
    """first grid fraction (FR order) at which cl[key] < lev; bracket, convergence of its ends, the log-linear
    interpolation in FR, and whether the sequence is monotone after it."""
    for i, r in enumerate(cl):
        if r[key] < lev:
            break
    else:
        return dict(reached=False, above_FR=cl[-1]['FR'], above_FR_real=cl[-1]['FR_real'],
                    value_at_last=cl[-1][key], conv_last=cl[-1]['conv'], printed='>%.2f' % cl[-1]['FR'])
    later_above = [q['FR'] for q in cl[i:] if q[key] >= lev]
    if i == 0:
        return dict(reached=True, below_grid=True, FR=cl[0]['FR'], printed='<=%.2f' % cl[0]['FR'],
                    later_above=later_above)
    a, b = cl[i - 1], cl[i]
    t = math.log(a[key] / lev) / math.log(a[key] / b[key])
    ca, cb = a['conv'], b['conv']
    mark = 'n.c.' if 'NOT-CONVERGED' in (ca, cb) else ('' if ca == cb == 'CONVERGED' else '*')
    return dict(reached=True, lo=a['FR'], hi=b['FR'], lo_FR_real=a['FR_real'], hi_FR_real=b['FR_real'],
                lo_value=a[key], hi_value=b[key], conv_lo=ca, conv_hi=cb, mark=mark,
                FR_interp_loglin=a['FR'] + t * (b['FR'] - a['FR']),
                FR_real_interp_loglin=a['FR_real'] + t * (b['FR_real'] - a['FR_real']),
                rows=(a['idx'], b['idx']), later_above=later_above,
                printed=('n.c.' if mark == 'n.c.' else '[%.2f,%.2f]%s' % (a['FR'], b['FR'], mark)))


# =============================================================================================== 1. the cells
CELLS = []
for j, f in enumerate(FRONT):
    L, eta = f['L'], f['eta']
    cl = cell(L, eta)
    FRs = [r['FR'] for r in cl]
    Ts = [r['trivial'] for r in cl]
    c1 = frontier([r['B1'] for r in cl], FRs, Ts)
    c2 = frontier([r['B2'] for r in cl], FRs, Ts)
    assert abs(c1[0] - f['FR_star']) < 1e-12 and c1[1] == f['bracket_lo'] and c1[2] == f['bracket_hi']
    lo2 = [r for r in cl if abs(r['FR'] - c2[1]) < 1e-12][0]
    hi2 = [r for r in cl if abs(r['FR'] - c2[2]) < 1e-12][0]
    lo1 = [r for r in cl if abs(r['FR'] - c1[1]) < 1e-12][0]
    hi1 = [r for r in cl if abs(r['FR'] - c1[2]) < 1e-12][0]
    out = dict(L=L, eta=eta, src='c3_frontier/C3_grid.json /frontier[%d]' % j,
               FRstar_old=c1[0], FRstar_old_stored=f['FR_star'], br_old=[c1[1], c1[2]],
               FRstar_new=c2[0], br_new=[c2[1], c2[2]], how_new=c2[3],
               br_new_FR_real=[lo2['FR_real'], hi2['FR_real']], br_new_conv=[lo2['conv'], hi2['conv']],
               br_old_conv=[lo1['conv'], hi1['conv']], br_new_rows=[lo2['idx'], hi2['idx']])
    for lab, frs in (('old', c1[0]), ('new', c2[0])):
        eV, pr = interp_V(frs, cl)
        eC = interp_C(frs, cl)
        wV, _ = interp_V(frs, cl, 'w', log=False)
        wC = interp_C(frs, cl, 'w', log=False)
        omwV, _ = interp_V(frs, cl, 'one_minus_w', log=True)
        omwC = interp_C(frs, cl, 'one_minus_w', log=True)
        out['err_%s' % lab] = eV                                   # primary: verify_jcp.py convention
        out['err_%s_rows' % lab] = pr
        out['err_%s_calib' % lab] = eC                             # calib_oos.py convention
        out['ratio_%s' % lab] = (1.0 + wC) / eV                    # bound (= 1+w there) over the true error
        out['ratio_%s_V' % lab] = (1.0 + wV) / eV
        out['omw_%s_linFRconv' % lab] = 1.0 - wC                   # 1-w at the crossing, four conventions
        out['omw_%s_linFRrealAll' % lab] = 1.0 - wV
        out['omw_%s_logFRconv' % lab] = omwC
        out['omw_%s_logFRrealAll' % lab] = omwV
        out['drift_PUB'] = DRIFT_PUB[L][EK(eta)]
        out['drift_margin_%s_PUB' % lab] = eV / DRIFT_PUB[L][EK(eta)]
        a, b = pr
        dg = max(PT[(L, round(rows[a]['FR'], 4))]['drift'][EK(eta)], PT[(L, round(rows[b]['FR'], 4))]['drift'][EK(eta)])
        out['drift_margin_%s_gridfile' % lab] = eV / dg
    # least loose certified grid point (smallest certificate/true error among rows where the branch beats T)
    for lab, bk, ck in (('old', 'B1', 'cert_old'), ('new', 'B2', 'cert_new')):
        nt = [r for r in cl if r[bk] < r['trivial'] and r['relL1']]
        ntc = [r for r in nt if r['conv'] == 'CONVERGED']
        out['least_loose_%s_all' % lab] = min((r[ck] / r['relL1'] for r in nt), default=None)
        out['least_loose_%s_conv' % lab] = min((r[ck] / r['relL1'] for r in ntc), default=None)
        if ntc:
            q = min(ntc, key=lambda r: r[ck] / r['relL1'])
            out['least_loose_%s_conv_row' % lab] = q['idx']
    # the fraction where the TRUE error first falls below 0.05 (stored FR_true['0.05']) and the branches there
    ft = f['FR_true']['0.05']
    sel2 = [r for r in cl if r['relL1'] is not None]
    chk = frontier([r['relL1'] for r in sel2], [r['FR'] for r in sel2], [0.05] * len(sel2))
    assert abs(chk[0] - ft['FR']) < 1e-12
    ra = [r for r in cl if abs(r['FR'] - ft['lo']) < 1e-12][0]
    rb = [r for r in cl if abs(r['FR'] - ft['hi']) < 1e-12][0]
    t = (ft['FR'] - ft['lo']) / (ft['hi'] - ft['lo'])
    w05 = (1 - t) * ra['w'] + t * rb['w']
    x05 = math.exp((1 - t) * math.log(ra['x']) + t * math.log(rb['x']))
    out.update(FR05=ft['FR'], FR05_br=[ft['lo'], ft['hi']], FR05_src=out['src'] + "/FR_true/0.05",
               B1_over_T_at_FR05=B1(w05, x05) / (1 + w05), B2_over_T_at_FR05=B2(w05, x05) / (1 + w05),
               cert_new_over_T_at_FR05=min(1.0, B2(w05, x05) / (1 + w05)),
               B1_over_T_upper_row=rb['B1_over_T'], B2_over_T_upper_row=rb['B2_over_T'], conv_upper_row=rb['conv'],
               B2_over_T_lower_row=ra['B2_over_T'], conv_lower_row=ra['conv'],
               FRstar_new_gt_FR05=c2[0] > ft['FR'], FRstar_old_gt_FR05=c1[0] > ft['FR'])
    # where the certificate itself first falls below 0.50, 0.20, 0.05 (Table S1 cols 6-7; the 5 % statements)
    for lab, ck in (('old', 'cert_old'), ('new', 'cert_new')):
        for lev in (0.5, 0.2, 0.05):
            out['cert_%s_below_%s' % (lab, lev)] = first_below(cl, ck, lev)
    CELLS.append(out)

byc = {(c['L'], c['eta']): c for c in CELLS}


def summ(v):
    v = sorted(v)
    return dict(n=len(v), min=v[0], max=v[-1], median=statistics.median(v), spread=v[-1] / v[0],
                geo_mid=math.sqrt(v[0] * v[-1]), half_spread=math.sqrt(v[-1] / v[0]),
                spread_extremes_dropped=v[-2] / v[1])


CROSS = {}
for lab in ('old', 'new'):
    e = [c['err_%s' % lab] for c in CELLS]
    eC = [c['err_%s_calib' % lab] for c in CELLS]
    s = summ(e)
    s['calib_convention'] = summ(eC)
    s['max_rel_diff_between_conventions'] = max(rel(a, b) for a, b in zip(e, eC))
    s['ratio_at_crossing'] = summ([c['ratio_%s' % lab] for c in CELLS])
    s['least_loose_all'] = summ([c['least_loose_%s_all' % lab] for c in CELLS])
    s['least_loose_conv'] = summ([c['least_loose_%s_conv' % lab] for c in CELLS])
    s['drift_margin_PUB'] = min(((c['drift_margin_%s_PUB' % lab], c['L'], c['eta']) for c in CELLS))
    s['drift_margin_gridfile'] = min(((c['drift_margin_%s_gridfile' % lab], c['L'], c['eta']) for c in CELLS))
    by = {(c['L'], c['eta']): c['err_%s' % lab] for c in CELLS}
    s['falls_with_L_every_eta'] = all(by[(L1, x)] > by[(L2, x)] for x in ETAS for L1, L2 in ((6, 8), (8, 10), (10, 12)))
    s['rises_with_eta_every_L'] = all(by[(L, e1)] < by[(L, e2)] for L in LS
                                      for e1, e2 in ((0.10, 0.15), (0.15, 0.18), (0.18, 0.25)))
    f612 = [by[(6, x)] / by[(12, x)] for x in ETAS]
    fe = [by[(L, 0.25)] / by[(L, 0.10)] for L in LS]
    s['L6_over_L12'] = dict(min=min(f612), max=max(f612), per_eta=f612)
    s['eta025_over_eta010'] = dict(min=min(fe), max=max(fe), per_L=fe)
    s['eta018_L6_to_L12'] = [by[(L, 0.18)] for L in LS]
    for conv in ('linFRconv', 'linFRrealAll', 'logFRconv', 'logFRrealAll'):
        v = [c['omw_%s_%s' % (lab, conv)] for c in CELLS]
        s['omw_' + conv] = summ(v)
    # the "extreme cell at each end" read as the cells extreme in the ERROR (an alternative reading of sm_frontier)
    order = sorted(CELLS, key=lambda c: c['err_%s' % lab])[1:-1]
    v = [c['omw_%s_linFRconv' % lab] for c in order]
    s['omw_linFRconv_drop_error_extremes'] = max(v) / min(v)
    s['FRstar_gt_FR05_cells'] = sum(c['FRstar_%s_gt_FR05' % lab] for c in CELLS)
    CROSS[lab] = s
# FR*(L) log-log local slopes at each eta (sm_frontier "varying by a factor 4.1 across the four sizes")
SLOPES = {}
for lab in ('old', 'new'):
    SLOPES[lab] = {}
    for eta in ETAS:
        fr = [byc[(L, eta)]['FRstar_%s' % lab] for L in LS]
        sl = [math.log(fr[k + 1] / fr[k]) / math.log(LS[k + 1] / LS[k]) for k in range(3)]
        SLOPES[lab][EK(eta)] = dict(FRstar=fr, slopes=sl, factor=max(abs(q) for q in sl) / min(abs(q) for q in sl))
# bracket widths of Table S1 columns 4 (true error 0.05) and 5 (FR*), and of FR* over all sixteen cells
WIDTHS = {}
for lab in ('old', 'new'):
    wd5 = [byc[(L, e)]['br_%s' % lab][1] - byc[(L, e)]['br_%s' % lab][0] for L in LS for e in ETAS]
    wd4 = [byc[(L, e)]['FR05_br'][1] - byc[(L, e)]['FR05_br'][0] for L in LS for e in (0.18, 0.15)]
    wdt = [byc[(L, e)]['br_%s' % lab][1] - byc[(L, e)]['br_%s' % lab][0] for L in LS for e in (0.18, 0.15)]
    WIDTHS[lab] = dict(FRstar_all16=[min(wd5), max(wd5)], table_col5=[min(wdt), max(wdt)],
                       table_col4=[min(wd4), max(wd4)], table_cols45=[min(wd4 + wdt), max(wd4 + wdt)])

# gap between what the reconstruction needs (FR05) and FR* (Table S1 caption: "the gap widens")
GAP = {}
for lab in ('old', 'new'):
    GAP[lab] = {EK(e): dict(difference=[byc[(L, e)]['FRstar_%s' % lab] - byc[(L, e)]['FR05'] for L in LS],
                            ratio=[byc[(L, e)]['FRstar_%s' % lab] / byc[(L, e)]['FR05'] for L in LS]) for e in (0.18, 0.15)}
    for e in (0.18, 0.15):
        g = GAP[lab][EK(e)]
        g['difference_increasing_in_L'] = all(a < b for a, b in zip(g['difference'], g['difference'][1:]))
        g['ratio_increasing_in_L'] = all(a < b for a, b in zip(g['ratio'], g['ratio'][1:]))

# ties at the cut (sec_5_app, sm_leakage: "move the bound by at most 7.8e-4 relative")
# Stored: c3_frontier/ties.json [blk]/points[]/f32_w_S, f64_w_S, f32_eta*/B, f64_eta*/B (first-order branch of the two
# tie-breakings); x of each is derived from its B and w.
TIES = dict(src='c3_frontier/ties.json [L block]/points[]/{f32,f64}_w_S, {f32,f64}_eta*/B', B1_max_rel=0.0, B2_max_rel=0.0,
            cert_new_max_rel=0.0)
for blk in J('c3_frontier', 'ties.json'):
    for p in blk['points']:
        for k in p:
            if k.startswith('f32_eta') and isinstance(p[k], dict):
                k64 = k.replace('f32', 'f64')
                w32, w64 = p['f32_w_S'], p['f64_w_S']
                b32, b64 = p[k]['B'], p[k64]['B']
                c32 = B2(w32, x_from_B1(b32, w32))
                c64 = B2(w64, x_from_B1(b64, w64))
                TIES['B1_max_rel'] = max(TIES['B1_max_rel'], abs(b32 - b64) / b64)
                TIES['B2_max_rel'] = max(TIES['B2_max_rel'], abs(c32 - c64) / c64)
                TIES['cert_new_max_rel'] = max(TIES['cert_new_max_rel'],
                                               abs(min(1 + w32, c32) - min(1 + w64, c64)) / min(1 + w64, c64))

# =============================================================================================== 2. published subspaces
# Stored fields used:
#   c3_frontier/published_fraction/PUB_rows.json [i]: w_S, B (first-order branch), triv, vac (= B/triv), relL1, conv
#       x is derived as B/(1+sqrt w)-sqrt(1-w) and cross-checked against the LamHat_over_eta of the deepest rows
#       of the file named in 'src' (published_fraction/pt_L*_PUB.json)
#   c3_frontier/published_fraction/L14_FR008.json [0]: w_S; rows[] at the deepest depth n=250: LamHat_over_eta,
#       bound_leak, bound_trivial, vacuity_factor, relL1, in_over_out
PUB = J('c3_frontier', 'published_fraction', 'PUB_rows.json')
L14 = J('c3_frontier', 'published_fraction', 'L14_FR008.json')[0]
deep14 = max(r['n'] for r in L14['rows'])
publ = []
xcheck = 0.0
for i, r in enumerate(PUB):
    w = r['w_S']
    x = x_from_B1(r['B'], w)
    pf = json.load(open(os.path.join(DATA, 'c3_frontier', 'published_fraction', r['src']), encoding='utf-8'))
    pf = pf[0] if isinstance(pf, list) else pf
    dr = [q for q in pf['rows'] if q['n'] == pf['nl'] and abs(q['eta'] - r['eta']) < 1e-12]
    if dr:
        xcheck = max(xcheck, rel(x, dr[0]['LamHat_over_eta']))
    publ.append(dict(L=r['L'], eta=r['eta'], S=r['S'], FR_real=r['FR_real'], w=w, x=x, B1=r['B'], B2=B2(w, x),
                     trivial=r['triv'], vac_stored=r['vac'], relL1=r['relL1'], conv=r['conv'],
                     src='c3_frontier/published_fraction/PUB_rows.json [%d] (w_S, B, triv, vac)' % i,
                     x_src='derived from B; = %s rows[n=%d,eta=%g]/LamHat_over_eta' % (r['src'], pf['nl'], r['eta'])))
for i, r in enumerate(L14['rows']):
    if r['n'] == deep14:
        w = L14['w_S']
        x = r['LamHat_over_eta']
        assert rel(B1(w, x), r['bound_leak']) < 1e-12
        publ.append(dict(L=14, eta=r['eta'], S=L14['S'], FR_real=L14['FR_real'], w=w, x=x, B1=r['bound_leak'],
                         B2=B2(w, x), trivial=r['bound_trivial'], vac_stored=r['vacuity_factor'],
                         relL1=r.get('relL1'), conv='indicator-only (in/out=%.4g)' % r['in_over_out'],
                         src='c3_frontier/published_fraction/L14_FR008.json [0]/rows[%d] (n=%d)' % (i, deep14),
                         x_src='c3_frontier/published_fraction/L14_FR008.json [0]/rows[%d]/LamHat_over_eta' % i))
for p in publ:
    p['B1_over_T'] = p['B1'] / p['trivial']
    p['B2_over_T'] = p['B2'] / p['trivial']
    p['smaller_branch'] = 'second-order' if p['B2'] < p['B1'] else 'first-order'
    p['vac_new_smaller_branch'] = min(p['B1'], p['B2']) / p['trivial']           # D7
    p['tau_old'] = tau_old(p['w'])
    p['tau_new'] = tau_new(p['w'])
    p['x_over_tau_old'] = p['x'] / p['tau_old']
    p['x_over_tau_new'] = p['x'] / p['tau_new']
    p['uncertified_new'] = p['B2'] >= p['trivial'] and p['B1'] >= p['trivial']
publ.sort(key=lambda p: (p['L'], p['eta']))
pubtab = {(p['L'], p['eta']): p for p in publ}
mono = dict(
    rising_with_L=all(pubtab[(L1, e)]['vac_new_smaller_branch'] < pubtab[(L2, e)]['vac_new_smaller_branch']
                      for e in ETAS for L1, L2 in ((6, 8), (8, 10), (10, 12), (12, 14))),
    falling_with_eta=all(pubtab[(L, e1)]['vac_new_smaller_branch'] > pubtab[(L, e2)]['vac_new_smaller_branch']
                         for L in (6, 8, 10, 12, 14) for e1, e2 in ((0.10, 0.15), (0.15, 0.18), (0.18, 0.25))))
mono_B2 = dict(
    rising_with_L=all(pubtab[(L1, e)]['B2_over_T'] < pubtab[(L2, e)]['B2_over_T']
                      for e in ETAS for L1, L2 in ((6, 8), (8, 10), (10, 12), (12, 14))),
    falling_with_eta=all(pubtab[(L, e1)]['B2_over_T'] > pubtab[(L, e2)]['B2_over_T']
                         for L in (6, 8, 10, 12, 14) for e1, e2 in ((0.10, 0.15), (0.15, 0.18), (0.18, 0.25))))
# L=12 published, eta=0.18: depth 250 -> 300 (sm_frontier "the bound by 1.4e-5")
pf12 = json.load(open(os.path.join(DATA, 'c3_frontier', 'published_fraction', 'pt_L12_FR0.17_nl300_PUB.json'), encoding='utf-8'))
pf12 = pf12[0] if isinstance(pf12, list) else pf12
r250 = [q for q in pf12['rows'] if q['n'] == 250 and abs(q['eta'] - 0.18) < 1e-12][0]
r300 = [q for q in pf12['rows'] if q['n'] == 300 and abs(q['eta'] - 0.18) < 1e-12][0]
w12 = pf12['w_S']
L12DEPTH = dict(src='c3_frontier/published_fraction/pt_L12_FR0.17_nl300_PUB.json rows n=250,300 eta=0.18',
                Lout_rel_change=rel(r300['Lambda_K_out'], r250['Lambda_K_out']),
                B1_rel_change=rel(B1(w12, r300['LamHat_over_eta']), B1(w12, r250['LamHat_over_eta'])),
                B1_stored_rel_change=rel(r300['bound_leak'], r250['bound_leak']),
                B2_rel_change=rel(B2(w12, r300['LamHat_over_eta']), B2(w12, r250['LamHat_over_eta'])),
                in_over_out_n300=r300['in_over_out'])
# L=14, eta=0.10: depth 200 -> 250 (Table III caption: 8.627 -> 8.613, 0.17 %)
l14d = {}
for n in (200, 250):
    q = [r for r in L14['rows'] if r['n'] == n and abs(r['eta'] - 0.1) < 1e-12][0]
    w = L14['w_S']
    l14d[n] = dict(B1_over_T=q['vacuity_factor'], B2_over_T=B2(w, q['LamHat_over_eta']) / q['bound_trivial'])
L14DEPTH = dict(src='c3_frontier/published_fraction/L14_FR008.json [0]/rows n=200,250 eta=0.10 (vacuity_factor, LamHat_over_eta)',
                depths=l14d, B1_fall_pct=100 * (1 - l14d[250]['B1_over_T'] / l14d[200]['B1_over_T']),
                B2_fall_pct=100 * (1 - l14d[250]['B2_over_T'] / l14d[200]['B2_over_T']),
                smaller_branch_both_depths='first-order' if all(v['B1_over_T'] < v['B2_over_T'] for v in l14d.values()) else 'mixed')

# =============================================================================================== 3. tau facts
WS = np.linspace(1e-6, 1.0, 200001)
TN = np.sqrt(1 + WS) - np.sqrt(1 - WS)
TO = (1 + WS) / (1 + np.sqrt(WS)) - np.sqrt(1 - WS)
w12pub = [p for p in publ if p['L'] == 12][0]['w']
TAU = dict(
    tau_old_1=tau_old(1.0), tau_new_1=tau_new(1.0),
    tau_old_L12pub=tau_old(w12pub), tau_new_L12pub=tau_new(w12pub), w_L12pub=w12pub,
    tau_old_0995=tau_old(0.995), tau_new_0995=tau_new(0.995),
    move_old_pct=100 * (1 - tau_old(0.995) / tau_old(1.0)), move_new_pct=100 * (1 - tau_new(0.995) / tau_new(1.0)),
    tau_new_increasing_on_grid=bool(np.all(np.diff(TN) > 0)), tau_new_positive_on_grid=bool(np.all(TN > 0)),
    tau_old_positive_everywhere_on_grid=bool(np.all(TO > 0)),
    tau_new_ge_tau_old_on_grid=bool(np.all(TN >= TO - 1e-15)),
    grid='w = linspace(1e-6, 1, 200001), arithmetic only')

# =============================================================================================== 4. Table S2 (tab:certified)
CERTPTS = [(6, 0.18, 0.97), (8, 0.18, 0.95), (10, 0.18, 0.90), (12, 0.18, 0.60), (12, 0.15, 0.60), (12, 0.25, 0.70)]
TS2 = []
for L, eta, fr in CERTPTS:
    q = [r for r in rows if r['L'] == L and abs(r['eta'] - eta) < 1e-9 and abs(r['FR'] - fr) < 1e-9][0]
    ptf = PT[(L, round(fr, 4))]
    drift = ptf['drift'][EK(eta)]
    cl = cell(L, eta)
    conv = [r for r in cl if r['conv'] == 'CONVERGED']
    best_old = min(conv, key=lambda r: r['cert_old'])
    best_new = min(conv, key=lambda r: r['cert_new'])
    TS2.append(dict(L=L, eta=eta, FR=fr, FR_real=q['FR_real'], S=q['S'], src=q['src'],
                    cert_old=q['cert_old'], cert_new=q['cert_new'],
                    cert_new_printed_up_3sf=rnd(q['cert_new'], 3, 'up'), cert_old_printed_up_3sf=rnd(q['cert_old'], 3, 'up'),
                    relL1=q['relL1'], drift_own_file=drift, drift_file=ptf['file'],
                    true_error_quotable=(q['relL1'] / drift) > 20, conv=q['conv'], in_over_out=q['in_over_out'],
                    gap_old=q['cert_old'] / q['relL1'], gap_new=q['cert_new'] / q['relL1'],
                    cap_over_cert_old=q['trivial'] / q['cert_old'], cap_over_cert_new=q['trivial'] / q['cert_new'],
                    still_tightest_converged_in_cell_new=best_new['idx'] == q['idx'],
                    still_tightest_converged_in_cell_old=best_old['idx'] == q['idx']))
gq = [t['gap_new'] for t in TS2 if t['true_error_quotable']]
go = [t['gap_old'] for t in TS2 if t['true_error_quotable']]
TS2_SUM = dict(gap_new_quotable=[min(gq), max(gq)], gap_old_quotable=[min(go), max(go)],
               gap_new_all=[min(t['gap_new'] for t in TS2), max(t['gap_new'] for t in TS2)],
               log10_gap_new=[math.log10(min(gq)), math.log10(max(gq))],
               log10_gap_old=[math.log10(min(go)), math.log10(max(go))])

# =============================================================================================== 5. Table S1 (tab:thresholds)
TS1 = []
for eta in (0.18, 0.15):
    for L in LS:
        c = byc[(L, eta)]
        TS1.append(dict(L=L, eta=eta, FR05_br=c['FR05_br'], FR05=c['FR05'],
                        FRstar_old='[%.2f,%.2f]' % tuple(c['br_old']), FRstar_new='[%.2f,%.2f]' % tuple(c['br_new']),
                        FRstar_new_conv=c['br_new_conv'],
                        cert_lt_050_old=c['cert_old_below_0.5']['printed'], cert_lt_050_new=c['cert_new_below_0.5']['printed'],
                        cert_lt_020_old=c['cert_old_below_0.2']['printed'], cert_lt_020_new=c['cert_new_below_0.2']['printed'],
                        cert_lt_005_old=c['cert_old_below_0.05']['printed'], cert_lt_005_new=c['cert_new_below_0.05']['printed']))

# =============================================================================================== 6. looseness
LOOSE = {}
for lab, fn in (('B1_branch', lambda r: r['B1'] / r['relL1']), ('B2_branch', lambda r: r['B2'] / r['relL1']),
                ('cert_old', lambda r: r['cert_old'] / r['relL1']), ('cert_new', lambda r: r['cert_new'] / r['relL1'])):
    conv = [r for r in rows if r['conv'] == 'CONVERGED' and r['relL1']]
    byP, byC = {}, {}
    for r in conv:
        byP.setdefault((r['L'], round(r['FR'], 6)), []).append(fn(r))
        byC.setdefault((r['L'], r['eta']), []).append(fn(r))
    q = {}
    for (L, fr), v in byP.items():
        if len(v) >= 2:
            q.setdefault(L, []).append(max(v) / min(v))
    cc = {k: max(v) / min(v) for k, v in byC.items()}
    allv = [fn(r) for r in conv]
    LOOSE[lab] = dict(fixedP_over_eta_max_L68=max(max(q[6]), max(q[8])), fixedP_over_eta_max_all=max(max(v) for v in q.values()),
                      fixedP_over_eta_min=min(min(v) for v in q.values()),
                      across_subspaces_of_cell=[min(cc.values()), max(cc.values())],
                      across_subspaces_argmin=min(cc, key=cc.get), across_subspaces_argmax=max(cc, key=cc.get),
                      per_cell={'L%d_eta%g' % k: v for k, v in sorted(cc.items())},
                      all_converged_rows=[min(allv), max(allv)])
ntB2 = [r['B2'] / r['relL1'] for r in rows if r['conv'] == 'CONVERGED' and r['relL1'] and r['nontrivial_new']]
LOOSE['B2_branch']['converged_nontrivial_rows'] = [min(ntB2), max(ntB2)]
assert LOOSE['B1_branch']['fixedP_over_eta_max_L68'] > 0

# =============================================================================================== 7. oracle ranking (z_rank_L6.json)
# Stored fields: [i] ranking, FR, w_S, Lamhat_over_eta (x), relR5 (first-order branch, checked), relR0 (true), trivial_rel
ZR = J('c3_frontier', 'adversarial', 'z_rank_L6.json')
zr = []
for i, r in enumerate(ZR):
    b1, b2 = B1(r['w_S'], r['Lamhat_over_eta']), B2(r['w_S'], r['Lamhat_over_eta'])
    assert rel(b1, r['relR5']) < 1e-9
    zr.append(dict(i=i, ranking=r['ranking'], FR=r['FR'], w=r['w_S'], x=r['Lamhat_over_eta'], B1=b1, B2=b2,
                   smaller=min(b1, b2), smaller_is=('B2' if b2 < b1 else 'B1'), trivial=r['trivial_rel'],
                   cert_new=min(r['trivial_rel'], b2), true=r['relR0'],
                   src='c3_frontier/adversarial/z_rank_L6.json [%d]' % i))
ORACLE = dict(rows=zr)
for key in ('B1', 'B2', 'smaller', 'cert_new'):
    born = {round(z['FR'], 4): z[key] for z in zr if z['ranking'] == 'repo'}
    orc = {round(z['FR'], 4): z[key] for z in zr if z['ranking'] == 'krylov'}
    fr50 = [f for f in sorted(born) if f >= 0.5 - 1e-9]
    ORACLE['oracle_worse_from_050_' + key] = all(orc[f] > born[f] for f in fr50) if key != 'cert_new' else \
        all(orc[f] >= born[f] for f in fr50)
    ORACLE['oracle_better_at_030_' + key] = orc[0.3] < born[0.3]
    ORACLE['pairs_' + key] = {('%.2f' % f): dict(born=born[f], oracle=orc[f]) for f in sorted(born)}

# =============================================================================================== 8. S* = supp(phi), w = 1
# Stored fields: support_witness.json /results[k]: L, w_S, frac; /table[j]: eta, relL1, Lambda_hat_over_eta (x),
#   bound_theoremB (first-order, checked), bound_trivial
SW = J('support_witness.json')
star = []
for k, res in enumerate(SW['results']):
    for j, t in enumerate(res['table']):
        w = min(res['w_S'], 1.0)
        b1, b2 = B1(w, t['Lambda_hat_over_eta']), B2(w, t['Lambda_hat_over_eta'])
        assert rel(b1, t['bound_theoremB']) < 1e-9
        star.append(dict(L=res['L'], eta=t['eta'], frac=res['frac'], w_S_stored=res['w_S'], relL1=t['relL1'],
                         x=t['Lambda_hat_over_eta'], B1=b1, B2=b2, trivial=t['bound_trivial'],
                         vacuous_new=b2 >= t['bound_trivial'], B2_over_true=b2 / t['relL1'],
                         src='support_witness.json /results[%d]/table[%d]' % (k, j)))

# =============================================================================================== 9. pooled sweep
# Stored fields (the five sources of src/recovered/build_n3.py):
#   thm1iii_violation/thm3_results.json [i]: LHS (error), relL1, wS, norm_phi2, Lam (absolute Lambda), eta,
#       RHS_new (first-order, checked), RHS_triv
#   thm1iii_violation/thm3_L8b.json [i]: LHS, relL1, wS, LB (= ||phi||^2 (1-w)), Lam, Lam_over_eta (ABSOLUTE
#       Lambda/eta here: eta = Lam/Lam_over_eta = 0.18), RHS_new
#   cert_stress.json, cert_akw.json, cert_teqsci.json /rows[i]: L1_true, relL1_true, w_S, norm_phi, norm_phiS,
#       norm_phiQ, Lambda (absolute), eta, upper_unclipped (first-order, checked), upper_trivial, lower, k_over_pi
POOL = []
for i, r in enumerate(J('thm1iii_violation', 'thm3_results.json')):
    n2 = r['norm_phi2']; nphi = math.sqrt(n2); w = r['wS']
    nQ, nS = nphi * math.sqrt(max(1 - w, 0)), nphi * math.sqrt(w)
    b1 = (nphi + nS) * (nQ + r['Lam'] / r['eta'])
    assert rel(b1, r['RHS_new']) < 1e-8
    POOL.append(dict(src='thm1iii_violation/thm3_results.json [%d]' % i, kind=r['kind'], err=r['LHS'], rel=r['relL1'],
                     triv=r['RHS_triv'], b1=r['RHS_new'], b2=(nQ + r['Lam'] / r['eta']) ** 2, w=w, eta=r['eta'],
                     frac=r['frac'], real=True, lam=r['Lam']))
for i, r in enumerate(J('thm1iii_violation', 'thm3_L8b.json')):
    w = r['wS']; n2 = r['LB'] / (1 - w); nQ = math.sqrt(n2) * math.sqrt(1 - w)
    eta = r['Lam'] / r['Lam_over_eta']
    POOL.append(dict(src='thm1iii_violation/thm3_L8b.json [%d]' % i, kind='L8b', err=r['LHS'], rel=r['relL1'],
                     triv=n2 * (1 + w), b1=r['RHS_new'], b2=(nQ + r['Lam_over_eta']) ** 2, w=w, eta=eta,
                     frac=r['FR'], real=True, lam=r['Lam']))
for fn in ('cert_stress.json', 'cert_akw.json', 'cert_teqsci.json'):
    for i, r in enumerate(J(fn)['rows']):
        nQ = r['norm_phiQ']
        b1 = (r['norm_phi'] + r['norm_phiS']) * (nQ + r['Lambda'] / r['eta'])
        assert rel(b1, r['upper_unclipped']) < 1e-8
        kpi = r.get('k_over_pi')
        cplx = fn == 'cert_akw.json' and not (abs(kpi) < 1e-9 or abs(kpi - 1.0) < 1e-9)
        POOL.append(dict(src='%s /rows[%d]' % (fn, i), kind=fn, err=r['L1_true'], rel=r['relL1_true'],
                         triv=r['upper_trivial'], b1=r['upper_unclipped'], b2=(nQ + r['Lambda'] / r['eta']) ** 2,
                         w=r['w_S'], eta=r['eta'], frac=r['frac'], real=not cplx, akw=fn == 'cert_akw.json', kpi=kpi,
                         lower=r['lower'], L=r.get('L'), branch=r.get('branch'), lam=r['Lambda'],
                         lower_saturated=r.get('lower_saturated')))
live = [p for p in POOL if p['rel'] > 1e-10]
for p in live:
    p['cert_old'] = min(p['triv'], p['b1'])
    p['cert_new'] = min(p['triv'], p['b2']) if p['real'] else p['cert_old']                         # D8 (k = 0, pi real)
    p['cert_new_akw_all_first'] = p['cert_old'] if p.get('akw') else min(p['triv'], p['b2'])         # D8 read literally
    p['ratio_old'] = p['err'] / p['cert_old']
    p['ratio_new'] = p['err'] / p['cert_new']
    p['ratio_new_akw_all_first'] = p['err'] / p['cert_new_akw_all_first']
real_live = [p for p in live if p['real']]
cplx_live = [p for p in live if not p['real']]
akw_live = [p for p in live if p.get('akw')]
pmin = min(live, key=lambda p: p['ratio_new'])
pmax = max(live, key=lambda p: p['ratio_new'])
POOLSUM = dict(
    n_pooled=len(POOL), n_live=len(live), n_real_live=len(real_live), n_complex_live=len(cplx_live),
    n_akw_live=len(akw_live), n_akw_real_k0_or_pi_live=sum(1 for p in akw_live if p['real']),
    leak_wins_old=sum(p['b1'] < p['triv'] for p in live),
    real_rows_B2_wins=sum(p['b2'] < p['triv'] for p in real_live),
    complex_rows_B1_wins=sum(p['b1'] < p['triv'] for p in cplx_live),
    real_rows_B1_wins=sum(p['b1'] < p['triv'] for p in real_live),
    akw_all_first_order=dict(nonakw_B2_wins=sum(p['b2'] < p['triv'] for p in live if not p.get('akw')),
                             akw_B1_wins=sum(p['b1'] < p['triv'] for p in akw_live)),
    complex_rows_B2_wins_if_Ltilde_eq_L_INFO_ONLY=sum(p['b2'] < p['triv'] for p in cplx_live),
    ratio_old=[min(p['ratio_old'] for p in live), max(p['ratio_old'] for p in live)],
    ratio_new=[pmin['ratio_new'], pmax['ratio_new']], ratio_new_min_at=pmin['src'], ratio_new_max_at=pmax['src'],
    ratio_new_akw_all_first=[min(p['ratio_new_akw_all_first'] for p in live), max(p['ratio_new_akw_all_first'] for p in live)],
    violations_new=sum(p['ratio_new'] > 1 + 1e-9 for p in live),
    violations_real_B2_alone=sum(p['err'] > p['b2'] * (1 + 1e-9) for p in real_live),
    rows_at_max_ratio=[p['src'] + ' w=%g' % p['w'] for p in live if p['ratio_new'] > 1 - 1e-9])
POOLSUM['nontrivial_new_total'] = POOLSUM['real_rows_B2_wins'] + POOLSUM['complex_rows_B1_wins']
POOLSUM['nontrivial_new_total_akw_all_first'] = (POOLSUM['akw_all_first_order']['nonakw_B2_wins']
                                                 + POOLSUM['akw_all_first_order']['akw_B1_wins'])
# the 378 deposited rows (cert_*.json): both sides with the new certificate
dep = [p for p in POOL if p['kind'] in ('cert_stress.json', 'cert_akw.json', 'cert_teqsci.json')]
for p in dep:
    p.setdefault('cert_new', min(p['triv'], p['b2']) if p['real'] else min(p['triv'], p['b1']))
dep_def = [p for p in dep if p['rel'] > 1e-10]
# tolerance 1e-7 relative on the lower side: 39 cert_stress rows are 'lower_saturated' (error = missed weight, Lambda = 0)
# and their stored quadrature of the true L1 lies below ||phi_Q||^2 by up to 1.6e-8 relative
DEP378 = dict(n=len(dep), lower_holds_tol1e7=sum(p['lower'] <= p['err'] * (1 + 1e-7) + 1e-15 for p in dep),
              lower_holds_tol1e9=sum(p['lower'] <= p['err'] * (1 + 1e-9) + 1e-15 for p in dep),
              lower_saturated_rows=sum(1 for p in dep if p.get('lower_saturated')),
              max_rel_shortfall_lower=max((p['lower'] / p['err'] - 1) for p in dep if p['err'] > 0),
              upper_new_holds_all=sum(p['err'] <= p['cert_new'] * (1 + 1e-9) + 1e-15 for p in dep),
              upper_new_holds_ratio_defined=sum(p['err'] <= p['cert_new'] * (1 + 1e-9) for p in dep if p['rel'] > 1e-10),
              upper_fails_listed=[(p['src'], p['err'], p['cert_new'], p['rel']) for p in dep
                                  if not p['err'] <= p['cert_new'] * (1 + 1e-9) + 1e-15],
              n_ratio_defined=len(dep_def), max_ratio_new=max(p['err'] / p['cert_new'] for p in dep_def),
              max_ratio_old=max(p['err'] / min(p['triv'], p['b1']) for p in dep_def))

# =============================================================================================== 10. A(k,omega), cert_akw.json
# Stored fields: rows[i] L, FR, branch, k_over_pi, w_S, Lambda_hat_over_eta, rel_upper (= min(1+w, B1)), slack
#   (= upper / L1_true), relL1_true, bound_nontrivial.  Complex probes: only the first-order values are evaluable (D8).
CA = J('cert_akw.json')['rows']
ak = []
for i, r in enumerate(CA):
    if r['L'] == 8 and abs(r['FR'] - 0.85) < 1e-9:
        b1 = min(1 + r['w_S'], B1(r['w_S'], r['Lambda_hat_over_eta']))
        assert rel(b1, r['rel_upper']) < 1e-9
        realk = abs(r['k_over_pi']) < 1e-9 or abs(r['k_over_pi'] - 1) < 1e-9
        ak.append(dict(i=i, branch=r['branch'], k_over_pi=r['k_over_pi'], rel_upper=r['rel_upper'], slack=r['slack'],
                       true=r['relL1_true'], nontrivial=r['bound_nontrivial'], real_probe=realk,
                       B2_rel=(min(1 + r['w_S'], B2(r['w_S'], r['Lambda_hat_over_eta'])) if realk else None)))
AKW = dict(src='cert_akw.json /rows[] with L=8, FR=0.85 (rel_upper, slack, relL1_true, k_over_pi)', n=len(ak),
           rel_upper=[min(a['rel_upper'] for a in ak), max(a['rel_upper'] for a in ak)],
           slack=[min(a['slack'] for a in ak), max(a['slack'] for a in ak)],
           all_nontrivial_first_order=all(a['nontrivial'] for a in ak),
           real_channels_B2=[(a['branch'], a['k_over_pi'], a['B2_rel'], a['B2_rel'] / a['true']) for a in ak if a['real_probe']],
           channels=ak)

# =============================================================================================== 11. k-block counts
KB = {}
for L in (6, 8):
    d = J('c3_frontier', 'adversarial', 'z_kblock_L%d.json' % L)
    rr = []
    for i, r in enumerate(d['rows']):
        w = min(r['w_S'], 1.0)
        x = max(x_from_B1(r['bound'], w), 0.0)
        rr.append(dict(k_index=r['k_index'], block_size=r['block_size'], first_nonvacuous_first_order=r['first_nonvacuous'],
                       B1_at_that_count=r['bound'], B2_at_that_count=B2(w, x)))
    KB[L] = dict(src='c3_frontier/adversarial/z_kblock_L%d.json /rows[]/first_nonvacuous, bound, w_S' % L,
                 counts=[min(q['first_nonvacuous_first_order'] for q in rr), max(q['first_nonvacuous_first_order'] for q in rr)],
                 rows=rr)
# "a factor 5-7 fewer than the site determinant basis requires" (site FR* at eta=0.18 times D over the block counts)
KFAC = {}
for lab in ('old', 'new'):
    f6 = byc[(6, 0.18)]['FRstar_%s' % lab] * 300
    f8 = byc[(8, 0.18)]['FRstar_%s' % lab] * 3920
    KFAC[lab] = dict(L6=[f6 / KB[6]['counts'][1], f6 / KB[6]['counts'][0]], L8=[f8 / KB[8]['counts'][1], f8 / KB[8]['counts'][0]])

# =============================================================================================== 12. Fig. 1 (L=8, eta=0.18)
# Stored: C3_grid.json rows (L=8, eta=0.18) and c3_frontier/adversarial/z_fine_L8.txt (|S|, w_S, relR0 = true error,
# relR5 = first-order branch, Lamhat/eta = x); support_witness.json L=8 eta=0.18.
zf = T('c3_frontier', 'adversarial', 'z_fine_L8.txt')
fine = []
for m in re.finditer(r"\|S\|=(\d+) FR=([\d.]+) eta=0\.18 \| w_S=([\d.e+-]+) relR0=([\d.e+-]+) relR5=([\d.e+-]+) "
                     r"Lamhat/eta=([\d.e+-]+)", zf):
    S_, w, r0, r5, x = int(m.group(1)), float(m.group(3)), float(m.group(4)), float(m.group(5)), float(m.group(6))
    assert rel(B1(w, x), r5) < 1e-5
    fine.append(dict(S=S_, FR=S_ / 3920.0, w=w, relL1=r0, B1=r5, B2=B2(w, x), trivial=1 + w, missed=1 - w,
                     B2_over_true=B2(w, x) / r0, src='c3_frontier/adversarial/z_fine_L8.txt |S|=%d' % S_))
c8 = cell(8, 0.18, 'FR_real')
fig = [dict(FR=r['FR_real'], relL1=r['relL1'], B1=r['B1'], B2=r['B2'], trivial=r['trivial'], missed=r['one_minus_w'],
            src=r['src']) for r in c8] + [dict(FR=q['FR'], relL1=q['relL1'], B1=q['B1'], B2=q['B2'], trivial=q['trivial'],
                                               missed=q['missed'], src=q['src']) for q in fine]
st8 = [s for s in star if s['L'] == 8 and abs(s['eta'] - 0.18) < 1e-12][0]
FIG1 = dict(rows=fig, star=st8, FRstar_old=byc[(8, 0.18)]['FRstar_old'], FRstar_new=byc[(8, 0.18)]['FRstar_new'],
            FR05=byc[(8, 0.18)]['FR05'], B2_first=fig[0]['B2'], B2_last=fig[-1]['B2'], B1_first=fig[0]['B1'],
            B1_last=fig[-1]['B1'], B2_max=max(q['B2'] for q in fig), ymax_present=80,
            cert_new_below_005=byc[(8, 0.18)]['cert_new_below_0.05'],
            fine_points_B2=[q['B2'] for q in fine],
            fine_points_still_locate_005=any(q['B2'] >= 0.05 for q in fine))

# =============================================================================================== 13. never below the true error
NB = []


def nb(source, key, true, b2, cert=None, extra='', lam_pos=True):
    NB.append(dict(source=source, key=key, true=true, B2=b2, ratio=b2 / true if true else None,
                   cert_ratio=(cert / true) if (cert is not None and true) else None, note=extra, leakage_nonzero=lam_pos))


for r in rows:
    if r['relL1']:
        nb('C3_grid.json', r['src'], r['relL1'], r['B2'], r['cert_new'], 'L=%d FR=%.2f eta=%.2f %s' % (r['L'], r['FR'], r['eta'], r['conv']),
           r['x'] > 0)
for f in sorted(glob.glob(os.path.join(DATA, 'c3_frontier', 'published_fraction', 'pt_L*_PUB.json'))):
    d = json.load(open(f, encoding='utf-8'))
    d = d[0] if isinstance(d, list) else d
    for i, q in enumerate(d['rows']):
        if q.get('relL1'):
            b2 = B2(d['w_S'], q['LamHat_over_eta'])
            nb(os.path.basename(f), 'rows[%d] n=%d eta=%g' % (i, q['n'], q['eta']), q['relL1'], b2, min(b2, q['bound_trivial']),
               '', q['LamHat_over_eta'] > 0)
for i, q in enumerate(L14['rows']):
    if q.get('relL1'):
        b2 = B2(L14['w_S'], q['LamHat_over_eta'])
        nb('L14_FR008.json', 'rows[%d] n=%d eta=%g' % (i, q['n'], q['eta']), q['relL1'], b2, min(b2, q['bound_trivial']),
           '', q['LamHat_over_eta'] > 0)
for q in fine:
    nb('z_fine_L8.txt', q['src'], q['relL1'], q['B2'], min(q['B2'], q['trivial']))
for z in zr:
    nb('z_rank_L6.json', z['src'], z['true'], z['B2'], z['cert_new'], z['ranking'], z['x'] > 0)
for s in star:
    nb('support_witness.json', s['src'], s['relL1'], s['B2'], min(s['B2'], s['trivial']), '', s['x'] > 0)
for p in POOL:
    if p['rel'] > 1e-10 and p['real']:
        nb('pool (real rows)', p['src'], p['err'], p['b2'], min(p['b2'], p['triv']), 'absolute L1, w=%.3g' % p['w'], p['lam'] > 0)
FULLSECTOR = []          # S = whole sector: A_S = A by algebra, both branches 0; a nonzero stored error is reference drift
for fn in ('cert_frontier_6-8.json', 'cert_frontier_10.json'):
    d = J(fn)
    for Lk, cv in d['curves'].items():
        for i, q in enumerate(cv['points']):
            if q.get('S_is_full_sector'):
                FULLSECTOR.append(dict(src='%s curves/%s/points[%d]' % (fn, Lk, i), relL1_true=q.get('relL1_true'),
                                       upper=q['upper'], mode=q['mode']))
                continue
            if q.get('relL1_true') is not None and q['relL1_true'] > 1e-10:
                b2 = (q['norm_phiQ'] + q['Lambda'] / q['eta']) ** 2
                assert rel((q['norm_phi'] + q['norm_phiS']) * (q['norm_phiQ'] + q['Lambda'] / q['eta']), q['upper_unclipped']) < 1e-8
                nb(fn, 'curves/%s/points[%d] frac=%g' % (Lk, i, q['frac']), q['L1_true'], b2, min(b2, q['upper_trivial']), 'absolute L1',
                   q['Lambda'] > 0)
for Lf in (4, 6, 8):
    d = J('c3_frontier', 'ladder_L%d.json' % Lf)
    for ci, c in enumerate(d['cases']):
        for ri, q in enumerate(c['rows']):
            if q['R0'] > 1e-10 * d['norm_phi2']:
                assert rel((math.sqrt(c['norm_phi2']) + c['norm_phiS']) * (c['norm_phiQ'] + q['Lam_over_eta']), q['R5']) < 1e-6
                b2 = (c['norm_phiQ'] + q['Lam_over_eta']) ** 2          # Lam_over_eta is ABSOLUTE Lambda/eta here
                nb('ladder_L%d.json' % Lf, 'cases[%d] (%s) rows[%d] eta=%g' % (ci, c['tag'], ri, q['eta']), q['R0'], b2,
                   min(b2, q['trivial']), 'absolute L1', q['Lam'] > 0)
info_cplx = [p['err'] / p['b2'] for p in POOL if p['rel'] > 1e-10 and not p['real']]
nbv = [q for q in NB if q['ratio'] is not None]
NBSUM = dict(n_rows=len(nbv), violations=sum(q['ratio'] < 1 - 1e-9 for q in nbv),
             smallest=min(nbv, key=lambda q: q['ratio']),
             smallest_cert_ratio=min((q for q in nbv if q['cert_ratio'] is not None), key=lambda q: q['cert_ratio']),
             per_source={})
for q in nbv:
    s = NBSUM['per_source'].setdefault(q['source'], dict(n=0, min_ratio=float('inf'), at=None))
    s['n'] += 1
    if q['ratio'] < s['min_ratio']:
        s['min_ratio'], s['at'] = q['ratio'], q['key'] + ' ' + q['note']
NBSUM['complex_akw_rows_INFO_ONLY_err_over_B2_if_Ltilde_eq_L'] = max(info_cplx) if info_cplx else None
NBSUM['excluded_full_sector_rows'] = FULLSECTOR
lp = [q for q in nbv if q['leakage_nonzero']]
z0 = [q for q in nbv if not q['leakage_nonzero']]
NBSUM['n_rows_zero_leakage'] = len(z0)
NBSUM['zero_leakage_ratio_range'] = [min(q['ratio'] for q in z0), max(q['ratio'] for q in z0)] if z0 else None
NBSUM['smallest_with_nonzero_leakage'] = min(lp, key=lambda q: q['ratio'])
NBSUM['smallest_with_nonzero_leakage_per_source'] = {}
for q in lp:
    s = NBSUM['smallest_with_nonzero_leakage_per_source'].setdefault(q['source'], dict(n=0, min_ratio=float('inf'), at=None))
    s['n'] += 1
    if q['ratio'] < s['min_ratio']:
        s['min_ratio'], s['at'] = q['ratio'], q['key'] + ' ' + q['note']
grid_nt = [r for r in rows if r['relL1'] and r['nontrivial_new']]
NBSUM['grid_B2_over_true_where_nontrivial'] = [min(r['B2_over_true'] for r in grid_nt), max(r['B2_over_true'] for r in grid_nt)]
NBSUM['grid_min_cert_new_over_true'] = min((r['cert_new_over_true'], r['src']) for r in rows if r['relL1'])

# =============================================================================================== 14. D13 re-verification
D13 = dict(
    a_cert_no_better_than_trivial_at_FR05_every_cell=all(c['B2_over_T_at_FR05'] > 1 and c['B1_over_T_at_FR05'] > 1 for c in CELLS),
    a_min_B2_over_T_at_FR05=min(c['B2_over_T_at_FR05'] for c in CELLS),
    a_max_B2_over_T_at_FR05=max(c['B2_over_T_at_FR05'] for c in CELLS),
    a_min_B1_over_T_at_FR05=min(c['B1_over_T_at_FR05'] for c in CELLS),
    a_upper_bracket_row_B2_over_T_min=min(c['B2_over_T_upper_row'] for c in CELLS),
    a_FRstar_new_above_FR05_cells=CROSS['new']['FRstar_gt_FR05_cells'],
    b_every_published_subspace_uncertified=all(p['uncertified_new'] for p in publ),
    b_min_smaller_branch_over_T=min(p['vac_new_smaller_branch'] for p in publ),
    b_min_B2_over_T=min(p['B2_over_T'] for p in publ),
    c_Sstar_eta018_vacuous=all(s['vacuous_new'] for s in star if abs(s['eta'] - 0.18) < 1e-12),
    c_Sstar_all_eta_vacuous=all(s['vacuous_new'] for s in star),
    c_Sstar_nonvacuous_rows=[(s['L'], s['eta'], s['B2']) for s in star if not s['vacuous_new']],
    c_Sstar_error_range=[min(s['relL1'] for s in star if abs(s['eta'] - 0.18) < 1e-12),
                         max(s['relL1'] for s in star if abs(s['eta'] - 0.18) < 1e-12)],
    d_falls_with_L=CROSS['new']['falls_with_L_every_eta'], d_rises_with_eta=CROSS['new']['rises_with_eta_every_L'],
    lemma_B2_le_B1_where_B1_nontrivial_violations=lemma_viol)

# =============================================================================================== 15. present literals (guards)
NS = json.load(open(NSRC, encoding='utf-8'))
WANT = [('sm_frontier.tex', '8.544'), ('sm_frontier.tex', '17.21'), ('sm_frontier.tex', '11.18'), ('sm_frontier.tex', '9.21'),
        ('sm_frontier.tex', '6.50'), ('sm_frontier.tex', '0.9625'), ('sm_frontier.tex', '0.9280'), ('sm_frontier.tex', '0.296'),
        ('sm_frontier.tex', '0.067'), ('sm_frontier.tex', '0.503'), ('sm_frontier.tex', '0.220'), ('sm_frontier.tex', '0.137'),
        ('sm_frontier.tex', '0.718'), ('sm_frontier.tex', '0.878'), ('sm_frontier.tex', '0.219'), ('sm_frontier.tex', '0.060'),
        ('sm_frontier.tex', '0.039'), ('sm_frontier.tex', '3830'), ('sm_frontier.tex', '3855'), ('sm_frontier.tex', '21.76'),
        ('sm_frontier.tex', '24.52'), ('sm_frontier.tex', '12.97'), ('sm_frontier.tex', '4.27'), ('sm_frontier.tex', '5.82'),
        ('sm_frontier.tex', '0.911'), ('sm_frontier.tex', '0.02'), ('sm_frontier.tex', '455'), ('sm_frontier.tex', '483'),
        ('sec_5_body.tex', '0.9625'), ('sec_5_body.tex', '0.718'), ('sec_5_body.tex', '2.8'), ('sec_5_body.tex', '0.219'),
        ('sec_5_body.tex', '8.627'), ('sec_5_body.tex', '8.613'), ('sec_5_body.tex', '0.17'), ('sec_5_body.tex', '4.236'),
        ('sec_5_body.tex', '50'), ('sec_5_body.tex', '455'), ('sec_5_body.tex', '483'), ('sec_10.tex', '455'), ('sec_10.tex', '483'),
        ('sec_3_app.tex', '0.9999999997'), ('figs/caption_thm1iii_violation.tex', '0.9999999997')]
REC = []
for tex, lit in WANT:
    for rr in NS:
        if rr['tex'] == tex and rr['lit'] == lit:
            REC.append(dict(tex=tex, lit=lit, line=rr['line'], file=rr['file'], key=rr['key'], value=rr['value'],
                            expr=rr.get('expr'), round=rr.get('round'), note=rr.get('note')))
# a guard record that points to a field other than the printed one (found 2026-11-15): sm_frontier 9.21
l14r16 = L14['rows'][16]
l14r38 = L14['rows'][38]
MISSRC = dict(printed='9.21 (sm_frontier.tex, "leakage branch of 17.21, 11.18, 9.21 and 6.50", the eta=0.18 value)',
              record='L14_FR008.json [0]/rows[16]/LamHat_over_eta = %r (n=%d, eta=%g)' % (l14r16['LamHat_over_eta'], l14r16['n'], l14r16['eta']),
              actual_quantity='L14_FR008.json [0]/rows[38]/bound_leak = %r (n=%d, eta=%g)' % (l14r38['bound_leak'], l14r38['n'], l14r38['eta']))

# =============================================================================================== 16. printed candidates
c_old, c_new = CROSS['old'], CROSS['new']
PRINT = dict(
    crossing_error_range_new=dict(full=[c_new['min'], c_new['max']],
                                  outward_2sf=[sci(c_new['min'], 2, 'down'), sci(c_new['max'], 2, 'up')],
                                  nearest_2sf=[sci(c_new['min'], 2), sci(c_new['max'], 2)]),
    crossing_error_range_old=dict(full=[c_old['min'], c_old['max']], outward_2sf=[sci(c_old['min'], 2, 'down'), sci(c_old['max'], 2, 'up')],
                                  nearest_2sf=[sci(c_old['min'], 2), sci(c_old['max'], 2)]),
    crossing_median_new=(c_new['median'], sci(c_new['median'], 2)),
    crossing_spread_new=(c_new['spread'], rnd(c_new['spread'], 2)),
    crossing_geo_mid_new=(c_new['geo_mid'], round(c_new['geo_mid'] * 1e3), c_new['half_spread']),
    ratio_at_crossing_new=dict(full=[c_new['ratio_at_crossing']['min'], c_new['ratio_at_crossing']['max']],
                               outward_2sf=[sci(c_new['ratio_at_crossing']['min'], 2, 'down'), sci(c_new['ratio_at_crossing']['max'], 2, 'up')],
                               nearest_2sf=[sci(c_new['ratio_at_crossing']['min'], 2), sci(c_new['ratio_at_crossing']['max'], 2)]),
    ratio_at_crossing_old=dict(full=[c_old['ratio_at_crossing']['min'], c_old['ratio_at_crossing']['max']]),
    least_loose_conv_new=dict(full=[c_new['least_loose_conv']['min'], c_new['least_loose_conv']['max']],
                              outward_2sf=[sci(c_new['least_loose_conv']['min'], 2, 'down'), sci(c_new['least_loose_conv']['max'], 2, 'up')]),
    least_loose_all_new=dict(full=[c_new['least_loose_all']['min'], c_new['least_loose_all']['max']]),
    drift_margin_new_floor=(c_new['drift_margin_PUB'], math.floor(c_new['drift_margin_PUB'][0])),
    drift_margin_old_floor=(c_old['drift_margin_PUB'], math.floor(c_old['drift_margin_PUB'][0])),
    eta018_new=[sci(v, 3) for v in c_new['eta018_L6_to_L12']],
    L6_over_L12_new=[rnd(c_new['L6_over_L12']['min'], 2), rnd(c_new['L6_over_L12']['max'], 2)],
    eta_factor_new=[rnd(c_new['eta025_over_eta010']['min'], 2), rnd(c_new['eta025_over_eta010']['max'], 2)],
    fixedP_new_up1=[rnd(LOOSE['B2_branch']['fixedP_over_eta_max_L68'], 1, 'up', 'dp'), rnd(LOOSE['B2_branch']['fixedP_over_eta_max_all'], 1, 'up', 'dp')],
    fixedP_old_up1=[rnd(LOOSE['B1_branch']['fixedP_over_eta_max_L68'], 1, 'up', 'dp'), rnd(LOOSE['B1_branch']['fixedP_over_eta_max_all'], 1, 'up', 'dp')],
    across_cell_new=[LOOSE['B2_branch']['across_subspaces_of_cell'], rnd(LOOSE['B2_branch']['across_subspaces_of_cell'][0], 3, 'down'),
                     rnd(LOOSE['B2_branch']['across_subspaces_of_cell'][1], 3, 'up')],
    tau_new_0995_floor=[rnd(TAU['tau_new_0995'], 4, 'down', 'dp'), rnd(TAU['tau_new_0995'], 3, 'down', 'dp')],
    tau_new_L12pub=rnd(TAU['tau_new_L12pub'], 4, 'near', 'dp'),
    move_new_up1=rnd(TAU['move_new_pct'], 1, 'up', 'dp'),
    slopes_factor_new_eta018=rnd(SLOPES['new']['0.18']['factor'], 2),
    slopes_factor_old_eta018=rnd(SLOPES['old']['0.18']['factor'], 2),
    FRstar_L8_018_new=rnd(byc[(8, 0.18)]['FRstar_new'], 2),
    vac_range_new=[rnd(min(p['vac_new_smaller_branch'] for p in publ), 3), rnd(max(p['vac_new_smaller_branch'] for p in publ), 3)],
    vac_range_eta018_new=[rnd(min(p['vac_new_smaller_branch'] for p in publ if abs(p['eta'] - 0.18) < 1e-12), 3),
                          rnd(max(p['vac_new_smaller_branch'] for p in publ if abs(p['eta'] - 0.18) < 1e-12), 3)],
    pooled_min_ratio_new=sci(POOLSUM['ratio_new'][0], 3, 'down'),
)

# =============================================================================================== 17. needs a run
NEEDS = [
    dict(item='z_kblock first non-vacuous counts with the second-order branch (L=6: 48-50 of 50; L=8: 455-483 of 490)',
         why='only the first-order count and the bound at that count are stored; the count at which B2 first beats 1+w needs '
             'the bound at smaller block subspaces', policy='D9 (keep the stored first-order counts, labelled; add that the '
             'second-order branch can only lower them)'),
    dict(item='A(k,omega) certificate with the second-order branch for complex momentum probes k not in {0, pi} '
              '(0.77-1.75, 350-900x; 80 of the 112 cert_akw rows)',
         why='Lambda~_P is not stored for complex probes; Lambda~(k) = Lambda(-k) on the same S is not established '
             '(the sets S_k were not stored)', policy='D8 (first-order values kept, labelled, one sentence why)'),
    dict(item='the rung anatomy of the threshold for the second-order chain (44/38/11/4.6 %, 52/31/11/4.1 %, 1.10/3.29, 15 of 18 cells, 84 ladder rows)',
         why='the stored ladders (verdict_frontier.json, ladder_L*.json) are the rungs of the first-order proof; the '
             'second-order proof has no pointwise step', policy='D10 (keep, labelled as the first-order chain)'),
    dict(item='fitted-rule comparison with the new crossing: COMPUTED since 2026-11-16 (D20) by null_ws7_second_order.py '
              '(null_ws7_second_order.json); the entry is kept as the record of 2026-11-15',
         why='null_ws7.py fits c and p by a grid search (an optimisation), which this script does not run; an '
             'unverified replica reported 13/14 and 12/14', policy='none of D8-D11; D5 forbids a new run -- decided on '
                                                                   '2026-11-16 by D20 (the copy null_ws7_second_order.py)'),
]

# =============================================================================================== outputs
def dump(name, obj):
    with open(os.path.join(OUT, name), 'w', encoding='utf-8') as fh:
        json.dump(obj, fh, indent=1, default=lambda o: o if not isinstance(o, (np.bool_, np.floating, np.integer)) else o.item(),
                  sort_keys=False)


dump('grid_rows.json', dict(note='C3_grid.json, 204 rows; x = Lam_over_eta; B1/B2/T relative to ||phi||^2', checks=dev,
                            rows=rows))
dump('cells.json', dict(note='16 (L, eta) cells; err_* in the verify_jcp.py convention (all rows, FR_real, log-linear), '
                             'err_*_calib in the calib_oos.py one', cells=CELLS, crossing=CROSS, slopes=SLOPES, widths=WIDTHS,
                        gap_FRstar_minus_FR05=GAP, ties=TIES))
dump('published.json', dict(rows=publ, x_crosscheck_max_rel_dev=xcheck, monotone_smaller_branch=mono, monotone_B2=mono_B2,
                            L12_depth=L12DEPTH, L14_depth=L14DEPTH))
dump('five_percent.json', dict(per_cell={'L%d_eta%g' % (c['L'], c['eta']): dict(new=c['cert_new_below_0.05'], old=c['cert_old_below_0.05'],
                                                                                 FR05_true=c['FR05']) for c in CELLS},
                               fine_L8=fine, fig1=FIG1['cert_new_below_005']))
dump('supplement.json', dict(table_S1=TS1, table_S2=TS2, table_S2_summary=TS2_SUM, tau=TAU, looseness=LOOSE, oracle=ORACLE,
                             Sstar=star, akw=AKW, kblock=KB, kblock_factor=KFAC, fig1=FIG1))
dump('pooled.json', dict(summary=POOLSUM, deposited_378=DEP378,
                         rows=[{k: v for k, v in p.items()} for p in live]))
dump('checks.json', dict(never_below=NBSUM, D13=D13, guard_records=REC, misattributed_guard_record=MISSRC))
dump('summary.json', dict(printed_candidates=PRINT, needs_a_run=NEEDS, runtime_s=time.time() - T0,
                          inputs=dict(data=_rel(DATA), number_sources=_rel(NSRC))))

# =============================================================================================== 18. v2 (2026-11-16)
# Rules D17-D23 of the change (docs/KNOWN_DISCREPANCIES.md section 36), written after an independent check of the
# numbers. Same rules as above: stored fields only, arithmetic, no fit and no solver. The fitted-rule comparison of
# D20 is the separate copy null_ws7_second_order.py (run after this script; it reads FRstar_new from cells.json).
# Output: v2.json.
T0_V2 = time.time()
V2 = dict(note='rules D17-D23 of the change (2026-11-16; docs/KNOWN_DISCREPANCIES.md section 36); arithmetic on '
               'stored fields only')


def smaller_of(b1, b2):
    """D17: the smaller of the two leakage branches and its name (the first-order one on an exact tie)."""
    return (b2, 'second-order') if b2 < b1 else (b1, 'first-order')


def loose_stats(fn):
    """section 6 (LOOSE) for one branch function fn(row): CONVERGED grid rows with a true error."""
    conv = [r for r in rows if r['conv'] == 'CONVERGED' and r['relL1']]
    byP, byC = {}, {}
    for r in conv:
        byP.setdefault((r['L'], round(r['FR'], 6)), []).append(fn(r) / r['relL1'])
        byC.setdefault((r['L'], r['eta']), []).append(fn(r) / r['relL1'])
    q = {}
    for (L_, fr_), v in byP.items():
        if len(v) >= 2:
            q.setdefault(L_, []).append(max(v) / min(v))
    cc = {k: max(v) / min(v) for k, v in byC.items()}
    allv = [fn(r) / r['relL1'] for r in conv]
    return dict(fixedP_over_eta_max_L68=max(max(q[6]), max(q[8])),
                fixedP_over_eta_max_all=max(max(v) for v in q.values()),
                across_subspaces_of_cell=[min(cc.values()), max(cc.values())],
                across_subspaces_argmin='L%d_eta%g' % min(cc, key=cc.get),
                across_subspaces_argmax='L%d_eta%g' % max(cc, key=cc.get),
                all_converged_rows=[min(allv), max(allv)], n_rows=len(conv))


# ---- D17: every value quoted above the trivial bound -> both branches, the smaller marked (that is what is printed)
D17 = {}
tab3 = []
for p in publ:
    v, nm = smaller_of(p['B1'], p['B2'])
    tab3.append(dict(L=p['L'], eta=p['eta'], B1=p['B1'], B2=p['B2'], T=p['trivial'], B1_over_T=p['B1_over_T'],
                     B2_over_T=p['B2_over_T'], smaller=nm, printed_full=v / p['trivial'],
                     printed_2dp=rnd(v / p['trivial'], 2, 'near', 'dp'), present_2dp='%.2f' % p['vac_stored'],
                     src=p['src']))
D17['table_III'] = tab3
D17['table_III_summary'] = dict(
    n_cells=len(tab3),
    n_first_order_smaller=sum(t['smaller'] == 'first-order' for t in tab3),
    n_second_order_smaller=sum(t['smaller'] == 'second-order' for t in tab3),
    second_order_cells=[[t['L'], t['eta']] for t in tab3 if t['smaller'] == 'second-order'],
    range_all=[min(t['printed_full'] for t in tab3), max(t['printed_full'] for t in tab3)],
    range_eta018=[min(t['printed_full'] for t in tab3 if abs(t['eta'] - 0.18) < 1e-12),
                  max(t['printed_full'] for t in tab3 if abs(t['eta'] - 0.18) < 1e-12)],
    cells_whose_print_changes=[[t['L'], t['eta'], t['present_2dp'], t['printed_2dp']] for t in tab3
                               if t['present_2dp'] != t['printed_2dp']],
    all_above_one=all(t['printed_full'] > 1 for t in tab3))
p12 = pubtab[(12, 0.18)]
D17['L12_published_eta018'] = dict(B1=p12['B1'], B2=p12['B2'], T=p12['trivial'], smaller=smaller_of(p12['B1'], p12['B2'])[1],
                                   x=p12['x'], w=p12['w'], tau_new=p12['tau_new'], tau_old=p12['tau_old'],
                                   B1_3dp=rnd(p12['B1'], 3, 'near', 'dp'), B2_2dp=rnd(p12['B2'], 2, 'near', 'dp'),
                                   src=p12['src'])
D17['L14_published_depth250'] = [dict(eta=e, B1=pubtab[(14, e)]['B1'], B2=pubtab[(14, e)]['B2'], T=pubtab[(14, e)]['trivial'],
                                      smaller=smaller_of(pubtab[(14, e)]['B1'], pubtab[(14, e)]['B2'])[1],
                                      B1_2dp=rnd(pubtab[(14, e)]['B1'], 2, 'near', 'dp'),
                                      B2_2dp=rnd(pubtab[(14, e)]['B2'], 2, 'near', 'dp'), src=pubtab[(14, e)]['src'])
                                 for e in ETAS]
orc = []
for fr_ in sorted(set(round(z['FR'], 4) for z in zr)):
    zb = [z for z in zr if z['ranking'] == 'repo' and abs(z['FR'] - fr_) < 1e-9]
    zo = [z for z in zr if z['ranking'] == 'krylov' and abs(z['FR'] - fr_) < 1e-9]
    if not zb or not zo:
        continue
    zb, zo = zb[0], zo[0]
    vb, nbn = smaller_of(zb['B1'], zb['B2'])
    vo, non = smaller_of(zo['B1'], zo['B2'])
    orc.append(dict(FR=fr_, born=dict(B1=zb['B1'], B2=zb['B2'], smaller=nbn, value=vb, T=zb['trivial'],
                                      beats_T=vb < zb['trivial'], src=zb['src']),
                    oracle=dict(B1=zo['B1'], B2=zo['B2'], smaller=non, value=vo, T=zo['trivial'],
                                beats_T=vo < zo['trivial'], src=zo['src']),
                    oracle_worse_smaller_branch=vo > vb))
D17['oracle'] = dict(rows=orc,
                     oracle_worse_at_every_FR_from_050=all(r['oracle_worse_smaller_branch'] for r in orc if r['FR'] >= 0.5 - 1e-9),
                     FR_checked_from_050=[r['FR'] for r in orc if r['FR'] >= 0.5 - 1e-9],
                     oracle_better_at_030=[not r['oracle_worse_smaller_branch'] for r in orc if abs(r['FR'] - 0.3) < 1e-9][0])
tie_rows = []
for blk in J('c3_frontier', 'ties.json'):
    for p in blk['points']:
        for k in sorted(p):
            if k.startswith('f32_eta') and isinstance(p[k], dict):
                k64 = k.replace('f32', 'f64')
                w32, w64 = p['f32_w_S'], p['f64_w_S']
                b32, b64 = p[k]['B'], p[k64]['B']
                c32, c64 = B2(w32, x_from_B1(b32, w32)), B2(w64, x_from_B1(b64, w64))
                s32, n32 = smaller_of(b32, c32)
                s64, n64 = smaller_of(b64, c64)
                tie_rows.append(dict(L=blk['L'], FR=p.get('FR'), eta=float(k.replace('f32_eta', '')),
                                     B1_rel=abs(b32 - b64) / b64, B2_rel=abs(c32 - c64) / c64,
                                     smaller_f32=n32, smaller_f64=n64, smaller_value_f64=s64, T=1.0 + w64,
                                     smaller_beats_T=s64 < 1.0 + w64, smaller_rel=abs(s32 - s64) / s64,
                                     cert_rel=abs(min(1 + w32, c32) - min(1 + w64, c64)) / min(1 + w64, c64)))
tmax = max(tie_rows, key=lambda t: t['smaller_rel'])
D17['ties'] = dict(src=TIES['src'], max_smaller_rel=tmax['smaller_rel'],
                   at='L=%d FR=%g eta=%g (smaller: %s, %r against T=%r)' % (tmax['L'], tmax['FR'], tmax['eta'],
                                                                            tmax['smaller_f64'], tmax['smaller_value_f64'], tmax['T']),
                   printed_up_2sf=sci(tmax['smaller_rel'], 2, 'up'),
                   max_B1_rel=TIES['B1_max_rel'], printed_B1_up_2sf=sci(TIES['B1_max_rel'], 2, 'up'),
                   max_B2_rel=TIES['B2_max_rel'], max_cert_new_rel=TIES['cert_new_max_rel'],
                   n_points=len(tie_rows), rows=tie_rows)
x250, x300 = r250['LamHat_over_eta'], r300['LamHat_over_eta']
s250 = smaller_of(B1(w12, x250), B2(w12, x250))
s300 = smaller_of(B1(w12, x300), B2(w12, x300))
D17['depth_L12_published_eta018'] = dict(src=L12DEPTH['src'], smaller_n250=s250[1], smaller_n300=s300[1],
                                         smaller_rel_change=rel(s300[0], s250[0]), printed_2sf=sci(rel(s300[0], s250[0]), 2),
                                         B2_rel_change=L12DEPTH['B2_rel_change'], Lout_rel_change=L12DEPTH['Lout_rel_change'])
D17['depth_L14_eta010'] = dict(src=L14DEPTH['src'], smaller_both_depths=L14DEPTH['smaller_branch_both_depths'],
                               depths=L14DEPTH['depths'], B1_fall_pct=L14DEPTH['B1_fall_pct'], B2_fall_pct=L14DEPTH['B2_fall_pct'])
f1 = []
for q in fig:
    v, nm = smaller_of(q['B1'], q['B2'])
    f1.append(dict(FR=q['FR'], B1=q['B1'], B2=q['B2'], T=q['trivial'], smaller=nm, leak_D17=v, relL1=q['relL1'],
                   missed=q['missed'], src=q['src']))
fs = FIG1['FRstar_new']
s8v, s8n = smaller_of(st8['B1'], st8['B2'])
D17['fig1'] = dict(
    rule='D17: the plotted leakage curve is the smaller of the two branches at each row (D14: the certificate is min of it and 1+w)',
    rows=f1, first=f1[0]['leak_D17'], last=f1[-1]['leak_D17'], first_branch=f1[0]['smaller'], last_branch=f1[-1]['smaller'],
    first_printed_2sf=sci(f1[0]['leak_D17'], 2), last_printed_2sf=sci(f1[-1]['leak_D17'], 2),
    first_FR_where_second_order_smaller=min(q['FR'] for q in f1 if q['smaller'] == 'second-order'),
    last_FR_where_first_order_smaller=max(q['FR'] for q in f1 if q['smaller'] == 'first-order'),
    max_leak_D17=max(q['leak_D17'] for q in f1), ymax_present=80, ymax_80_still_holds=max(q['leak_D17'] for q in f1) < 80,
    max_B2_alone=FIG1['B2_max'], ymax_if_B2_alone_must_exceed=FIG1['B2_max'],
    min_plotted=min(min(q['relL1'] for q in f1), min(q['missed'] for q in f1), min(q['leak_D17'] for q in f1)),
    ymin_present=5e-9,
    FRstar_new=fs, FRstar_new_6dp=rnd(fs, 6, 'near', 'dp'), FRstar_new_2dp=rnd(fs, 2, 'near', 'dp'),
    shaded_right_edge=fs,
    label_x_same_offsets_as_now=[rnd(fs - (0.794049 - 0.785), 3, 'near', 'dp'), rnd(fs + (0.805 - 0.794049), 3, 'near', 'dp')],
    max_leak_D17_on_rows_FR_055_to_FRstar=max(q['leak_D17'] for q in f1 if 0.55 <= q['FR'] <= fs),
    star=dict(B1=st8['B1'], B2=st8['B2'], smaller=s8n, leak_D17=s8v, stored_bound_theoremB=st8['B1'], src=st8['src']))
D17['looseness'] = dict(B1=loose_stats(lambda r: r['B1']), B2=loose_stats(lambda r: r['B2']),
                        smaller=loose_stats(lambda r: min(r['B1'], r['B2'])), cert_new=loose_stats(lambda r: r['cert_new']))
assert rel(D17['looseness']['B2']['fixedP_over_eta_max_all'], LOOSE['B2_branch']['fixedP_over_eta_max_all']) < 1e-12
assert rel(D17['looseness']['B1']['across_subspaces_of_cell'][0], LOOSE['B1_branch']['across_subspaces_of_cell'][0]) < 1e-12
V2['D17'] = D17

# ---- D18: momentum-resolved probes, L=8, FR=0.85 (cert_akw.json): real channels k=0, pi use the certificate of D6,
#      complex channels the first-order branch (D8)
akr, akc = [], []
for a in ak:
    r = CA[a['i']]
    w, x = r['w_S'], r['Lambda_hat_over_eta']
    if a['real_probe']:
        b2 = B2(w, x)
        c = min(1.0 + w, b2)
        akr.append(dict(i=a['i'], branch=a['branch'], k_over_pi=a['k_over_pi'], w=w, x=x, cert=c, B2=b2,
                        nontrivial=b2 < 1.0 + w, true=a['true'], cert_over_true=c / a['true'],
                        first_order=a['rel_upper'], first_order_over_true=a['slack']))
    else:
        akc.append(dict(i=a['i'], branch=a['branch'], k_over_pi=a['k_over_pi'], w=w, x=x, rel_upper=a['rel_upper'],
                        slack=a['slack'], ratio_recomputed=a['rel_upper'] / a['true'], nontrivial=a['nontrivial'],
                        true=a['true']))
assert max(rel(q['ratio_recomputed'], q['slack']) for q in akc) < 1e-9
akw_live_real = [p for p in live if p.get('akw') and p['real']]
V2['D18'] = dict(
    src=AKW['src'] + '; w_S, Lambda_hat_over_eta for the second-order branch',
    real=dict(n_channels=len(akr), k_over_pi=sorted(set(q['k_over_pi'] for q in akr)),
              n_symmetry_distinct=len(set((round(q['w'], 8), rnd(q['x'], 4)) for q in akr)),
              cert=[min(q['cert'] for q in akr), max(q['cert'] for q in akr)],
              cert_over_true=[min(q['cert_over_true'] for q in akr), max(q['cert_over_true'] for q in akr)],
              all_nontrivial=all(q['nontrivial'] for q in akr),
              first_order_same_channels=[min(q['first_order'] for q in akr), max(q['first_order'] for q in akr)],
              first_order_over_true_same_channels=[min(q['first_order_over_true'] for q in akr),
                                                   max(q['first_order_over_true'] for q in akr)],
              rows=akr),
    complex=dict(n_channels=len(akc), k_over_pi=sorted(set(q['k_over_pi'] for q in akc)),
                 n_symmetry_distinct=len(set((round(q['w'], 8), rnd(q['x'], 4)) for q in akc)),
                 first_order=[min(q['rel_upper'] for q in akc), max(q['rel_upper'] for q in akc)],
                 first_order_over_true=[min(q['slack'] for q in akc), max(q['slack'] for q in akc)],
                 all_nontrivial=all(q['nontrivial'] for q in akc), rows=akc),
    all16_first_order=dict(rel_upper=AKW['rel_upper'], slack=AKW['slack']),
    pooled=dict(total=POOLSUM['nontrivial_new_total'], real_rows=POOLSUM['n_real_live'], real_rows_B2_wins=POOLSUM['real_rows_B2_wins'],
                complex_rows=POOLSUM['n_complex_live'], complex_rows_B1_wins=POOLSUM['complex_rows_B1_wins'],
                akw_real_rows=len(akw_live_real), akw_real_rows_B2_wins=sum(p['b2'] < p['triv'] for p in akw_live_real),
                complex_count_uses='stored Lambda_P with the first-order branch (conservative: Lambda^ = min(Lambda, Lambda~) <= Lambda)',
                old_count=POOLSUM['leak_wins_old']))

# ---- D19: where the certificate first falls below 0.05, every stored (L, eta), both branches
d19 = {}
for c in CELLS:
    row = {}
    for lab in ('old', 'new'):
        fb = c['cert_%s_below_0.05' % lab]
        if fb['reached']:
            row[lab] = dict(reached=True, bracket=[fb['lo'], fb['hi']], values=[fb['lo_value'], fb['hi_value']],
                            conv=[fb['conv_lo'], fb['conv_hi']], FR_loglin=fb['FR_interp_loglin'],
                            FR_real=[fb['lo_FR_real'], fb['hi_FR_real']], back_above_later=fb['later_above'],
                            printed=fb['printed'], rows=list(fb['rows']))
        else:
            row[lab] = dict(reached=False, last_FR=fb['above_FR'], value_at_last=fb['value_at_last'], conv_last=fb['conv_last'],
                            printed=fb['printed'])
    d19['L%d_eta%g' % (c['L'], c['eta'])] = row
V2['D19'] = dict(cells=d19,
                 L6_new_not_reached_at_eta=[e for e in ETAS if not byc[(6, e)]['cert_new_below_0.05']['reached']],
                 L6_new_reached_at_eta=[e for e in ETAS if byc[(6, e)]['cert_new_below_0.05']['reached']],
                 L6_old_reached_at_eta=[e for e in ETAS if byc[(6, e)]['cert_old_below_0.05']['reached']],
                 eta018_new={str(L_): d19['L%d_eta0.18' % L_]['new'].get('printed') for L_ in LS},
                 eta015_new={str(L_): d19['L%d_eta0.15' % L_]['new'].get('printed') for L_ in LS})

# ---- D21: "never below the true error", every stored row, with the qualifier and the excluded rows counted
CEN = []


def cen(source, key, rel_err, err, b2, b1, real=True, full=False, lam_pos=True):
    if full:
        cls = 'excluded: full sector (A_S = A by algebra; the stored error is reference drift)'
    elif rel_err is None:
        cls = 'no stored error'
    elif rel_err <= 1e-10:
        cls = 'excluded: true relative error <= 1e-10'
    elif not real:
        cls = 'complex probe: first-order branch only (D8)'
    else:
        cls = 'tested'
    CEN.append(dict(source=source, key=key, rel=rel_err, err=err, B2=b2, B1=b1, cls=cls, lam_pos=lam_pos,
                    ratio_B2=(b2 / err) if (err and b2 is not None) else None,
                    ratio_B1=(b1 / err) if (err and b1 is not None) else None))


for r in rows:
    cen('C3_grid.json', r['src'], r['relL1'], r['relL1'], r['B2'], r['B1'], lam_pos=r['x'] > 0)
for f in sorted(glob.glob(os.path.join(DATA, 'c3_frontier', 'published_fraction', 'pt_L*_PUB.json'))):
    d = json.load(open(f, encoding='utf-8'))
    d = d[0] if isinstance(d, list) else d
    for i, q in enumerate(d['rows']):
        if 'relL1' in q:
            cen(os.path.basename(f), 'rows[%d] n=%d eta=%g' % (i, q['n'], q['eta']), q['relL1'], q['relL1'],
                B2(d['w_S'], q['LamHat_over_eta']), B1(d['w_S'], q['LamHat_over_eta']), lam_pos=q['LamHat_over_eta'] > 0)
for i, q in enumerate(L14['rows']):
    if 'relL1' in q:
        cen('L14_FR008.json', 'rows[%d] n=%d eta=%g' % (i, q['n'], q['eta']), q['relL1'], q['relL1'],
            B2(L14['w_S'], q['LamHat_over_eta']), B1(L14['w_S'], q['LamHat_over_eta']), lam_pos=q['LamHat_over_eta'] > 0)
for q in fine:
    cen('z_fine_L8.txt', q['src'], q['relL1'], q['relL1'], q['B2'], q['B1'])
ZF6 = J('c3_frontier', 'adversarial', 'z_fine_L6.json')
for i, q in enumerate(ZF6):
    b1z = (1.0 + math.sqrt(min(q['w_S'], 1.0))) * (q['phiQ_rel'] + q['Lamhat_over_eta'])
    assert rel(b1z, q['relR5']) < 1e-6 or q['relR5'] < 1e-12
    cen('z_fine_L6.json', '[%d] k=%d FR=%g eta=%g' % (i, q['k'], q['FR'], q['eta']), q['relR0'], q['relR0'],
        (q['phiQ_rel'] + q['Lamhat_over_eta']) ** 2, q['relR5'], lam_pos=q['Lamhat_over_eta'] > 0)
for z in zr:
    cen('z_rank_L6.json', z['src'], z['true'], z['true'], z['B2'], z['B1'], lam_pos=z['x'] > 0)
for s in star:
    cen('support_witness.json', s['src'], s['relL1'], s['relL1'], s['B2'], s['B1'], lam_pos=s['x'] > 0)
for p in POOL:
    cen('pool:' + p['kind'] if p['kind'] in ('cert_stress.json', 'cert_akw.json', 'cert_teqsci.json') else
        'pool:' + p['src'].split(' ')[0], p['src'], p['rel'], p['err'], p['b2'] if p['real'] else None, p['b1'],
        real=p['real'], lam_pos=p['lam'] > 0)
for fn in ('cert_frontier_6-8.json', 'cert_frontier_10.json'):
    d = J(fn)
    for Lk, cv in d['curves'].items():
        for i, q in enumerate(cv['points']):
            full = bool(q.get('S_is_full_sector'))
            b2q = None if full else (q['norm_phiQ'] + q['Lambda'] / q['eta']) ** 2
            b1q = None if full else q['upper_unclipped']
            cen(fn, 'curves/%s/points[%d] frac=%g' % (Lk, i, q['frac']), q.get('relL1_true'), q.get('L1_true'),
                b2q, b1q, full=full, lam_pos=(not full) and q['Lambda'] > 0)
for Lf in (4, 6, 8):
    d = J('c3_frontier', 'ladder_L%d.json' % Lf)
    for ci, c in enumerate(d['cases']):
        n2 = c.get('norm_phi2', d.get('norm_phi2'))
        for ri, q in enumerate(c['rows']):
            cen('ladder_L%d.json' % Lf, 'cases[%d] (%s) rows[%d] eta=%g' % (ci, c['tag'], ri, q['eta']), q['R0'] / n2,
                q['R0'], (c['norm_phiQ'] + q['Lam_over_eta']) ** 2, q['R5'], lam_pos=q['Lam'] > 0)
tested = [q for q in CEN if q['cls'] == 'tested']
cplx = [q for q in CEN if q['cls'].startswith('complex')]
excl_r = [q for q in CEN if q['cls'] == 'excluded: true relative error <= 1e-10']
excl_f = [q for q in CEN if q['cls'].startswith('excluded: full sector')]
noerr = [q for q in CEN if q['cls'] == 'no stored error']
lpos = [q for q in tested if q['lam_pos']]
lzero = [q for q in tested if not q['lam_pos']]
qmin = min(tested, key=lambda q: q['ratio_B2'])
qminp = min(lpos, key=lambda q: q['ratio_B2'])
per = {}
for q in CEN:
    s = per.setdefault(q['source'], dict(n=0, tested=0, excluded_roundoff=0, excluded_full_sector=0, no_error=0,
                                         complex_first_order=0, min_ratio_B2=None))
    s['n'] += 1
    if q['cls'] == 'tested':
        s['tested'] += 1
        s['min_ratio_B2'] = q['ratio_B2'] if s['min_ratio_B2'] is None else min(s['min_ratio_B2'], q['ratio_B2'])
    elif q['cls'].startswith('excluded: true'):
        s['excluded_roundoff'] += 1
    elif q['cls'].startswith('excluded: full'):
        s['excluded_full_sector'] += 1
    elif q['cls'] == 'no stored error':
        s['no_error'] += 1
    else:
        s['complex_first_order'] += 1
V2['D21'] = dict(
    statement='B2 (the second-order branch) is never below the stored true error on every stored real-probe row whose '
              'true relative L1 error exceeds 1e-10; rows with zero error (round-off, or the full sector, where both '
              'sides vanish up to reference drift) are excluded and counted',
    n_rows_all=len(CEN), n_tested=len(tested),
    violations_strict=sum(q['err'] > q['B2'] for q in tested),
    violations_tol_1e9=sum(q['err'] > q['B2'] * (1 + 1e-9) for q in tested),
    min_ratio=qmin['ratio_B2'], min_ratio_at=qmin['source'] + ' ' + qmin['key'],
    n_tested_zero_leakage=len(lzero),
    zero_leakage_ratio_range=[min(q['ratio_B2'] for q in lzero), max(q['ratio_B2'] for q in lzero)] if lzero else None,
    min_ratio_nonzero_leakage=qminp['ratio_B2'], min_ratio_nonzero_leakage_at=qminp['source'] + ' ' + qminp['key'],
    n_excluded_roundoff=len(excl_r), n_excluded_full_sector=len(excl_f), n_no_stored_error=len(noerr),
    excluded_roundoff_B2_below_stored_error=sum(1 for q in excl_r if q['B2'] is not None and q['err'] > q['B2']),
    excluded_roundoff_B1_below_stored_error=sum(1 for q in excl_r if q['B1'] is not None and q['err'] > q['B1']),
    excluded_full_sector=[q['source'] + ' ' + q['key'] + ' stored relL1=%r' % q['rel'] for q in excl_f],
    complex_rows=len(cplx), complex_violations_B1=sum(q['err'] > q['B1'] for q in cplx),
    complex_min_ratio_B1=min(q['ratio_B1'] for q in cplx) if cplx else None,
    per_source=per)
dep_r = [q for q in CEN if q['source'] in ('pool:cert_stress.json', 'pool:cert_akw.json', 'pool:cert_teqsci.json')]
V2['D21']['deposited_378'] = dict(
    n=len(dep_r), tested=sum(q['cls'] == 'tested' for q in dep_r),
    complex_first_order=sum(q['cls'].startswith('complex') for q in dep_r),
    excluded_roundoff=sum(q['cls'].startswith('excluded: true') for q in dep_r),
    each_leakage_branch_holds_on_ratio_defined=sum(1 for q in dep_r if q['rel'] > 1e-10 and q['err'] <= min(
        x_ for x_ in (q['B2'], q['B1']) if x_ is not None)),
    n_ratio_defined=sum(1 for q in dep_r if q['rel'] > 1e-10),
    roundoff_rows=[dict(key=q['key'], rel=q['rel'], err=q['err'], B2=q['B2'], B1=q['B1']) for q in dep_r
                   if q['cls'].startswith('excluded: true')])

# ---- D22: the facts the new wording must state, recomputed
c80 = byc[(8, 0.18)]
V2['D22'] = dict(
    looseness_second_order_branch=D17['looseness']['B2'], looseness_smaller_branch=D17['looseness']['smaller'],
    looseness_first_order_branch=D17['looseness']['B1'],
    B2_over_true_converged_nontrivial=LOOSE['B2_branch']['converged_nontrivial_rows'],
    cert_over_true_least_loose_per_cell_conv=[CROSS['new']['least_loose_conv']['min'], CROSS['new']['least_loose_conv']['max']],
    ratio_at_crossing=[CROSS['new']['ratio_at_crossing']['min'], CROSS['new']['ratio_at_crossing']['max']],
    prop2=dict(tau_new_1=TAU['tau_new_1'], tau_old_1=TAU['tau_old_1'],
               threshold='W/(2 sqrt(2) eta) for real problems (tau_new(1) = sqrt 2); W/(2 eta) stays for the first-order branch'),
    fig1_fine_points_B2=FIG1['fine_points_B2'], fig1_fine_points_still_locate_005=FIG1['fine_points_still_locate_005'],
    fig1_cert_005=dict(bracket_FR=[c80['cert_new_below_0.05']['lo'], c80['cert_new_below_0.05']['hi']],
                       bracket_FR_real=[c80['cert_new_below_0.05']['lo_FR_real'], c80['cert_new_below_0.05']['hi_FR_real']],
                       values=[c80['cert_new_below_0.05']['lo_value'], c80['cert_new_below_0.05']['hi_value']]),
    five_percent_eta018=V2['D19']['eta018_new'], five_percent_eta015=V2['D19']['eta015_new'],
    shift_FR=c80['cert_new_below_0.05']['FR_interp_loglin'] - c80['FR05'],
    shift_FR_bracket=[c80['cert_new_below_0.05']['lo'] - c80['FR05'], c80['cert_new_below_0.05']['hi'] - c80['FR05']])

# ---- D23: the guard records of the L=14 branch values (sm_frontier.tex, "17.21, 11.18, 9.21 and 6.50")
lmap = {}
for i, r in enumerate(L14['rows']):
    if r['n'] == deep14:
        lmap['%.2f' % r['bound_leak']] = (i, r)
d23 = []
for rr in NS:
    if rr['tex'] == 'sm_frontier.tex' and rr['lit'] in ('17.21', '11.18', '9.21', '6.50'):
        i, r = lmap[rr['lit']]
        right = '[0]/rows[%d]/bound_leak' % i
        mm = re.search(r'rows\[(\d+)\]/(\w+)', rr['key'])
        pr = L14['rows'][int(mm.group(1))]
        d23.append(dict(lit=rr['lit'], line=rr['line'], record_key=rr['key'], record_value=rr['value'],
                        record_points_to='n=%d eta=%g field %s' % (pr['n'], pr['eta'], mm.group(2)),
                        printed_quantity='first-order branch (bound_leak) at the deepest depth n=%d, eta=%g' % (r['n'], r['eta']),
                        right_key=right, right_value=r['bound_leak'], record_is_right=rr['key'] == right))
V2['D23'] = dict(file='c3_frontier/published_fraction/L14_FR008.json', records=d23)

# ---- guard records whose printed literal changes or goes (for the guards); present literals of the number sheet
CHANGED_LITS = {
    'main.tex': ['1.2', '7.3'],
    'sec_1.tex': ['1.66', '0.98', '1.2', '7.3', '5.9', '0.77', '1.75', '350', '900'],
    'sec_3_body.tex': ['1.8', '2.0', '23', '148', '1.6', '2.4'],
    'sec_3_app.tex': ['61', '376', '378'],
    'sec_5_body.tex': ['1.2', '7.3', '1.66', '0.9625', '61', '0.718', '2.8', '0.219', '0.77', '1.75', '350', '900', '0.98',
                       '0.977', '0.983', '2.7', '5.9', '46', '0.9280', '7.2'],
    'sec_5_app.tex': ['46', '7.8'],
    'sec_9.tex': ['3'],
    'sec_10.tex': ['1.2', '7.3', '0.77', '1.75', '350', '900', '1.66', '0.98'],
    'sm_frontier.tex': ['0.98', '0.45', '0.9625', '0.92', '0.82', '0.90', '0.75', '0.80', '0.60', '0.70', '0.40', '0.45',
                        '0.95', '0.97', '0.85', '0.503', '0.220', '0.137', '0.718', '0.878', '0.219', '0.084', '0.060',
                        '0.039', '3802', '3830', '3855', '0.977', '0.983', '2.8', '9.2', '1.2', '7.3', '2.7', '5.9', '46',
                        '5.65', '3.50', '2.52', '1.60', '3.7', '1.8', '2.0', '23', '148', '1.6', '2.4', '0.067', '0.296',
                        '0.9280', '0.928', '7.2', '4.1', '45', '9.4', '3.8', '0.911', '61', '0.77', '1.75', '350', '900'],
    'sm_leakage.tex': ['7.8'],
    'figs/caption_thm1iii_violation.tex': ['1.18', '61'],
    'figs/fig_results_native.tex': ['20', '0.039', '0.79', '0.977', '0.983'],
}
gr = []
for rr in NS:
    if rr['lit'] in CHANGED_LITS.get(rr['tex'], []):
        gr.append(dict(tex=rr['tex'], line=rr['line'], lit=rr['lit'], exp=rr['exp'], file=rr['file'], key=rr.get('key'),
                       expr=rr.get('expr'), round=rr.get('round'), value=rr.get('value')))
V2['guard_records_touched'] = gr

# ---- smaller facts the sheet needs
f3830 = [q for q in fine if q['S'] == 3830][0]
f3802 = [q for q in fine if q['S'] == 3802][0]
f3855 = [q for q in fine if q['S'] == 3855][0]
V2['presheet'] = dict(
    ceiling_0060=dict(stored_first_order=f3830['B1'], printed='0.060', ceiling_3dp=rnd(f3830['B1'], 3, 'up', 'dp'),
                      nearest_3dp=rnd(f3830['B1'], 3, 'near', 'dp'), src=f3830['src']),
    fine_first_order=[f3802['B1'], f3830['B1'], f3855['B1']],
    figS10_label=dict(label_y=1.1e-5, ymin=1e-6, min_ratio_new=POOLSUM['ratio_new'][0], min_ratio_new_at=POOLSUM['ratio_new_min_at'],
                      rows_below_1e4=sum(1 for p in live if p['ratio_new'] < 1e-4),
                      min_ratio_old=POOLSUM['ratio_old'][0]),
    L8_eta018_cert_005=dict(lo=c80['cert_new_below_0.05']['lo'], hi=c80['cert_new_below_0.05']['hi'],
                            lo_S=[q['S'] for q in rows if q['idx'] == c80['cert_new_below_0.05']['rows'][0]][0],
                            hi_S=[q['S'] for q in rows if q['idx'] == c80['cert_new_below_0.05']['rows'][1]][0],
                            lo_value_up_2sf=rnd(c80['cert_new_below_0.05']['lo_value'], 2, 'up'),
                            hi_value_up_2sf=rnd(c80['cert_new_below_0.05']['hi_value'], 2, 'up'),
                            FR05=c80['FR05']))
V2['runtime_s'] = time.time() - T0_V2
dump('v2.json', V2)
print('v2: D17 table III %d/%d first/second-order smaller; ties %s; D18 real %d (%.4g..%.4g, x%.4g..%.4g), complex %d (%.4g..%.4g, x%.4g..%.4g)'
      % (D17['table_III_summary']['n_first_order_smaller'], D17['table_III_summary']['n_second_order_smaller'],
         D17['ties']['printed_up_2sf'], len(akr), V2['D18']['real']['cert'][0], V2['D18']['real']['cert'][1],
         V2['D18']['real']['cert_over_true'][0], V2['D18']['real']['cert_over_true'][1], len(akc),
         V2['D18']['complex']['first_order'][0], V2['D18']['complex']['first_order'][1],
         V2['D18']['complex']['first_order_over_true'][0], V2['D18']['complex']['first_order_over_true'][1]))
print('v2: D21 %d rows, %d tested, %d violations, min %.10g; excluded %d round-off + %d full sector; %d complex; %d no error'
      % (len(CEN), len(tested), V2['D21']['violations_strict'], V2['D21']['min_ratio'], len(excl_r), len(excl_f),
         len(cplx), len(noerr)))

# ----------------------------------------------------------------------------------------------- console
print('grid: %d rows; B1 reproduces bound_leak to %.1e; lemma violations %d' % (len(rows), dev['B1_vs_bound_leak'], len(lemma_viol)))
print('non-trivial rows: old %d, new %d' % (sum(r['nontrivial_old'] for r in rows), sum(r['nontrivial_new'] for r in rows)))
for lab in ('old', 'new'):
    s = CROSS[lab]
    print('%s crossing error %.4e .. %.4e median %.4e spread %.4f geo %.3e half %.3f ; (1+w)/err %.1f..%.1f ; drift >= %.2f (L=%d eta=%g)'
          % (lab, s['min'], s['max'], s['median'], s['spread'], s['geo_mid'], s['half_spread'], s['ratio_at_crossing']['min'],
             s['ratio_at_crossing']['max'], s['drift_margin_PUB'][0], s['drift_margin_PUB'][1], s['drift_margin_PUB'][2]))
print('printed:', json.dumps(PRINT, default=str)[:3000])
print('pool:', json.dumps(POOLSUM, default=str)[:1500])
print('D13:', json.dumps(D13, default=str))
print('never below: %d rows, %d violations, smallest %s' % (NBSUM['n_rows'], NBSUM['violations'], json.dumps(NBSUM['smallest'], default=str)))
print('runtime %.2f s' % (time.time() - T0))
