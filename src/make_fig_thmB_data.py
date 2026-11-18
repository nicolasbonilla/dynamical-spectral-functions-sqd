# -*- coding: utf-8 -*-
"""make_fig_thmB_data.py -- the blue series of Fig. S10 (fig:thm1iii-violation) of the Supplement of the manuscript
prepared for JCP: for each pooled subspace, the measured error ||A - A_S||_1 divided by the certificate of Theorem 1,
as the two pgfplots tables paper/figs/n3_thmB_frac.dat (x = subspace fraction) and paper/figs/n3_thmB_eta.dat
(x = eta).

2026-11-16, second-order branch of Theorem 1 (rules D6, D8, D18 of the change, docs/KNOWN_DISCREPANCIES.md
section 36).  Until this date both tables were written by src/recovered/build_n3.py with the first-order certificate
min(||phi||^2 (1+w_S), (||phi||+||phi_S||)(||phi_Q|| + Lambda_S/eta)).  This script writes them with the
certificate the theorem now gives:
  * real probes (every row except the complex ones below):
        min(||phi||^2 (1+w_S), (||phi_Q|| + Lambda_S/eta)^2)                   (D6; Lambda~_S = Lambda_S)
  * complex momentum-resolved probes (cert_akw.json rows with k_over_pi not 0 or 1):
        min(||phi||^2 (1+w_S), (||phi||+||phi_S||)(||phi_Q|| + Lambda_S/eta))   (D8: Lambda~_S was not stored)
It reads the same five stored files as build_n3.py, in the same order (thm3_results.json, cert_stress.json,
cert_akw.json, cert_teqsci.json, then thm3_L8b.json), keeps the same live rows (relative L1 > 1e-10) and writes
the same format ('x y' header, '%.6g %.6g', written in text mode, so with the platform's line ends, as
build_n3.py writes its tables).  build_n3.py still writes the other eight n3_*.dat.
NO SIMULATION: arithmetic on stored fields only (one thread, well under a second).

Checks, all asserted before anything is written:
  * --old rebuilds both tables with build_n3.py's own certificate and compares them BYTE FOR BYTE, line ends
    normalised to LF on both sides, with the tables of deposit 3.4.1, kept at _superseded/arxiv_v4/paper/figs/n3_thmB_{frac,eta}.dat (so the rows, their
    order and the format are build_n3.py's); it writes nothing.
  * every new ratio equals data/c3_frontier/second_order/pooled.json rows[].ratio_new (matched by source row) to
    1e-12 relative; 606 live rows; no ratio above 1 (+1e-9); the count of rows where the leakage entry beats the
    trivial cap equals the re-scoring's 85 = 57 real + 28 complex.

Usage (from the repository root):  python src/make_fig_thmB_data.py [--old] [--dry]

Deposit 3.5.0: the maker of the JCP package (2026-11-16) with its paths taken from its place in the repository;
nothing else changed.  src/check_figures.py runs it after src/recovered/build_n3.py and compares both tables byte
for byte.
"""
import io
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))                      # <repo>/src
REL = os.path.dirname(HERE)                                            # <repo>
DATA = os.path.join(REL, 'data')
OUT = os.path.join(REL, 'paper', 'figs')
OLD = os.path.join(REL, '_superseded', 'arxiv_v4', 'paper', 'figs')    # the first-order tables of deposit 3.4.1
POOLED = os.path.join(DATA, 'c3_frontier', 'second_order', 'pooled.json')
NOISE = 1e-10                                                          # build_n3.py: no ratio at or below this


def load(*p):
    with io.open(os.path.join(*p), encoding='utf-8') as fh:
        return json.load(fh)


def rows():
    """the pooled rows, in build_n3.py's order: (src key of second_order/pooled.json, frac, eta, err, relerr,
    old certificate as build_n3.py takes it, new certificate, real?, leakage entry beats the trivial cap?)"""
    out = []
    for i, r in enumerate(load(DATA, 'thm1iii_violation', 'thm3_results.json')):
        w = r['wS']
        nq = math.sqrt(r['norm_phi2']) * math.sqrt(max(1.0 - w, 0.0))
        b2 = (nq + r['Lam'] / r['eta']) ** 2
        out.append(dict(src='thm1iii_violation/thm3_results.json [%d]' % i, x=r['frac'], eta=r['eta'],
                        err=r['LHS'], rel=r['relL1'], old=min(r['RHS_new'], r['RHS_triv']),
                        new=min(r['RHS_triv'], b2), real=True, wins=b2 < r['RHS_triv']))
    for fn in ('cert_stress.json', 'cert_akw.json', 'cert_teqsci.json'):
        for i, r in enumerate(load(DATA, fn)['rows']):
            nq, triv, b1 = r['norm_phiQ'], r['upper_trivial'], r['upper_unclipped']
            assert abs(min(triv, b1) - r['upper']) <= 1e-12 * r['upper'], (fn, i)
            b2 = (nq + r['Lambda'] / r['eta']) ** 2
            kpi = r.get('k_over_pi')
            real = not (fn == 'cert_akw.json' and not (abs(kpi) < 1e-9 or abs(kpi - 1.0) < 1e-9))
            leak = b2 if real else b1
            out.append(dict(src='%s /rows[%d]' % (fn, i), x=r['frac'], eta=r['eta'], err=r['L1_true'],
                            rel=r['relL1_true'], old=r['upper'], new=min(triv, leak), real=real, wins=leak < triv))
    src = io.open(os.path.join(REL, 'src', 'recovered', 'thm3_L8b.py'), encoding='utf-8').read()
    m = re.search(r"^L,U,eta\s*=\s*([0-9.]+)\s*,\s*([0-9.]+)\s*,\s*([0-9.]+)\s*$", src, re.M)
    assert m, 'could not read L,U,eta from thm3_L8b.py'
    eta8 = float(m.group(3))
    for i, r in enumerate(load(DATA, 'thm1iii_violation', 'thm3_L8b.json')):
        w = r['wS']
        n2 = r['LB'] / (1.0 - w)                         # LB = ||phi||^2 (1 - w_S)
        triv = n2 * (1.0 + w)
        b2 = (math.sqrt(n2) * math.sqrt(1.0 - w) + r['Lam_over_eta']) ** 2   # Lam_over_eta: absolute Lambda/eta
        assert abs(r['Lam'] / r['Lam_over_eta'] - eta8) < 1e-12
        out.append(dict(src='thm1iii_violation/thm3_L8b.json [%d]' % i, x=r['FR'], eta=eta8, err=r['LHS'],
                        rel=r['relL1'], old=min(r['RHS_new'], triv), new=min(triv, b2), real=True, wins=b2 < triv))
    return out


def table(pts):
    """build_n3.py's w(): header 'x y', '%.6g %.6g' per point; LF here, written in text mode (platform line ends)"""
    return 'x y\n' + ''.join('%.6g %.6g\n' % p for p in pts)


def main(argv):
    live = [r for r in rows() if r['rel'] > NOISE]
    assert len(live) == 606, len(live)
    key = 'old' if '--old' in argv else 'new'
    for r in live:
        assert r[key] > 0.0, r
        r['ratio'] = r['err'] / r[key]
    files = {'n3_thmB_frac.dat': table([(r['x'], r['ratio']) for r in live]),
             'n3_thmB_eta.dat': table([(r['eta'], r['ratio']) for r in live])}
    if key == 'old':
        for fn, b in files.items():
            dep = io.open(os.path.join(OLD, fn), 'rb').read().replace(b'\r\n', b'\n')
            assert b.encode('ascii') == dep, '%s: the first-order rebuild differs from the deposited file' % fn
            print('make_fig_thmB_data --old: %s rebuilt with build_n3.py\'s certificate is byte-identical, line '
                  'ends normalised to LF, to _superseded/arxiv_v4/paper/figs/%s (%d bytes)' % (fn, fn, len(dep)))
        return 0
    pooled = {p['src']: p for p in load(POOLED)['rows']}
    assert len(pooled) == len(live)
    for r in live:
        q = pooled[r['src']]
        assert q['real'] == r['real'] and abs(r['ratio'] - q['ratio_new']) <= 1e-12 * q['ratio_new'], (r, q)
    n_real, n_cplx = sum(r['real'] for r in live), sum(not r['real'] for r in live)
    w_real = sum(r['wins'] for r in live if r['real'])
    w_cplx = sum(r['wins'] for r in live if not r['real'])
    summ = load(POOLED)['summary']
    assert (n_real, n_cplx, w_real, w_cplx) == (summ['n_real_live'], summ['n_complex_live'],
                                                summ['real_rows_B2_wins'], summ['complex_rows_B1_wins'])
    viol = sum(r['ratio'] > 1.0 + 1e-9 for r in live)
    assert viol == 0
    lo = min(live, key=lambda r: r['ratio'])
    hi = max(live, key=lambda r: r['ratio'])
    near = [r['src'] for r in live if hi['ratio'] - r['ratio'] < 1e-9]
    print('make_fig_thmB_data: %d live rows (%d real-probe, %d complex); violations %d; ratio %.17g (%s) to %.17g (%s)'
          % (len(live), n_real, n_cplx, viol, lo['ratio'], lo['src'], hi['ratio'], hi['src']))
    print('make_fig_thmB_data: leakage entry below the trivial cap on %d rows = %d real + %d complex; '
          'rows within 1e-9 of the maximum: %s' % (w_real + w_cplx, w_real, w_cplx, near))
    if '--dry' in argv:
        print('make_fig_thmB_data --dry: nothing written')
        return 0
    for fn, b in files.items():
        with io.open(os.path.join(OUT, fn), 'w', encoding='ascii') as fh:
            fh.write(b)
        back = io.open(os.path.join(OUT, fn), 'r', encoding='ascii').read()
        assert back == b
        print('make_fig_thmB_data: wrote %s (%d points, %d characters, platform line ends)'
              % (os.path.join(OUT, fn), len(live), len(b)))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
