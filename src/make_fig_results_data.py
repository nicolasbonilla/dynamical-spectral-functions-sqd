# -*- coding: utf-8 -*-
"""make_fig_results_data.py -- the data table of Fig. 1 of the manuscript prepared for JCP (fig:results).

2026-11-16, second-order branch of Theorem 1 (rules D14, D17, D22 of the change, docs/KNOWN_DISCREPANCIES.md
section 36).  The leakage curve is the SMALLER of the two leakage branches at each row (D17),
    B1 = (1+sqrt w)(sqrt(1-w)+x)   the first-order branch (stored: bound_leak, relR5),
    B2 = (sqrt(1-w)+x)^2           the second-order branch (a real probe: Lambda~ = Lambda),
with w = w_S and x = Lambda^/eta (stored: Lam_over_eta, Lamhat/eta), so that the certificate min(1+w, B2) is
the lower of the two plotted upper curves; B1 is the smaller up to FR 0.60, B2 from FR 0.70 on.  The crossing
FR* of the trivial bound is the re-scored one.  The only arithmetic is B2 and min(B1, B2) on stored fields;
every other value is COPIED.  Every plotted leakage value and FR* are asserted equal to the re-scoring outputs
data/c3_frontier/second_order/v2.json (D17.fig1) and cells.json (written by src/frontier/rescore_second_order.py).

  * Born-ranked (infinite-shot) subspaces, L = 8, eta = 0.18 t, sector fractions 0.20--0.95:
    data/c3_frontier/C3_grid.json, key 'grid', rows with L == 8 and eta == 0.18
    (columns FR_real, relL1 = true relative L1 error, bound_leak = B1, Lam_over_eta = x,
    bound_trivial = 1 + w_S, w_S).
  * The same ranking at |S| = 3802, 3830, 3855 (FR 0.970--0.983), reconstructed with the dense
    Rayleigh-Ritz spectrum: data/c3_frontier/adversarial/z_fine_L8.txt (relR0 = true
    relative L1 error, relR5 = B1, Lamhat/eta = x, w_S), parsed with the regular expression
    verify.py uses on that file, extended to Lamhat/eta.  With B2 these three points lie far
    below 0.05; the certificate's first fall below 0.05 is bracketed by the grid rows at FR 0.90
    and 0.92.
  * The probe's own determinant support S* = supp(phi), w_S = 1 exactly, L = 8, eta = 0.18 t:
    data/support_witness.json (results[L=8].frac, table[eta=0.18].relL1,
    bound_theoremB = B1, the smaller there; that column is not plotted).
  * The fraction at which the true error falls below 0.05, L = 8, eta = 0.18 t: C3_grid.json,
    key 'frontier' (FR_true['0.05'].FR); the crossing FR* of the certificate and the trivial
    bound: data/c3_frontier/second_order/cells.json (cells[L=8, eta=0.18].FRstar_new; the stored FR_star is the
    crossing of the first-order branch).

Writes  paper/figs/fig_results_L8.dat      (ranked subspaces, both sources; the columns
                                           branch, B1 and B2 record which one is plotted)
        paper/figs/fig_results_star.dat    (the one S* point)
        paper/figs/fig_results_marks.tex   (\\def's of FR* and FR_true(0.05))
and re-reads each file to assert that every written number equals its source to 1e-9.

Usage (from the repository root):  python src/make_fig_results_data.py

Deposit 3.5.0: the maker of the JCP package (2026-11-16) with its paths taken from its place in the repository and
the provenance comments it writes into the three files naming deposit paths; the numbers it writes are unchanged.
src/check_figures.py regenerates the three files in a scratch tree and compares them byte for byte.
"""
import io
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))                      # <repo>/src
REPO = os.path.dirname(HERE)                                           # <repo>
DATA = os.path.join(REPO, 'data')
OUT = os.path.join(REPO, 'paper', 'figs')
RESC = os.path.join(DATA, 'c3_frontier', 'second_order')              # the re-scored outputs
L, ETA = 8, 0.18


def b2(w, x):
    """second-order branch relative to ||phi||^2, real probe (Theorem 1, D1): (sqrt(1-w) + x)^2"""
    return (math.sqrt(max(1.0 - w, 0.0)) + x) ** 2


def load(rel):
    with io.open(os.path.join(DATA, rel), encoding='utf-8') as fh:
        return json.load(fh)


def main():
    grid = load('c3_frontier/C3_grid.json')
    rows = sorted((r for r in grid['grid'] if r['L'] == L and abs(r['eta'] - ETA) < 1e-12),
                  key=lambda r: r['FR_real'])
    assert len(rows) == 15, len(rows)
    tab = []
    for r in rows:
        w, x, b1 = r['w_S'], r['Lam_over_eta'], r['bound_leak']
        assert abs((1.0 + math.sqrt(w)) * (math.sqrt(max(1.0 - w, 0.0)) + x) - b1) <= 1e-12 * b1, r
        tab.append((r['FR_real'], r['relL1'], min(b1, b2(w, x)), r['bound_trivial'], 1.0 - w, 'grid',
                    'B1' if b1 <= b2(w, x) else 'B2', b1, b2(w, x)))
    with io.open(os.path.join(DATA, 'c3_frontier/adversarial/z_fine_L8.txt'), encoding='utf-8') as fh:
        zf = fh.read()
    fine = [(int(m.group(1)), float(m.group(2)), float(m.group(3)), float(m.group(4)), float(m.group(5)),
             float(m.group(6)))
            for m in re.finditer(r"\|S\|=(\d+) FR=([\d.]+) eta=0\.18 \| w_S=([\d.e+-]+) relR0=([\d.e+-]+) "
                                 r"relR5=([\d.e+-]+) Lamhat/eta=([\d.e+-]+)", zf)]
    assert [f[0] for f in fine] == [3802, 3830, 3855], fine
    for k, fr, w, r0, r5, x in fine:
        # relR5 is stored to 7 significant figures: it reproduces from (w, x) to 1e-5 relative
        assert abs((1.0 + math.sqrt(w)) * (math.sqrt(max(1.0 - w, 0.0)) + x) - r5) <= 1e-5 * r5, k
        tab.append((k / 3920.0, r0, min(r5, b2(w, x)), 1.0 + w, 1.0 - w, 'fine', 'B1' if r5 <= b2(w, x) else 'B2',
                    r5, b2(w, x)))
    front = [f for f in grid['frontier'] if f['L'] == L and abs(f['eta'] - ETA) < 1e-12]
    assert len(front) == 1
    fr_true05 = front[0]['FR_true']['0.05']['FR']
    # FR*: where the certificate min(1+w, B2) first meets the trivial bound (D6), as re-scored in
    # data/c3_frontier/second_order/cells.json; the stored FR_star (0.794049) is the crossing of the first-order
    # branch
    with io.open(os.path.join(RESC, 'cells.json'), encoding='utf-8') as fh:
        c8 = [c for c in json.load(fh)['cells'] if c['L'] == L and abs(c['eta'] - ETA) < 1e-12]
    assert len(c8) == 1 and abs(c8[0]['FRstar_old_stored'] - front[0]['FR_star']) <= 1e-15
    assert abs(c8[0]['FR05'] - fr_true05) <= 1e-15
    fr_star = c8[0]['FRstar_new']
    with io.open(os.path.join(RESC, 'v2.json'), encoding='utf-8') as fh:
        v2 = json.load(fh)['D17']['fig1']
    assert len(v2['rows']) == len(tab) == 18 and abs(v2['FRstar_new'] - fr_star) <= 1e-15
    for t, q in zip(tab, v2['rows']):
        assert abs(t[0] - q['FR']) <= 1e-12 and abs(t[2] - q['leak_D17']) <= 1e-11 * q['leak_D17'], (t, q)
        assert t[6] == {'first-order': 'B1', 'second-order': 'B2'}[q['smaller']], (t, q)
    sw = [r for r in load('support_witness.json')['results'] if r['L'] == L][0]
    st = [t for t in sw['table'] if abs(t['eta'] - ETA) < 1e-12][0]
    star = (sw['frac'], st['relL1'], st['bound_theoremB'], st['bound_trivial'])
    assert st['bound_theoremB'] <= b2(min(sw['w_S'], 1.0), st['Lambda_hat_over_eta'])   # D17: B1 is the smaller

    p_tab = os.path.join(OUT, 'fig_results_L8.dat')
    with io.open(p_tab, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('# Fig. 1 (fig:results): Born-ranked (infinite-shot) subspaces, L=8, eta=0.18t.\n')
        fh.write('# Written by src/make_fig_results_data.py from\n')
        fh.write('# data/c3_frontier/C3_grid.json (src=grid) and adversarial/z_fine_L8.txt (src=fine).\n')
        fh.write('# leak = min(B1, B2), the smaller leakage branch (2026-11-16, rule D17): B1 = the stored\n')
        fh.write('# first-order branch (bound_leak, relR5), B2 = (sqrt(1-w_S) + Lambda^/eta)^2; branch = the smaller.\n')
        fh.write('FR relL1 leak trivial missed src branch B1 B2\n')
        for t in tab:
            fh.write('%.10g %.10g %.10g %.10g %.10g %s %s %.10g %.10g\n' % t)
    p_star = os.path.join(OUT, 'fig_results_star.dat')
    with io.open(p_star, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('# S* = supp(phi), w_S = 1, L=8, eta=0.18t; data/support_witness.json\n')
        fh.write('FR relL1 leak trivial\n%.10g %.10g %.10g %.10g\n' % star)
    p_mk = os.path.join(OUT, 'fig_results_marks.tex')
    with io.open(p_mk, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('%% written by src/make_fig_results_data.py from data/ and data/c3_frontier/second_order/\n')
        fh.write('\\def\\FigResFRstar{%.6f}%% second_order/cells.json L=8 eta=0.18 FRstar_new (certificate)\n'
                 % fr_star)
        fh.write('\\def\\FigResFRtrue{%.6f}%% C3_grid.json frontier L=8 eta=0.18 FR_true[0.05]\n' % fr_true05)

    # round trip
    back = [l.split() for l in io.open(p_tab, encoding='utf-8').read().splitlines()
            if l and not l.startswith('#')][1:]
    assert len(back) == len(tab) == 18
    for b, t in zip(back, tab):
        assert b[5:7] == list(t[5:7]), (b, t)
        for x, y in zip(b[:5] + b[7:], t[:5] + t[7:]):
            assert abs(float(x) - y) <= 1e-9 * max(1.0, abs(y)), (x, y)
    print('fig_results: %d ranked points (15 grid + 3 fine), S* at FR=%.4f err=%.4f, FR*=%.4f, '
          'FR_true(0.05)=%.4f' % (len(tab), star[0], star[1], fr_star, fr_true05))
    print('fig_results: true error %.3e..%.3e; smaller leakage branch %.3g..%.3g; 1-w %.2e..%.2e'
          % (min(t[1] for t in tab), max(t[1] for t in tab), min(t[2] for t in tab), max(t[2] for t in tab),
             min(t[4] for t in tab), max(t[4] for t in tab)))
    print('fig_results: first-order branch the smaller up to FR %.2f, second-order from FR %.2f'
          % (max(t[0] for t in tab if t[6] == 'B1'), min(t[0] for t in tab if t[6] == 'B2')))
    return 0


if __name__ == '__main__':
    sys.exit(main())
