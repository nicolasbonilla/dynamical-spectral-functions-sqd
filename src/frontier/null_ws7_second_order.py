# -*- coding: utf-8 -*-
"""Robustness of the out-of-sample null test across every 2-vs-2 and 1-vs-3 split.

null_ws7_second_order.py -- COPY of src/frontier/null_ws7.py (rule D20 of the change to the second-order branch,
2026-11-16; docs/KNOWN_DISCREPANCIES.md section 36).  The manuscript prepared for JCP names it in Secs. V D and S4.
The computation is the deposit's, line for line.  Differences, each marked in the code:
  (1) CHANGED: paths. The grid is read from data/c3_frontier (or $C3_OUT, as the original allows).
  (2) CHANGED: FR* per cell is the RE-SCORED crossing of the certificate, cells[].FRstar_new of
      data/c3_frontier/second_order/cells.json (written by src/frontier/rescore_second_order.py), instead of the
      stored first-order FR_star.  With --stored the stored FR_star is used, as in the original (control: the
      printed table must then be the original's, line for line).
  (3) ADDED (record only, no effect on any number): the best (spread, c, p) of each fitted rule and its firing
      level c*eta^p are kept, and everything is written to data/c3_frontier/second_order/null_ws7_second_order.json
      (or, with --stored, null_ws7_first_order_copy.json there).
Arithmetic and a small grid search over the two constants on stored numbers; no simulation.  Runs in about 20 s on
one thread.  Usage, from the repository root:
    python src/frontier/null_ws7_second_order.py > data/c3_frontier/second_order/null_ws7_second_order.txt
    python src/frontier/null_ws7_second_order.py --stored > data/c3_frontier/second_order/null_ws7_first_order_copy.txt
Deposit 3.5.0: the script of the JCP package (2026-11-16) with its paths taken from its place in the repository
(and the two source strings it records naming deposit paths); nothing else changed.
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):   # ADDED: one thread
    os.environ.setdefault(_v, '1')
import json, math, os, itertools, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
# --- (1) CHANGED: paths --------------------------------------------------------------------------------------------
REPO = os.path.normpath(os.path.join(HERE, '..', '..'))                # <repo> (this file is <repo>/src/frontier)
OUT = os.environ.get('C3_OUT') or os.path.join(REPO, 'data', 'c3_frontier')
RES = os.path.join(REPO, 'data', 'c3_frontier', 'second_order')
STORED = '--stored' in sys.argv[1:]
# --------------------------------------------------------------------------------------------------------------------
G = json.load(open(os.path.join(OUT, 'C3_grid.json')))
grid, front = G['grid'], G['frontier']
cells = [(f['L'], f['eta']) for f in front]
# --- (2) CHANGED: the crossing fraction FR* of each cell -------------------------------------------------------------
if STORED:
    FRS = {(f['L'], f['eta']): f['FR_star'] for f in front}
    FRS_SRC = 'data/c3_frontier/C3_grid.json frontier[].FR_star (stored, first-order branch)'
else:
    _C = json.load(open(os.path.join(RES, 'cells.json'), encoding='utf-8'))['cells']
    FRS = {(c['L'], c['eta']): c['FRstar_new'] for c in _C}
    FRS_SRC = 'data/c3_frontier/second_order/cells.json cells[].FRstar_new (re-scored, certificate of D6 = second-order branch)'
assert sorted(FRS) == sorted(cells)
# --------------------------------------------------------------------------------------------------------------------
def series(L, eta):
    rows = sorted([r for r in grid if r['L'] == L and abs(r['eta']-eta) < 1e-9
                   and r['conv'] == 'CONVERGED' and r['relL1'] and r['relL1'] > 0],
                  key=lambda r: r['FR'])
    return (np.array([r['FR'] for r in rows]), np.array([1.0-r['w_S'] for r in rows]),
            np.array([r['relL1'] for r in rows]))
S = {c: series(*c) for c in cells}
def err_at(c, lvl):
    fr, om, rl = S[c]
    if lvl > om[0] or lvl < om[-1]: return None
    return math.exp(float(np.interp(float(np.interp(-math.log(lvl), -np.log(om), fr)), fr, np.log(rl))))
def cert(c):
    fr, om, rl = S[c]
    return math.exp(float(np.interp(FRS[c], fr, np.log(rl))))
sizes = [6, 8, 10, 12]
print("%-14s %-14s %-9s %-9s %-9s %s" % ("TRAIN L", "TEST L", "cert", "nullA", "nullB", "winner"))
tally = {'cert': 0, 'null': 0, 'na': 0}
RECORD = []                                                               # (3) ADDED: record only
for k in (1, 2, 3):
    for tr_s in itertools.combinations(sizes, k):
        te_s = [s for s in sizes if s not in tr_s]
        TRAIN = [c for c in cells if c[0] in tr_s]; TEST = [c for c in cells if c[0] in te_s]
        ct_tr = np.array([cert(c) for c in TRAIN]); ct_te = np.array([cert(c) for c in TEST])
        m_tr = float(np.median(ct_tr)); cs = ct_te.max()/ct_te.min()
        res = {}
        fit = {}                                                          # (3) ADDED: record only
        for name, ps in (('A', [None]), ('B', list(np.linspace(0, 5, 101)))):
            best = None
            for p in ps:
                for lg in np.arange(-6.0, -1.5, 0.02):
                    c0 = 10.0**lg
                    tr = [err_at(c, c0*(1.0 if p is None else c[1]**p)) for c in TRAIN]
                    if any(v is None for v in tr): continue
                    tr = np.array(tr); med = float(np.median(tr))
                    if not (m_tr/2 <= med <= m_tr*2): continue
                    sc = tr.max()/tr.min()
                    if best is None or sc < best[0]: best = (sc, c0, p)
            fit[name] = best                                              # (3) ADDED: record only
            if best is None: res[name] = None; continue
            te = [err_at(c, best[1]*(1.0 if best[2] is None else c[1]**best[2])) for c in TEST]
            res[name] = None if any(v is None for v in te) else (np.max(te)/np.min(te))
        vals = [v for v in res.values() if v is not None]
        win = 'cert' if (vals and cs < min(vals)) else ('null' if vals else 'n/a')
        tally['cert' if win == 'cert' else ('null' if win == 'null' else 'na')] += 1
        print("%-14s %-14s x%-8.2f %-9s %-9s %s" %
              (tr_s, tuple(te_s), cs,
               ("x%.2f" % res['A']) if res['A'] else "unreach",
               ("x%.2f" % res['B']) if res['B'] else "unreach", win.upper()))
        # (3) ADDED: record only -----------------------------------------------------------------------------------
        RECORD.append(dict(
            train=list(tr_s), test=list(te_s), cert_test_spread=float(cs), cert_errors_test=[float(v) for v in ct_te],
            nullA_test_spread=None if res['A'] is None else float(res['A']),
            nullB_test_spread=None if res['B'] is None else float(res['B']), winner=win,
            fitA=None if fit['A'] is None else dict(train_spread=float(fit['A'][0]), c=float(fit['A'][1]),
                                                     firing_level=float(fit['A'][1])),
            fitB=None if fit['B'] is None else dict(train_spread=float(fit['B'][0]), c=float(fit['B'][1]), p=float(fit['B'][2]),
                                                     firing_level_at_eta={('%g' % e): float(fit['B'][1] * e ** fit['B'][2])
                                                                          for e in (0.10, 0.15, 0.18, 0.25)})))
        # ---------------------------------------------------------------------------------------------------------
print("\nsplits where the ZERO-parameter certificate beats BOTH fitted w_S rules: %d / %d"
      % (tally['cert'], sum(tally.values())))
# --- (3) ADDED: summary and file (record only) -------------------------------------------------------------------------
n = len(RECORD)
lvB = [v for r in RECORD if r['fitB'] for v in r['fitB']['firing_level_at_eta'].values()]
lvA = [r['fitA']['firing_level'] for r in RECORD if r['fitA']]
SUM = dict(
    n_splits=n, FRstar_source=FRS_SRC, FRstar={'L%d_eta%g' % c: float(FRS[c]) for c in cells},
    cert_beats_both=tally['cert'], null_wins=tally['null'], na=tally['na'],
    two_constant_rule_tighter_than_cert=sum(1 for r in RECORD if r['nullB_test_spread'] is not None
                                            and r['nullB_test_spread'] < r['cert_test_spread']),
    eta_blind_rule_looser_than_cert=sum(1 for r in RECORD if r['nullA_test_spread'] is None
                                        or r['nullA_test_spread'] > r['cert_test_spread']),
    two_constant_firing_level_range=[min(lvB), max(lvB)] if lvB else None,
    two_constant_firing_level_range_per_eta={('%g' % e): [min(r['fitB']['firing_level_at_eta']['%g' % e] for r in RECORD if r['fitB']),
                                                          max(r['fitB']['firing_level_at_eta']['%g' % e] for r in RECORD if r['fitB'])]
                                             for e in (0.10, 0.15, 0.18, 0.25)} if lvB else None,
    eta_blind_firing_level_range=[min(lvA), max(lvA)] if lvA else None,
    floor_5e3_over_two_constant_levels=[5e-3 / max(lvB), 5e-3 / min(lvB)] if lvB else None)
print("two-constant rule tighter than the certificate on %d / %d splits; eta-blind rule looser on %d / %d"
      % (SUM['two_constant_rule_tighter_than_cert'], n, SUM['eta_blind_rule_looser_than_cert'], n))
if lvB:
    print("two-constant rule fires at 1-w_S = %.3e .. %.3e (all splits, eta = 0.10-0.25); floor 5e-3 is %.1f .. %.1f times larger"
          % (min(lvB), max(lvB), 5e-3 / max(lvB), 5e-3 / min(lvB)))
with open(os.path.join(RES, 'null_ws7_first_order_copy.json' if STORED else 'null_ws7_second_order.json'), 'w',
          encoding='utf-8') as fh:
    json.dump(dict(summary=SUM, splits=RECORD), fh, indent=1)
