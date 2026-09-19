# -*- coding: utf-8 -*-
r"""C3 aggregation: the grid, the frontier curve with its uncertainty, the
convergence census, and the verdict on the published fraction.

Reads out/pt_L*_FR*_nl*.json ; writes out/C3_grid.json and prints the tables.
"""
import os, sys, glob, json, math, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
# --- deposit paths (repaired 2026-09-19; see docs/C3_FRONTIER.md) ----------
SRC = os.path.dirname(HERE)                       # release/src
ROOT = os.path.dirname(SRC)                       # repository root
OUT = os.environ.get('C3_OUT') or os.path.join(ROOT, 'data', 'c3_frontier')
TMP = os.environ.get('C3_TMP') or os.path.join(ROOT, 'build')
# --------------------------------------------------------------------------
ETAS = [0.10, 0.15, 0.18, 0.25]


def load():
    pts = {}
    for f in sorted(glob.glob(os.path.join(OUT, 'pt_L*_nl*.json'))):
        # The 204-row grid is the FRONTIER sweep of 2026-09-18.  The four
        # published-fraction points measured on 2026-09-19 live in
        # OUT/published_fraction/ and carry a _PUB tag; they are summarised
        # by pubsum.py and are deliberately NOT aggregated here, because
        # adding them would silently change the row count the manuscript
        # quotes.  This skip is what keeps `python src/frontier/fsum.py`
        # reproducing the deposited C3_grid.json exactly.
        if f.endswith('_PUB.json'):
            continue
        d = json.load(open(f))
        key = (d['meta']['L'], round(d['FR'], 4))
        # if several depths exist for the same (L,FR), keep the deepest
        if key in pts and pts[key]['nl'] >= d['nl']:
            continue
        pts[key] = d
    return pts


def top_rows(d):
    n = d['nl']
    return {round(r['eta'], 4): r for r in d['rows'] if r['n'] == n}


def plateau_of(d, eta, klast=3):
    """Robust plateau test on Lambda_K,out: max relative spread over the last
    `klast` ladder depths.  Lambda_out is NOT monotone at small n (measured), so a
    two-point test is not enough."""
    seq = sorted([(r['n'], r['Lambda_K_out'], r['bound_leak_outonly'])
                  for r in d['rows'] if abs(r['eta'] - eta) < 1e-9])
    if len(seq) < 2:
        return dict(rel_spread=float('nan'), n_used=[], Lout=float('nan'))
    tail = seq[-klast:] if len(seq) >= klast else seq
    vals = np.array([x[1] for x in tail])
    ref = vals[-1]
    return dict(rel_spread=float(np.abs(vals - ref).max() / max(abs(ref), 1e-300)),
                n_used=[x[0] for x in tail], Lout=float(ref),
                Bout_spread=float(np.abs(np.array([x[2] for x in tail])
                                         - tail[-1][2]).max()))


def conv_flag(r, plateau_relspread):
    """Convergence of the CERTIFICATE in Krylov depth.
    in/out is the ratio of the truncation part of Lambda to the boundary part.
    Calibrated at L=8 against the dense Lambda_S: in/out = 0.05 gave K/S-1 = 1.3e-3."""
    io = r['in_over_out']
    if io < 0.05 and plateau_relspread < 1e-3:
        return 'CONVERGED'
    if io < 0.30:
        return 'MARGINAL'
    return 'NOT-CONVERGED'


def frontier(Bs, FRs, thresh):
    """FR at which the leakage bound B crosses the trivial bound.
    Returns (FR*, lo, hi) with the bracket as the uncertainty. log-linear interp."""
    order = np.argsort(FRs)
    FRs = np.asarray(FRs)[order]; Bs = np.asarray(Bs)[order]
    th = np.asarray(thresh)[order]
    below = Bs < th
    if below.all():
        return (float(FRs[0]), float('nan'), float(FRs[0]), 'below-grid')
    if not below.any():
        return (float('nan'), float(FRs[-1]), float('nan'), 'above-grid')
    i = int(np.argmax(below))            # first index where B < trivial
    if i == 0:
        return (float(FRs[0]), float('nan'), float(FRs[0]), 'below-grid')
    x0, x1 = FRs[i - 1], FRs[i]
    y0 = math.log(Bs[i - 1] / th[i - 1]); y1 = math.log(Bs[i] / th[i])
    xs = x0 + (x1 - x0) * (0.0 - y0) / (y1 - y0)
    return (float(xs), float(x0), float(x1), 'bracketed')


def main():
    pts = load()
    Ls = sorted(set(k[0] for k in pts))
    grid = []
    for L in Ls:
        FRs = sorted(k[1] for k in pts if k[0] == L)
        for FR in FRs:
            d = pts[(L, FR)]
            tr = top_rows(d)
            for eta in ETAS:
                r = tr.get(round(eta, 4))
                if r is None:
                    continue
                pl = plateau_of(d, eta)
                grid.append(dict(
                    L=L, FR=FR, FR_real=d['FR_real'], S=d['S'], D=d['meta']['D'],
                    eta=eta, nl=d['nl'], w_S=d['w_S'],
                    Lam_over_eta=r['LamHat_over_eta'],
                    Lam_in=r['Lambda_K_in'], Lam_out=r['Lambda_K_out'],
                    Lam_out_over_eta=r['LamHat_out_over_eta'],
                    bound_leak=r['bound_leak'], bound_trivial=r['bound_trivial'],
                    cert=r['cert'], bound_leak_limit=r['bound_leak_outonly'],
                    cert_limit=r['cert_outonly'],
                    relL1=r.get('relL1'), relL1_win=r.get('relL1_win'),
                    slack=r.get('slack'), slack_leakbranch=r.get('slack_leakbranch'),
                    in_over_out=r['in_over_out'],
                    plateau_relspread=pl['rel_spread'],
                    plateau_depths=pl['n_used'],
                    Bout_spread=pl.get('Bout_spread'),
                    conv=conv_flag(r, pl['rel_spread']),
                    nontrivial=bool(r['bound_leak'] < r['bound_trivial']),
                    nontrivial_limit=bool(r['bound_leak_outonly'] < r['bound_trivial']),
                    ties_in=d['subspace_diag']['n_tied_at_cut_inside'],
                    ties_tot=d['subspace_diag']['n_tied_at_cut_total'],
                    zeroBorn=d['subspace_diag']['n_zero_wc_in_S'],
                    t_point=d['t_point'], peak_GiB=d['peak_vec_GiB']))

    # ---------------- frontier curve ----------------
    fr = []
    for L in Ls:
        for eta in ETAS:
            sel = [g for g in grid if g['L'] == L and abs(g['eta'] - eta) < 1e-9]
            # the specified grid is FR in {0.20..0.90}; refinement points added around
            # the crossing and below it are used too -- they only tighten the bracket.
            if len(sel) < 2:
                continue
            FRs = [g['FR'] for g in sel]
            xs, lo, hi, how = frontier([g['bound_leak'] for g in sel], FRs,
                                       [g['bound_trivial'] for g in sel])
            xl, ll, hl, howl = frontier([g['bound_leak_limit'] for g in sel], FRs,
                                        [g['bound_trivial'] for g in sel])
            nconv = sum(1 for g in sel if g['conv'] == 'CONVERGED')
            # where does the TRUE error cross the repo's own threshold?  This is the
            # answer to "is the frontier physics or is it the Cauchy-Schwarz step?"
            sel2 = [g for g in sel if g['relL1'] is not None]
            true_thr = {}
            for thr in (0.05, 0.01):
                xs2, lo2, hi2, how2 = frontier([g['relL1'] for g in sel2],
                                               [g['FR'] for g in sel2],
                                               [thr] * len(sel2))
                true_thr[str(thr)] = dict(FR=xs2, lo=lo2, hi=hi2, how=how2)
            # slack of the certificate where it first becomes non-trivial
            nt = [g for g in sel if g['nontrivial'] and g['slack'] is not None]
            slack_at_frontier = min((g['slack'] for g in nt), default=float('nan'))
            fr.append(dict(L=L, eta=eta, FR_star=xs, bracket_lo=lo, bracket_hi=hi, how=how,
                           FR_star_limit=xl, bracket_limit=(ll, hl), how_limit=howl,
                           FR_true=true_thr, slack_at_frontier=slack_at_frontier,
                           n_points=len(sel), n_converged=nconv))

    json.dump(dict(grid=grid, frontier=fr), open(os.path.join(OUT, 'C3_grid.json'), 'w'),
              indent=1, default=float)

    # ---------------- tables ----------------
    print("\n" + "=" * 118)
    print("C3 GRID -- Theorem B certificate vs measured error.  Protocol K=18,dt=0.5,sub=16 "
          "(repo).  Full reorthogonalisation.")
    print("cert = min{leak bound, trivial 1+w_S}.  relL1 = ||A - A_K||_1/||phi||^2 (whole line).")
    print("=" * 118)
    hdr = ("  L    FR   |S|/D    w_S        eta   Lam/eta   B_leak  triv   cert    relL1     "
           "slack   in/out    conv")
    for L in Ls:
        print("\n" + hdr)
        print("  " + "-" * 114)
        for g in sorted([x for x in grid if x['L'] == L], key=lambda x: (x['FR'], x['eta'])):
            print("%3d  %.2f  %.3f  %.8f  %.2f  %8.4f %8.4f %5.3f %7.4f %.3e %8.1f %.2e  %s%s"
                  % (g['L'], g['FR'], g['FR_real'], g['w_S'], g['eta'], g['Lam_over_eta'],
                     g['bound_leak'], g['bound_trivial'], g['cert'],
                     g['relL1'] if g['relL1'] is not None else float('nan'),
                     g['slack'] if g['slack'] is not None else float('nan'),
                     g['in_over_out'], g['conv'],
                     "  *CERT*" if g['nontrivial'] else ""))

    print("\n" + "=" * 96)
    print("FRONTIER  FR*(L,eta): the fraction at which the certificate stops being trivial.")
    print("  measured  = from Lambda_K at the depth actually run (certifies A_K at that depth)")
    print("  limit     = from Lambda_K,out alone (the n_l -> inf extrapolation; NOT a bound)")
    print("=" * 96)
    print("  L   eta    FR*      bracket         +-grid   +-depth   FR*_limit   conv/pts")
    for f in sorted(fr, key=lambda x: (x['L'], x['eta'])):
        dg = (f['bracket_hi'] - f['bracket_lo']) / 2 if f['bracket_hi'] == f['bracket_hi']             and f['bracket_lo'] == f['bracket_lo'] else float('nan')
        dd = abs(f['FR_star'] - f['FR_star_limit'])             if f['FR_star'] == f['FR_star'] and f['FR_star_limit'] == f['FR_star_limit'] else float('nan')
        f['unc_grid'] = dg; f['unc_depth'] = dd
        print("%3d  %.2f   %6s   [%.2f,%.2f]   %7.3f   %7.4f   %9s   %d/%d"
              % (f['L'], f['eta'],
                 ("%.3f" % f['FR_star']) if f['FR_star'] == f['FR_star'] else ">grid",
                 f['bracket_lo'] if f['bracket_lo'] == f['bracket_lo'] else float('nan'),
                 f['bracket_hi'] if f['bracket_hi'] == f['bracket_hi'] else float('nan'),
                 dg, dd,
                 ("%.3f" % f['FR_star_limit']) if f['FR_star_limit'] == f['FR_star_limit'] else ">grid",
                 f['n_converged'], f['n_points']))

    print("\n" + "=" * 104)
    print("IS THE FRONTIER PHYSICS, OR IS IT THE BOUND?  (the warning of VE section 11)")
    print("  FR*        : where the CERTIFICATE stops being trivial")
    print("  FR(0.05)   : where the TRUE rel-L1 crosses the repo's own 5% threshold")
    print("  gap        : how much extra sector the certificate demands over the truth")
    print("=" * 104)
    print("  L   eta     FR*      FR(relL1=0.05)   FR(relL1=0.01)   gap(FR*-FR.05)   "
          "slack at frontier")
    for f in sorted(fr, key=lambda x: (x['L'], x['eta'])):
        a = f['FR_star']; b = f['FR_true']['0.05']['FR']; c2 = f['FR_true']['0.01']['FR']
        print("%3d  %.2f   %6s   %14s   %14s   %14s   %s"
              % (f['L'], f['eta'],
                 ("%.3f" % a) if a == a else ">0.97",
                 ("%.3f" % b) if b == b else "n/a",
                 ("%.3f" % c2) if c2 == c2 else "n/a",
                 ("%+.3f" % (a - b)) if (a == a and b == b) else "n/a",
                 ("%.0fx" % f['slack_at_frontier'])
                 if f['slack_at_frontier'] == f['slack_at_frontier'] else "n/a"))

    # ---------------- convergence census ----------------
    print("\n" + "=" * 70)
    print("CONVERGENCE CENSUS")
    print("=" * 70)
    from collections import Counter
    c = Counter((g['L'], g['conv']) for g in grid)
    for L in Ls:
        tot = sum(v for (l, _), v in c.items() if l == L)
        print("  L=%2d  total=%3d   " % (L, tot)
              + "  ".join("%s=%d" % (k, c[(L, k)]) for k in
                          ('CONVERGED', 'MARGINAL', 'NOT-CONVERGED') if c[(L, k)]))
    bad = [g for g in grid if g['conv'] != 'CONVERGED']
    if bad:
        print("\n  NOT fully converged (L, FR, eta, in/out, nl):")
        for g in bad:
            print("    L=%2d FR=%.2f eta=%.2f  in/out=%.3e  nl=%d  -> %s "
                  "(B=%.4f vs limit B_out=%.4f, trivial=%.4f)"
                  % (g['L'], g['FR'], g['eta'], g['in_over_out'], g['nl'], g['conv'],
                     g['bound_leak'], g['bound_leak_limit'], g['bound_trivial']))

    # ---------------- published fraction ----------------
    print("\n" + "=" * 70)
    print("VERDICT ON THE PUBLISHED FRACTION")
    print("=" * 70)
    pub = {4: 0.92, 6: 0.82, 8: 0.56, 10: 0.36, 12: 0.18, 14: 0.08}
    for L in Ls:
        for g in grid:
            if g['L'] == L and abs(g['FR'] - pub.get(L, -1)) < 5e-3:
                print("  L=%2d FR=%.2f eta=%.2f: B=%.4f vs trivial %.4f -> %s by x%.2f  "
                      "(relL1=%.3e, in/out=%.1e, %s)"
                      % (L, g['FR'], g['eta'], g['bound_leak'], g['bound_trivial'],
                         "NON-TRIVIAL" if g['nontrivial'] else "VACUOUS",
                         g['bound_leak'] / g['bound_trivial'],
                         g['relL1'] if g['relL1'] else float('nan'),
                         g['in_over_out'], g['conv']))

    # ---------------- cost ----------------
    print("\n" + "=" * 70)
    print("MEASURED COST")
    print("=" * 70)
    for L in Ls:
        sel = [pts[k] for k in pts if k[0] == L]
        tt = sum(d['t_point'] for d in sel)
        print("  L=%2d  %2d points, %7.1f s of subspace work (%.2f h), peak %.2f GiB/point, "
              "nl=%s" % (L, len(sel), tt, tt / 3600,
                         max(d['peak_vec_GiB'] for d in sel),
                         sorted(set(d['nl'] for d in sel))))


if __name__ == '__main__':
    main()
