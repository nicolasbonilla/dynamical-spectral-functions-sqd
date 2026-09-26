# -*- coding: utf-8 -*-
"""N3 appendix figure: violation ratio of the published Theorem 1(iii) bound vs. the
replacement (Theorem B), as a function of sector fraction and of eta.
EVERY number is read from JSON. Nothing is typed by hand.

DEPOSITED 2026-09-26 (src/recovered/README.md).  Written on 2026-09-18 in the author's session
scratch directory, where it wrote the ten paper/figs/n3_*.dat tables of Fig. S5
(fig:thm1iii-violation), and recovered from there.  Changes against the recovered copy, and
nothing else:
  * paths, derived from __file__: the two thm3 JSONs and the three certificate scans are read
    from data/thm1iii_violation/ (the scans from its cert_earlier_run/ -- the EARLIER run of
    the same scans that the plotted pool was built from, NOT the deposited data/cert_*.json;
    see that folder's README), thm3_L8b.py from src/recovered/, and the .dat tables are
    written to paper/figs/, where they are committed; n3_summary.json goes to build/;
  * the legend self-check at the end: two of its five search keys were updated to the legend
    wording the committed fragment now prints -- the two entries cite thm:leakage and
    fig:akwsampled by reference where they used to hard-code "Theorem~B" and "Fig.~5"; the
    legend was reworded by hand after this script ran.  The check itself (legend count ==
    .dat line count == builder count) is unchanged.
Run on those inputs it regenerates all ten n3_*.dat byte-identically (src/check_figures.py
guards this) and prints 621 pooled / 606 live / 468 violations / 363 of 363 at w_S >= 0.99 /
mildest factor 2.10.  Pointed at the deposited data/cert_*.json instead it gives 469 and 2.005,
the discrepancy the paper states in Sec. S9.
"""
import json, os, sys, re, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))          # <repo>/src/recovered
REPO = os.path.dirname(os.path.dirname(HERE))              # <repo>
SP = os.path.join(REPO, "data", "thm1iii_violation")       # thm3_results.json, thm3_L8b.json
OUT = os.path.join(REPO, "paper", "figs")                  # the ten n3_*.dat, and the fragment
BUILD = os.path.join(REPO, "build")                        # n3_summary.json (not deposited)
if not os.path.isdir(BUILD):
    os.makedirs(BUILD)
NOISE = 1e-10          # relative-L1 floor: below this the reconstruction is exact to machine noise

rows = []

# ---------- source 1: the original sweep (thm3_results.json) ----------
R = json.load(open(os.path.join(SP, "thm3_results.json")))
for r in R:
    rows.append(dict(src="thm3_results.json", family=r["tag"], sel=r["kind"],
                     frac=r["frac"], eta=r["eta"], wS=r["wS"],
                     err=r["LHS"], relerr=r["relL1"],
                     pub=r["RHS_pap"],
                     thmB=min(r["RHS_new"], r["RHS_triv"]),
                     leak_wins=bool(r["RHS_new"] < r["RHS_triv"])))
n_thm3 = len(R)

# ---------- source 2: flagship configuration L=8, eta=0.18 (thm3_L8b.json) ----------
L8 = json.load(open(os.path.join(SP, "thm3_L8b.json")))
sr = open(os.path.join(HERE, "thm3_L8b.py"), encoding="utf-8").read()
m = re.search(r"^L,U,eta\s*=\s*([0-9.]+)\s*,\s*([0-9.]+)\s*,\s*([0-9.]+)\s*$", sr, re.M)
assert m, "could not read L,U,eta from thm3_L8b.py"
L8L, L8U, L8eta = int(m.group(1)), float(m.group(2)), float(m.group(3))
flag = []
for r in L8:
    nphi2 = r["LB"] / (1.0 - r["wS"])          # LB = ||phi||^2 (1-w_S)
    triv = nphi2 * (1.0 + r["wS"])
    flag.append(dict(src="thm3_L8b.json", family="Hub L=%d U=%.1f" % (L8L, L8U),
                     sel="topTE", frac=r["FR"], eta=L8eta, wS=r["wS"],
                     err=r["LHS"], relerr=r["relL1"], pub=r["RHS_pap"],
                     thmB=min(r["RHS_new"], triv),
                     leak_wins=bool(r["RHS_new"] < triv)))

# ---------- sources 3-5: today's independent certificate run (certificate.py) ----------
CERT = os.path.join(SP, "cert_earlier_run")
cert_counts = {}
for f, famkey, selkey in [("cert_stress.json", "family", "selection"),
                          ("cert_akw.json", "branch", "FR"),
                          ("cert_teqsci.json", "seed", "K")]:
    d = json.load(open(os.path.join(CERT, f)))
    cert_counts[f] = len(d["rows"])
    for r in d["rows"]:
        rows.append(dict(src=f, family=str(r.get(famkey, f)), sel=str(r.get(selkey, "")),
                         L=r.get("L"), branch=r.get("branch"), kpi=r.get("k_over_pi"),
                         frac=r["frac"], eta=r["eta"], wS=r["w_S"],
                         err=r["L1_true"], relerr=r["relL1_true"],
                         pub=2.0 * r["lower"],          # 2||phi||^2 (1-w_S)
                         thmB=r["upper"],               # min(trivial, leakage), as stated
                         leak_wins=bool(r["bound_nontrivial"])))

ALL = rows + flag

# ---------- classification ----------
deg  = [r for r in ALL if r["relerr"] <= NOISE]
live = [r for r in ALL if r["relerr"] > NOISE]
infp = [r for r in live if r["pub"] <= 0.0]
finp = [r for r in live if r["pub"] > 0.0]
for r in finp:
    r["ratio_pub"] = r["err"] / r["pub"]
for r in live:
    r["ratio_B"] = r["err"] / r["thmB"] if r["thmB"] > 0 else float("inf")

viol_fin = [r for r in finp if r["ratio_pub"] > 1.0]
violB = [r for r in live if r["ratio_B"] > 1.0]


def sec(t):
    print("\n" + "=" * 72 + "\n" + t + "\n" + "=" * 72)


sec("PROVENANCE")
print("thm3_results.json          %4d rows" % n_thm3)
print("thm3_L8b.json              %4d rows   (flagship L=%d, U=%.1f, eta=%.2f, read from thm3_L8b.py)"
      % (len(L8), L8L, L8U, L8eta))
for k, v in cert_counts.items():
    print("%-26s %4d rows" % (k, v))
print("TOTAL pooled               %4d rows" % len(ALL))

sec("CLASSIFICATION")
print("exact to machine noise (relL1 <= %g), no ratio defined : %d" % (NOISE, len(deg)))
print("live subspaces                                          : %d" % len(live))
print("  published bound identically ZERO (w_S=1 numerically)  : %d   -> ratio = +inf" % len(infp))
print("  published bound > 0, finite ratio                     : %d" % len(finp))
print("     of which violate (ratio > 1)                       : %d" % len(viol_fin))
print("TOTAL violations of Theorem 1(iii)                      : %d / %d live" % (len(infp) + len(viol_fin), len(live)))
print("Theorem B violations                                    : %d / %d live" % (len(violB), len(live)))
print("Theorem B ratio range: %.4g .. %.6g" % (min(r["ratio_B"] for r in live), max(r["ratio_B"] for r in live)))
print("finite published ratio range: %.4g .. %.6g" % (min(r["ratio_pub"] for r in finp), max(r["ratio_pub"] for r in finp)))

sec("THE REGIME THE PAPER INVOKES (w_S >= 0.99)")
hi = [r for r in live if r["wS"] >= 0.99]
hiv = [r for r in hi if (r["pub"] <= 0) or (r["err"] / r["pub"] > 1)]
print("live subspaces with w_S >= 0.99 : %d" % len(hi))
print("   violating Theorem 1(iii)     : %d  (%.1f%%)" % (len(hiv), 100.0 * len(hiv) / len(hi)))
print("   violating Theorem B          : %d" % sum(1 for r in hi if r["ratio_B"] > 1))
print("   min violation factor there   : %.4g" % min((r["err"] / r["pub"]) for r in hi if r["pub"] > 0))

sec("THE PAPER'S OWN FIGURE-5 CONFIGURATION (cert_akw.json, L=8, eta=0.18t, momentum seeds)")
f5 = [r for r in finp if r["src"] == "cert_akw.json" and r.get("L") == 8]
print("  n = %d  (4 fractions x 8 momenta x 2 branches)" % len(f5))
for fr in sorted(set(r["frac"] for r in f5)):
    g = [r["ratio_pub"] for r in f5 if r["frac"] == fr]
    gb = [r["ratio_B"] for r in f5 if r["frac"] == fr]
    ge = [r["relerr"] for r in f5 if r["frac"] == fr]
    print("   FR=%.2f  n=%2d  violation ratio %.2f .. %.2f   relL1 %.3e .. %.3e   ThmB ratio %.4f .. %.4f"
          % (fr, len(g), min(g), max(g), min(ge), max(ge), min(gb), max(gb)))
print("  overall violation ratio: %.2f .. %.2f ; all > 1: %s"
      % (min(r["ratio_pub"] for r in f5), max(r["ratio_pub"] for r in f5),
         all(r["ratio_pub"] > 1 for r in f5)))

sec("SITE-SEEDED L=8 RUN (thm3_L8b.json) -- NOT the paper's Fig-5 configuration")
print("  thm3_L8b.py docstring states: 'Site seed instead of momentum seed'.")
print("  release/src/sampled_akw.py (Fig. 5) seeds with c^dag_{k,up}|GS>.  Different observable.")
for r in flag:
    print("  FR=%.2f  relL1=%.4e  pub bound=%.4e  ratio=%.4f   ThmB ratio=%.4g"
          % (r["frac"], r["relerr"], r["pub"], r["err"] / r["pub"], r["ratio_B"]))

sec("ETA COVERAGE")
for e in sorted(set(round(r["eta"], 4) for r in live)):
    rr = [r for r in live if abs(r["eta"] - e) < 1e-9]
    v = [r for r in rr if (r["pub"] <= 0) or (r["err"] / r["pub"] > 1)]
    fr = [r["err"] / r["pub"] for r in rr if r["pub"] > 0]
    print("  eta=%.2f  n=%4d  viol=%4d (%5.1f%%)  median finite ratio=%8.3g  max ThmB ratio=%.4g"
          % (e, len(rr), len(v), 100.0 * len(v) / len(rr),
             (np.median(fr) if fr else float("nan")), max(r["ratio_B"] for r in rr)))

sec("FRACTION COVERAGE")
print("  frac range: %.3f .. %.3f" % (min(r["frac"] for r in live), max(r["frac"] for r in live)))
print("  leakage branch beats trivial cap in %d of %d live rows" % (sum(1 for r in live if r["leak_wins"]), len(live)))

sec("EXACTLY-ATTAINED THEOREM-B ROWS (ratio_B == 1)")
sat = [r for r in live if abs(r["ratio_B"] - 1.0) < 1e-9]
print("  n = %d ; all have w_S = %s ; leakage branch used: %s"
      % (len(sat), sorted(set("%.3e" % r["wS"] for r in sat)),
         sorted(set(str(r["leak_wins"]) for r in sat))))
print("  -> these are subspaces ORTHOGONAL to phi: A_S == 0, so ||A-A_S||_1 = ||phi||^2 = the")
print("     trivial cap ||phi||^2(1+w_S). Exact algebraic identity on the branch that says nothing.")

# ---------- .dat tables for pgfplots ----------


def w(fn, pts):
    with open(os.path.join(OUT, fn), "w") as fh:
        fh.write("x y\n")
        for x, y in pts:
            fh.write("%.6g %.6g\n" % (x, y))
    return len(pts)


HIW = 0.99
pub_hi = [r for r in finp if r["wS"] >= HIW]
pub_lo = [r for r in finp if r["wS"] < HIW]
inf_hi = [r for r in infp if r["wS"] >= HIW]
inf_lo = [r for r in infp if r["wS"] < HIW]
HUB = ("thm3_results.json", "thm3_L8b.json", "cert_akw.json", "cert_teqsci.json")

counts = {}
counts["pub_frac_hi"] = w("n3_pub_frac_hi.dat", [(r["frac"], r["ratio_pub"]) for r in pub_hi])
counts["pub_frac_lo"] = w("n3_pub_frac_lo.dat", [(r["frac"], r["ratio_pub"]) for r in pub_lo])
counts["pub_eta_hi"] = w("n3_pub_eta_hi.dat", [(r["eta"], r["ratio_pub"]) for r in pub_hi])
counts["pub_eta_lo"] = w("n3_pub_eta_lo.dat", [(r["eta"], r["ratio_pub"]) for r in pub_lo])
counts["thmB_frac"] = w("n3_thmB_frac.dat", [(r["frac"], r["ratio_B"]) for r in live])
counts["thmB_eta"] = w("n3_thmB_eta.dat", [(r["eta"], r["ratio_B"]) for r in live])
counts["inf_frac"] = w("n3_inf_frac.dat", [(r["frac"], 1.0) for r in infp])
counts["inf_eta"] = w("n3_inf_eta.dat", [(r["eta"], 1.0) for r in infp])
fig5 = [r for r in finp if r["src"] == "cert_akw.json" and r.get("L") == 8]
counts["fig5_frac"] = w("n3_fig5_frac.dat", [(r["frac"], r["ratio_pub"]) for r in fig5])
counts["fig5_eta"] = w("n3_fig5_eta.dat", [(r["eta"], r["ratio_pub"]) for r in fig5])

sec("SPLIT BY CAPTURED WEIGHT (w_S >= %.2f)" % HIW)
print("  finite published ratio, w_S >= %.2f : %d   (all > 1? %s ; min %.4g)"
      % (HIW, len(pub_hi), all(r["ratio_pub"] > 1 for r in pub_hi),
         min(r["ratio_pub"] for r in pub_hi)))
print("  finite published ratio, w_S <  %.2f : %d   (of which > 1: %d)"
      % (HIW, len(pub_lo), sum(1 for r in pub_lo if r["ratio_pub"] > 1)))
print("  bound-identically-zero rows: %d high-weight, %d low-weight" % (len(inf_hi), len(inf_lo)))
print("  Hubbard-chain live rows: %d ; stress-family live rows: %d"
      % (sum(1 for r in live if r["src"] in HUB), sum(1 for r in live if r["src"] not in HUB)))
for e in sorted(set(round(r["eta"], 4) for r in live)):
    rr = [r for r in live if abs(r["eta"] - e) < 1e-9]
    print("     eta=%.2f : Hubbard %3d / stress %3d" % (e, sum(1 for r in rr if r["src"] in HUB),
                                                       sum(1 for r in rr if r["src"] not in HUB)))

sec("DAT FILES WRITTEN")
for k, v in counts.items():
    print("  %-12s %4d points" % (k, v))


json.dump(dict(n_pooled=len(ALL), n_live=len(live), n_deg=len(deg), n_inf=len(infp),
               n_fin=len(finp), n_viol=len(infp) + len(viol_fin), n_violB=len(violB),
               n_hi=len(hi), n_hiv=len(hiv),
               thmB_max=max(r["ratio_B"] for r in live),
               thmB_min=min(r["ratio_B"] for r in live),
               pub_min=min(r["ratio_pub"] for r in finp),
               pub_max=max(r["ratio_pub"] for r in finp),
               etas=sorted(set(round(r["eta"], 4) for r in live)),
               frac_min=min(r["frac"] for r in live), frac_max=max(r["frac"] for r in live),
               flagship=[dict(FR=r["frac"], ratio=r["err"] / r["pub"], relL1=r["relerr"],
                              ratioB=r["ratio_B"]) for r in flag],
               sources=dict(thm3_results=n_thm3, thm3_L8b=len(L8), **cert_counts),
               dat_counts=counts, hiw=HIW,
               n_pub_hi=len(pub_hi), n_pub_lo=len(pub_lo),
               n_inf_hi=len(inf_hi), n_inf_lo=len(inf_lo),
               pub_hi_min=min(r["ratio_pub"] for r in pub_hi),
               n_sat=len(sat),
               n_hub=sum(1 for r in live if r["src"] in HUB),
               n_stress=sum(1 for r in live if r["src"] not in HUB)),
          open(os.path.join(BUILD, "n3_summary.json"), "w"), indent=1)
print("\nsummary -> n3_summary.json")

# ---------- self-check: the counts written into the .tex legend must match the .dat files ----------
sec("SELF-CHECK  (.tex legend counts vs .dat line counts)")
tex = open(os.path.join(OUT, "fig_thm1iii_violation_native.tex"), encoding="utf-8").read()
expect = {"n3_pub_frac_lo.dat": ("$w_{\mathcal S}<0.99$", counts["pub_frac_lo"]),
          "n3_pub_frac_hi.dat": ("$w_{\mathcal S}\ge 0.99$", counts["pub_frac_hi"]),
          "n3_thmB_frac.dat": ("Theorem~\\ref{thm:leakage}, same $\\mathcal S$", counts["thmB_frac"]),
          "n3_inf_frac.dat": ("ratio $=\infty$", counts["inf_frac"]),
          "n3_fig5_frac.dat": ("Fig.~\\ref{fig:akwsampled} setting", counts["fig5_frac"])}
bad = 0
for dat, (key, n) in expect.items():
    real = sum(1 for _ in open(os.path.join(OUT, dat))) - 1
    idx = tex.find(key)
    seg = tex[idx:idx + 200] if idx >= 0 else ""
    inpar = re.search(r"\((\d+)\)", seg)
    shown = int(inpar.group(1)) if inpar else -1
    ok = (real == n == shown)
    bad += 0 if ok else 1
    print("  %-22s dat=%4d builder=%4d legend=%4d  %s" % (dat, real, n, shown, "OK" if ok else "MISMATCH"))
print("  ->", "ALL LEGEND COUNTS MATCH THE DATA" if bad == 0 else "%d MISMATCH(ES) -- FIX BEFORE SHIPPING" % bad)
assert bad == 0, "legend counts do not match the data"
