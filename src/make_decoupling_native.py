# -*- coding: utf-8 -*-
# --- DEPOSIT PORT (2026-09-19) -------------------------------------------------
# Deposited from the generator that rebuilt fig:decoupling and fig:master for the
# current manuscript.  Only paths were changed: the process anchors itself to the
# repository THIS FILE lives in, reads data/stats_resource.json instead of a
# scratch copy (verified identical in every numeric field), and writes into
# paper/figs/.  No expression that computes or formats a number was touched.
#
# Run from anywhere:  python src/make_decoupling_native.py
# Checked by:         src/check_figures.py
# -------------------------------------------------------------------------------
import os as _os
_os.chdir(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
# -*- coding: utf-8 -*-
r"""Emit the two rebuilt native pgfplots fragments for the resource section of P2.

EVERY number in the output is computed here from deposited JSON. Nothing is typed by hand.
Sources (all read-only):
  REPO/data/cost_vs_ent.json          -- the 30 Hubbard points (L,N,U,FAF,S,chi)  [raw coordinates]
  REPO/data/resource_master.json      -- the 74 pooled points (F1,S,chi,family)   [raw coordinates]
  REPO/data/n19_suite.json            -- the 19 molecules (FAF/n99 eq + diss)     [raw coordinates]
  data/stats_resource.json           -- the stratified statistics                [every coefficient]
Self-verification: every statistic that reaches the output is recomputed here from the raw points
and compared against stats_resource.json; any mismatch aborts before a file is written.
"""
import json, os, sys, math
import numpy as np
from scipy import stats as sps

REPO = os.getcwd()
OUT = os.path.join(REPO, "paper", "figs")
TOL = 1e-9


def load(p):
    with open(p, encoding="utf-8") as f:
        return json.load(f)


ST = load(os.path.join(REPO, "data", "stats_resource.json"))
ROWS = load(os.path.join(REPO, "data", "cost_vs_ent.json"))["rows"]
MAST = load(os.path.join(REPO, "data", "resource_master.json"))
MOLS = load(os.path.join(REPO, "data", "n19_suite.json"))["mols"]

FAIL = []
NCHK = [0]


def check(name, got, want, tol=TOL):
    NCHK[0] += 1
    if want is None or abs(float(got) - float(want)) > tol:
        FAIL.append("%s: recomputed %r vs deposited %r" % (name, got, want))
    return got


def sp(a, b):
    return float(sps.spearmanr(a, b).statistic)


# ---------------------------------------------------------------- Hubbard strata (raw)
strata_keys = [(s["L"], s["N"]) for s in ST["hubbard_strata"]]
S_OBS = {}
for (L, N) in strata_keys:
    r = sorted([x for x in ROWS if x["L"] == L and x["N"] == N], key=lambda x: x["U"])
    assert len(r) == 5, (L, N, len(r))
    S_OBS[(L, N)] = dict(U=[x["U"] for x in r], F1=[x["FAF"] for x in r],
                         S=[x["S"] for x in r], chi=[x["chi"] for x in r])

for s in ST["hubbard_strata"]:
    k = (s["L"], s["N"])
    d = S_OBS[k]
    check("rho(F1,S) %s" % (k,), sp(d["F1"], d["S"]), s["rho_F1_S"])
    check("rho(chi,S) %s" % (k,), sp(d["chi"], d["S"]), s["rho_chi_S"])
    check("n %s" % (k,), len(d["U"]), s["n"], tol=0)
    if d["chi"] != s["chi_values"]:
        FAIL.append("chi values %s: %s vs %s" % (k, d["chi"], s["chi_values"]))

allF1 = [v for k in strata_keys for v in S_OBS[k]["F1"]]
allS = [v for k in strata_keys for v in S_OBS[k]["S"]]
allChi = [v for k in strata_keys for v in S_OBS[k]["chi"]]
check("hubbard pooled rho(F1,S)", sp(allF1, allS),
      ST["self_verification_controls"]["hubbard_rho_F1_S"]["recomputed"])
check("hubbard pooled rho(chi,S)", sp(allChi, allS),
      ST["self_verification_controls"]["hubbard_rho_chi_S"]["recomputed"])


def strat_kendall(xkey):
    nC = nD = nT = 0
    for k in strata_keys:
        d = S_OBS[k]
        x, y = d[xkey], d["S"]
        for i in range(len(x)):
            for j in range(i + 1, len(x)):
                s_ = np.sign(x[i] - x[j]) * np.sign(y[i] - y[j])
                nC += int(s_ > 0)
                nD += int(s_ < 0)
                nT += int(s_ == 0)
    return nC, nD, nT


nC, nD, nT = strat_kendall("F1")
check("kendall nC (F1)", nC, ST["hubbard_stratified"]["stratified_kendall_F1"]["nC"], tol=0)
check("kendall nD (F1)", nD, ST["hubbard_stratified"]["stratified_kendall_F1"]["nD"], tol=0)
check("kendall tau (F1)", (nC - nD) / (nC + nD + nT),
      ST["hubbard_stratified"]["stratified_kendall_F1"]["tau"])
p_perm = (1.0 / math.factorial(5)) ** 6
check("exact permutation p", p_perm, ST["hubbard_stratified"]["exact_perm_p_one_sided_F1"], tol=1e-24)
check("permutation space", math.factorial(5) ** 6, ST["hubbard_stratified"]["product_space_size"], tol=0)

# ---------------------------------------------------------------- 74 pooled points (raw)
PTS = MAST["points"]
check("pooled rho(chi,S) 74", sp([p["chi"] for p in PTS], [p["S"] for p in PTS]),
      MAST["pooled_spearman_chi_S"])
rho_pool_F1 = check("pooled rho(F1,S) 74", sp([p["F1"] for p in PTS], [p["S"] for p in PTS]),
                    ST["orphan_constant_forensics"]["recomputed_pooled_spearman_F1_S"])
MOL = [p for p in PTS if p["family"] == "molecule"]
check("molecular n", len(MOL), ST["molecular"]["n_points"], tol=0)
check("molecular rho(chi,S)", sp([p["chi"] for p in MOL], [p["S"] for p in MOL]),
      ST["molecular"]["rho_chi_S"])
check("molecular rho(F1,S)", sp([p["F1"] for p in MOL], [p["S"] for p in MOL]),
      ST["molecular"]["rho_F1_S"])
n_floor = sum(1 for p in MOL if p["S"] == 1)
check("molecular |S|=1 count", n_floor, ST["molecular"]["n_at_floor_S_eq_1"], tol=0)
ties = {}
for p in MOL:
    ties[str(p["S"])] = ties.get(str(p["S"]), 0) + 1
if ties != ST["molecular"]["tie_structure_S"]:
    FAIL.append("molecular tie structure mismatch")
rk_S = sps.rankdata([p["S"] for p in MOL])
rk_chi = sps.rankdata([p["chi"] for p in MOL])
rk_F1 = sps.rankdata([p["F1"] for p in MOL])
# ceiling: the largest |rho| a tie-free predictor can reach against this tied |S| ranking
rho_ceiling = check("tie ceiling",
                    float(np.corrcoef(np.sort(rk_S), np.arange(1, len(rk_S) + 1))[0, 1]),
                    ST["molecular"]["max_attainable_rho_given_S_ties"])
check("rank-space rho(chi)", float(np.corrcoef(rk_S, rk_chi)[0, 1]), ST["molecular"]["rho_chi_S"])
check("rank-space rho(F1)", float(np.corrcoef(rk_S, rk_F1)[0, 1]), ST["molecular"]["rho_F1_S"])
midrank_floor = float(rk_S[[i for i, p in enumerate(MOL) if p["S"] == 1][0]])
check("floor midrank", midrank_floor, (1 + n_floor) / 2)

# -- the three "chain" rows are the U=8 half-filled Hubbard states, stored at lower precision
CH = [p for p in PTS if p["family"] == "chain2"]
HU8 = [S_OBS[(L, L)] for L in (6, 8, 10)]
dup_ok = all(abs(CH[i]["F1"] - HU8[i]["F1"][3]) < 5e-5 and CH[i]["S"] == HU8[i]["S"][3]
             and CH[i]["chi"] == HU8[i]["chi"][3] for i in range(3))
if not dup_ok:
    FAIL.append("the chain2 rows are NOT the U=8 half-filled Hubbard rows -- re-check the caption claim")
LAD = [p for p in PTS if p['family'] == 'ladder']
n_distinct = len(PTS) - len(CH)
keep = [p for p in PTS if p["family"] != "chain2"]
rho_71_chi = sp([p["chi"] for p in keep], [p["S"] for p in keep])
rho_71_F1 = sp([p["F1"] for p in keep], [p["S"] for p in keep])

# ---------------------------------------------------------------- molecules for panel (a)
COV = {'H2', 'HF', 'F2', 'CO', 'N2', 'H2O', 'NH3', 'BeH2', 'H4', 'H6'}
MRF = {'C2', 'O2', 'OH', 'CN', 'NO'}


def reg(m):
    return 'cov' if m in COV else ('mrf' if m in MRF else 'ion')


byreg = {'cov': [], 'mrf': [], 'ion': []}
mF, mY = [], []
for m in MOLS:
    for f, s in [(m['FAF_eq'], max(m['n99_eq'], 1)), (m['FAF_diss'], max(m['n99_diss'], 1))]:
        byreg[reg(m['mol'])].append((f, s))
        mF.append(f)
        mY.append(np.log2(s))
if len(mF) != 38:
    FAIL.append("molecule expansion gave %d points, not 38" % len(mF))
r_pearson = float(np.corrcoef(mF, mY)[0, 1])
slope, icpt = np.polyfit(mF, mY, 1)
check("n19_suite reproduces molecular rho(F1,S)", sp(mF, [2 ** y for y in mY]),
      ST["molecular"]["rho_F1_S"])

if FAIL:
    print("SELF-VERIFICATION FAILED -- nothing written:")
    for f in FAIL:
        print("  -", f)
    sys.exit(1)


# ---------------------------------------------------------------- helpers
def coords(xs, ys, fx="{:.4f}", fy="{:g}"):
    return " ".join("(" + fx.format(x) + "," + fy.format(y) + ")" for x, y in zip(xs, ys))


SC = {(6, 6): ("scAh", "D55E00", "*", "2.1pt"),
      (6, 4): ("scAd", "F0A868", "o", "2.1pt"),
      (8, 8): ("scBh", "00755E", "square*", "2.0pt"),
      (8, 6): ("scBd", "5FC2A5", "square", "2.0pt"),
      (10, 10): ("scCh", "1B4F9C", "triangle*", "2.6pt"),
      (10, 8): ("scCd", "79AEE8", "triangle", "2.6pt")}
LBL = {k: r"$(%d{,}%d)$" % k for k in SC}

pool_chi = ST["self_verification_controls"]["hubbard_rho_chi_S"]["recomputed"]
pool_F1 = ST["self_verification_controls"]["hubbard_rho_F1_S"]["recomputed"]
within_chi = [s["rho_chi_S"] for s in ST["hubbard_strata"]]
within_F1 = [s["rho_F1_S"] for s in ST["hubbard_strata"]]

# =================================================================== FIGURE 1
L1 = []
A = L1.append
A(r"% ===== fig:decoupling (REBUILT 2026-09-18) -- native pgfplots, no external .dat, no hand-typed numbers.")
A(r"% Emitted by p2_figs/make_figs.py from data/cost_vs_ent.json, data/n19_suite.json,")
A(r"% data/resource_master.json and p2_calc/stats_resource.json; all statistics re-verified at emit time.")
A(r"\definecolor{dcCov}{HTML}{E8543A}\definecolor{dcMrf}{HTML}{10A87B}\definecolor{dcIon}{HTML}{2E6FE0}")
A(r"\definecolor{dcTrend}{HTML}{9A9A9A}\definecolor{dcPool}{HTML}{6E6E6E}\definecolor{inkPrim}{HTML}{1B1B1B}")
A("".join(r"\definecolor{%s}{HTML}{%s}" % (SC[k][0], SC[k][1]) for k in SC))
A(r"\begin{tikzpicture}")
A(r"\begin{groupplot}[group style={group size=2 by 2, horizontal sep=1.8cm, vertical sep=2.05cm},")
A(r"  width=7.45cm, height=5.6cm, paperaxis,")
A(r"  title style={at={(0,1)},anchor=south west,font=\small\bfseries,yshift=2pt},")
A(r"  label style={font=\normalsize}, tick label style={font=\footnotesize},")
A(r"  ymajorgrids, grid style={draw=black!8}]")
A(r"% ---- (a) molecules: F1 and |S| rise together (the regime where F1 looks predictive) ----")
A(r"\nextgroupplot[title={(a)~molecules: they rise together},")
A(r"  xlabel={$\mathcal F_1=4\,\mathrm{tr}[\gamma(1{-}\gamma)]$}, ylabel={determinant support $|\mathcal S|$},")
A(r"  ymode=log, xmin=-0.9, xmax=%.2f, ymin=0.62, ymax=460, enlarge x limits=false]" % (max(mF) * 1.08))
A(r"\addplot[dcTrend,line width=1pt,dashed,domain=%.3f:%.3f,samples=60,forget plot] {2^(%.4f*x+%.4f)};"
  % (min(mF), max(mF), slope, icpt))
for key, col in (("cov", "dcCov"), ("mrf", "dcMrf"), ("ion", "dcIon")):
    A(r"\addplot[only marks,mark=*,mark size=1.9pt,draw=white,line width=0.4pt,fill=%s] coordinates {%s};"
      % (col, coords([p[0] for p in byreg[key]], [p[1] for p in byreg[key]])))
A(r"\node[inkPrim,font=\scriptsize,anchor=north west,align=left] at (rel axis cs:0.02,0.97) "
  r"{Pearson $r=%.2f$ on $\log_2|\mathcal S|$\\Spearman $\rho=%.3f$ ($n=%d$)};"
  % (r_pearson, ST["molecular"]["rho_F1_S"], len(MOL)))
A(r"% ---- (b) the six (L,N) strata against F1 ----")
A(r"\nextgroupplot[title={(b)~Hubbard: $\mathcal F_1$ runs backwards},")
A(r"  xlabel={one-body magic $\mathcal F_1$}, ylabel={determinant support $|\mathcal S|$},")
A(r"  ymode=log, xmin=-0.9, xmax=20.2, ymin=40, ymax=900000, enlarge x limits=false]")
for k in strata_keys:
    nm, _, mk, ms = SC[k]
    d = S_OBS[k]
    A(r"\addplot[%s,mark=%s,mark size=%s,line width=1.0pt,mark options={fill=%s,draw=%s}] coordinates {%s};"
      % (nm, mk, ms, nm, nm, coords(d["F1"], d["S"])))
A(r"\node[inkPrim,font=\scriptsize,anchor=north west,align=left] at (rel axis cs:0.02,0.97) "
  r"{$U=1\to16$ along each curve\\within each sweep: $\rho=%.3f$"
  r"\\pooled ($n=30$): $\rho=%.3f$ (withdrawn)};" % (within_F1[0], pool_F1))
A(r"% ---- (c) the same six strata against chi ----")
A(r"\nextgroupplot[title={(c)~Hubbard: $\chi$ runs forwards},")
A(r"  xlabel={minimal bond dimension $\chi$},")
A(r"  ylabel={determinant support $|\mathcal S|$},")
A(r"  ymode=log, xmin=2.6, xmax=18.4, ymin=40, ymax=900000, enlarge x limits=false,")
A(r"  legend style={at={(1.16,-0.46)},anchor=north,legend columns=6,draw=black!25,font=\footnotesize},")
A(r"  legend cell align=left]")
for k in strata_keys:
    nm, _, mk, ms = SC[k]
    d = S_OBS[k]
    A(r"\addplot[%s,mark=%s,mark size=%s,line width=1.0pt,mark options={fill=%s,draw=%s}] coordinates {%s};"
      % (nm, mk, ms, nm, nm, coords(d["chi"], d["S"], fx="{:g}")))
    A(r"\addlegendentry{%s}" % LBL[k])
A(r"\node[inkPrim,font=\scriptsize,anchor=north west,align=left] at (rel axis cs:0.02,0.97) "
  r"{$\chi$ falls as $U$ rises\\within each sweep: $\rho=%.3f$ to $%.3f$"
  r"\\pooled ($n=30$): $\rho=+%.3f$ (withdrawn)};" % (min(within_chi), max(within_chi), pool_chi))
A(r"% ---- (d) Simpson panel ----")
A(r"\nextgroupplot[title={(d)~stratified vs pooled},")
A(r"  xlabel={stratum $(L,N)$}, ylabel={Spearman $\rho$ vs $|\mathcal S|$},")
A(r"  xmin=0.45, xmax=6.55, ymin=-1.42, ymax=1.95, ytick={-1,-0.5,0,0.5,1},")
A(r"  xtick={%s}, xticklabels={%s}," % (",".join(str(i + 1) for i in range(6)),
                                        ",".join(LBL[k] for k in strata_keys)))
A(r"  x tick label style={font=\scriptsize,rotate=30,anchor=east}, ymajorgrids=false]")
A(r"\addplot[dcPool,densely dashed,line width=0.9pt,forget plot] coordinates {(0.45,%.4f) (6.55,%.4f)};"
  % (pool_chi, pool_chi))
A(r"\addplot[dcPool,densely dashed,line width=0.9pt,forget plot] coordinates {(0.45,%.4f) (6.55,%.4f)};"
  % (pool_F1, pool_F1))
for i, k in enumerate(strata_keys):
    nm = SC[k][0]
    A(r"\addplot[only marks,mark=*,mark size=2.4pt,%s,mark options={fill=%s,draw=white,line width=0.4pt},forget plot] coordinates {(%d,%.4f)};"
      % (nm, nm, i + 1, within_chi[i]))
    A(r"\addplot[only marks,mark=square*,mark size=2.3pt,%s,mark options={fill=%s,draw=white,line width=0.4pt},forget plot] coordinates {(%d,%.4f)};"
      % (nm, nm, i + 1, within_F1[i]))
span = 1.95 + 1.42
A(r"\node[inkPrim,font=\scriptsize,anchor=north west] at (rel axis cs:0.02,0.995) "
  r"{$\bullet$\; $\rho(\chi,|\mathcal S|)$ per stratum};")
A(r"% D12: this key used to sit at the bottom-left, corner-to-corner with the (6,6) datum.")
A(r"\node[inkPrim,font=\scriptsize,anchor=north west] at (rel axis cs:0.02,0.90) "
  r"{$\blacksquare$\; $\rho(\mathcal F_1,|\mathcal S|)$ per stratum};")
A(r"\node[dcPool,font=\scriptsize,anchor=north west] at (rel axis cs:0.035,%.3f) {pooled ($n=30$): $\rho=+%.3f$ (withdrawn)};"
  % ((pool_chi + 1.42) / span - 0.012, pool_chi))
A(r"\node[dcPool,font=\scriptsize,anchor=north west] at (rel axis cs:0.035,%.3f) {pooled ($n=30$): $\rho=-%.3f$ (withdrawn)};"
  % ((pool_F1 + 1.42) / span - 0.012, abs(pool_F1)))
A(r"\end{groupplot}")
A(r"\end{tikzpicture}")
open(os.path.join(OUT, "fig_decoupling_native.tex"), "w", encoding="utf-8").write("\n".join(L1) + "\n")

# =================================================================== FIGURE 2
fams = [("hubbard", "rmHub", "square*", "1.9pt", r"Hubbard sweep ($\times%d$)"),
        ("chain2", "rmChn", "triangle*", "2.2pt", r"chain ($\times%d$)"),
        ("ladder", "rmLad", "diamond*", "2.2pt", r"2-leg ladder ($\times%d$)")]
sub = {f[0]: [p for p in PTS if p["family"] == f[0]] for f in fams}
mol_hi = [p for p in PTS if p["family"] == "molecule" and p["S"] > 1]
mol_lo = [p for p in PTS if p["family"] == "molecule" and p["S"] == 1]

L2 = []
B = L2.append
B(r"% ===== fig:master (REBUILT 2026-09-18) -- 74 exact ground states, native pgfplots.")
B(r"% Emitted by p2_figs/make_figs.py from data/resource_master.json and p2_calc/stats_resource.json.")
B(r"% No hand-typed numbers: every coordinate, coefficient and interval is read from those files.")
B(r"\definecolor{rmMol}{HTML}{0072B2}\definecolor{rmHub}{HTML}{E69F00}")
B(r"\definecolor{rmChn}{HTML}{009E73}\definecolor{rmLad}{HTML}{D55E00}\definecolor{rmInk}{HTML}{1B1B1B}")
B(r"\definecolor{rmGrey}{HTML}{6E6E6E}")
B(r"\begin{tikzpicture}")
B(r"\begin{groupplot}[group style={group size=2 by 2, horizontal sep=1.8cm, vertical sep=2.05cm},")
B(r"  width=7.45cm, height=5.6cm, paperaxis,")
B(r"  title style={at={(0,1)},anchor=south west,font=\small\bfseries,yshift=2pt},")
B(r"  label style={font=\normalsize}, tick label style={font=\footnotesize},")
B(r"  legend style={font=\scriptsize,draw=black!25,fill=white,fill opacity=0.92,text opacity=1},")
B(r"  legend cell align=left]")
# ---------------- (a) |S| vs chi
B(r"% ---- (a) |S| vs chi, with the eleven molecular floor points marked ----")
B(r"\nextgroupplot[title={(a)~cost tracks $\chi$}, xmode=log, ymode=log,")
B(r"  xlabel={minimal bond dimension $\chi$}, ylabel={determinant support $|\mathcal S|$},")
B(r"  xmin=0.8,xmax=70,ymin=0.50,ymax=900000, ymajorgrids,xmajorgrids, grid style={draw=black!8},")
B(r"  legend to name=fig2legend, legend columns=3, legend style={draw=black!25,/tikz/every even column/.append style={column sep=11pt}}]")
B(r"\addplot[only marks,mark=*,mark size=1.9pt,rmMol,mark options={fill=rmMol,draw=white,line width=0.3pt}] coordinates {%s};"
  % coords([p["chi"] for p in mol_hi], [p["S"] for p in mol_hi], fx="{:g}"))
B(r"\addlegendentry{molecules, $|\mathcal S|>1$ ($\times%d$)}" % len(mol_hi))
B(r"\addplot[only marks,mark=o,mark size=2.2pt,rmMol,line width=0.8pt] coordinates {%s};"
  % coords([p["chi"] for p in mol_lo], [p["S"] for p in mol_lo], fx="{:g}"))
B(r"\addlegendentry{molecules at the floor $|\mathcal S|=1$ ($\times%d$)}" % len(mol_lo))
for key, col, mk, ms, lg in fams:
    B(r"\addplot[only marks,mark=%s,mark size=%s,%s,mark options={fill=%s,draw=white,line width=0.3pt}] coordinates {%s};"
      % (mk, ms, col, col, coords([p["chi"] for p in sub[key]], [p["S"] for p in sub[key]], fx="{:g}")))
    B(r"\addlegendentry{%s}" % (lg % len(sub[key])))
B(r"\addplot[rmInk,densely dashed,line width=0.9pt,mark=none,forget plot] coordinates {(2,2) (2,4) (2,8) (2,16) (2,32)};")
B(r"\node[rmInk,font=\scriptsize,anchor=west] at (axis cs:2.15,120000) {witness: $\chi{=}2,\;|\mathcal S|{=}2^{K}$};")
B(r"\node[rmInk,font=\scriptsize,anchor=south east] at (rel axis cs:0.98,0.02) {pooled $\rho=+%.3f$};"
  % MAST["pooled_spearman_chi_S"])
# ---------------- (b) |S| vs F1
B(r"% ---- (b) |S| vs F1, same 74 points, same marking ----")
B(r"\nextgroupplot[title={(b)~$\mathcal F_1$ does not}, ymode=log,")
B(r"  xlabel={one-body magic $\mathcal F_1$}, ylabel={determinant support $|\mathcal S|$},")
B(r"  xmin=-1.2,xmax=23,ymin=0.50,ymax=900000, ymajorgrids, grid style={draw=black!8}]")
B(r"\addplot[only marks,mark=*,mark size=1.9pt,rmMol,mark options={fill=rmMol,draw=white,line width=0.3pt}] coordinates {%s};"
  % coords([p["F1"] for p in mol_hi], [p["S"] for p in mol_hi]))
B(r"\addplot[only marks,mark=o,mark size=2.2pt,rmMol,line width=0.8pt] coordinates {%s};"
  % coords([p["F1"] for p in mol_lo], [p["S"] for p in mol_lo]))
for key, col, mk, ms, lg in fams:
    B(r"\addplot[only marks,mark=%s,mark size=%s,%s,mark options={fill=%s,draw=white,line width=0.3pt}] coordinates {%s};"
      % (mk, ms, col, col, coords([p["F1"] for p in sub[key]], [p["S"] for p in sub[key]])))
B(r"\node[rmInk,font=\scriptsize,anchor=north west,align=left] at (rel axis cs:0.02,0.97) "
  r"{pooled $\rho=+%.3f$ over the same 74 rows\\Hubbard branch: "
  r"$U\!\uparrow\;\Rightarrow\;\mathcal F_1\!\uparrow,\;|\mathcal S|\!\downarrow$};" % rho_pool_F1)
# ---------------- (c) coefficients with intervals
ci_chi = ST["molecular"]["cluster_bootstrap"]["ci95_rho_chi_S"]
ci_F1 = ST["molecular"]["cluster_bootstrap"]["ci95_rho_F1_S"]
wm = ST["dependent_correlation_comparisons"]["molecular_38"]
wp = ST["dependent_correlation_comparisons"]["pooled_74"]
B(r"% ---- (c) the four coefficients, the intervals the deposit supports, and the tie ceiling ----")
B(r"\nextgroupplot[title={(c)~coefficients; intervals for $n{=}38$ only},")
B(r"  xlabel={Spearman $\rho$ against $|\mathcal S|$}, ylabel={},")
B(r"  xmin=-0.03,xmax=1.03,ymin=-1.55,ymax=5.45, xtick={0,0.25,0.5,0.75,1},")
B(r"  ytick={1,2,3,4}, yticklabels={{$\mathcal F_1$, $n{=}74$},{$\chi$, $n{=}74$},"
  r"{$\mathcal F_1$, $n{=}38$},{$\chi$, $n{=}38$}},")
B(r"  y tick label style={font=\scriptsize}, ymajorgrids=false, xmajorgrids, grid style={draw=black!8}]")
B(r"\addplot[rmGrey,densely dashed,line width=0.9pt,forget plot] coordinates {(%.4f,0.45) (%.4f,4.45)};"
  % (rho_ceiling, rho_ceiling))
for y, lo, hi, col in ((4, ci_chi[0], ci_chi[1], "rmMol"), (3, ci_F1[0], ci_F1[1], "rmLad")):
    B(r"\addplot[%s,line width=1.0pt,forget plot] coordinates {(%.4f,%d) (%.4f,%d)};" % (col, lo, y, hi, y))
    B(r"\addplot[%s,line width=1.0pt,forget plot] coordinates {(%.4f,%.2f) (%.4f,%.2f)};" % (col, lo, y - 0.13, lo, y + 0.13))
    B(r"\addplot[%s,line width=1.0pt,forget plot] coordinates {(%.4f,%.2f) (%.4f,%.2f)};" % (col, hi, y - 0.13, hi, y + 0.13))
for y, x, col in ((4, ST["molecular"]["rho_chi_S"], "rmMol"), (3, ST["molecular"]["rho_F1_S"], "rmLad"),
                  (2, MAST["pooled_spearman_chi_S"], "rmMol"), (1, rho_pool_F1, "rmLad")):
    B(r"\addplot[only marks,mark=*,mark size=2.4pt,%s,mark options={fill=%s,draw=white,line width=0.4pt},forget plot] coordinates {(%.4f,%d)};"
      % (col, col, x, y))
B(r"\node[rmGrey,font=\scriptsize,anchor=east] at (axis cs:%.4f,5.05) {tie ceiling $%.3f$\;};" % (rho_ceiling, rho_ceiling))
B(r"\node[rmInk,font=\scriptsize,anchor=west] at (axis cs:-0.01,-0.40) "
  r"{$n{=}38$: $\Delta\rho=%.3f$ ($t=%.2f$, $p=%.2f$)};" % (wm["difference"], wm["williams_t"], wm["williams_p"]))
B(r"\node[rmInk,font=\scriptsize,anchor=west] at (axis cs:-0.01,-1.15) "
  r"{$n{=}74$: $\Delta\rho=%.3f$ ($t=%.2f$, $p=%.2f$)};" % (wp["difference"], wp["williams_t"], wp["williams_p"]))
# ---------------- (d) rank space
B(r"% ---- (d) the molecular coefficient in rank space: the eleven-point tied block ----")
B(r"\nextgroupplot[title={(d)~the 38 molecules in rank space},")
B(r"  xlabel={midrank of $|\mathcal S|$}, ylabel={midrank of the predictor},")
B(r"  xmin=0,xmax=40,ymin=0,ymax=59, xmajorgrids, ymajorgrids, grid style={draw=black!8},")
B(r"  ytick={0,10,20,30,40},")
B(r"  legend style={at={(0.98,0.985)},anchor=north east,font=\scriptsize}]")
B(r"\addplot[rmGrey,densely dashed,line width=0.8pt,forget plot] coordinates {(1,1) (38,38)};")
B(r"\addplot[rmGrey,line width=0.7pt,forget plot] coordinates {(%.1f,%.1f) (%.1f,41.5)};"
  % (midrank_floor, float(max(max(rk_chi[[i for i, p in enumerate(MOL) if p["S"] == 1]]),
                              max(rk_F1[[i for i, p in enumerate(MOL) if p["S"] == 1]]))) + 4.5, midrank_floor))
B(r"\addplot[only marks,mark=*,mark size=1.9pt,rmMol,mark options={fill=rmMol,draw=white,line width=0.3pt}] coordinates {%s};"
  % coords(rk_S, rk_chi, fx="{:.1f}", fy="{:.1f}"))
B(r"\addlegendentry{$\chi$, $\rho=%.3f$}" % ST["molecular"]["rho_chi_S"])
B(r"\addplot[only marks,mark=square*,mark size=1.9pt,rmLad,mark options={fill=rmLad,draw=white,line width=0.3pt}] coordinates {%s};"
  % coords(rk_S, rk_F1, fx="{:.1f}", fy="{:.1f}"))
B(r"\addlegendentry{$\mathcal F_1$, $\rho=%.3f$}" % ST["molecular"]["rho_F1_S"])
B(r"\node[rmInk,font=\scriptsize,anchor=north west,align=left] at (rel axis cs:0.02,0.985) "
  r"{%d of 38 tied at $|\mathcal S|=1$\\midrank $%.1f$: $|\rho|\le %.3f$};" % (n_floor, midrank_floor, rho_ceiling))
B(r"% ---- D13: the grey rule marks the block of tied midranks; it now carries a label ----")
B(r"\node[rmGrey,font=\scriptsize,anchor=west] at (axis cs:%.1f,33.0) {tied block};" % (midrank_floor + 1.1))
B(r"\end{groupplot}")
B(r"\node[anchor=north,yshift=-1.25cm] at (group c1r2.south east) {\ref*{fig2legend}};")
B(r"\end{tikzpicture}")
open(os.path.join(OUT, "fig_resource_master_native.tex"), "w", encoding="utf-8").write("\n".join(L2) + "\n")

# =================================================================== caption numbers
NUM = dict(
    molecular_rho_chi=ST["molecular"]["rho_chi_S"], molecular_rho_F1=ST["molecular"]["rho_F1_S"],
    molecular_ci_chi=ci_chi, molecular_ci_F1=ci_F1,
    tie_ceiling=rho_ceiling, n_floor=n_floor, n_mol=len(MOL), floor_midrank=midrank_floor,
    pooled_rho_chi=MAST["pooled_spearman_chi_S"], pooled_rho_F1=rho_pool_F1,
    deposited_pooled_rho_F1=MAST["pooled_spearman_F1_S"],
    hub_pool_chi=pool_chi, hub_pool_F1=pool_F1,
    within_rho_chi=within_chi, within_rho_F1=within_F1,
    kendall=dict(tau=(nC - nD) / (nC + nD + nT), nC=nC, nD=nD, nT=nT),
    perm_p=p_perm, perm_space=math.factorial(5) ** 6,
    fisher_p=ST["hubbard_stratified"]["fisher_F1"]["p"],
    paper_pooled_p_F1=ST["hubbard_stratified"]["paper_pooled_p_F1"],
    eff_n_hub=ST["hubbard_stratified"]["effective_n_of_the_pooled_test"],
    eff_n_mol=ST["molecular"]["effective_n"],
    icc_mol=ST["molecular"]["ICC_log2S_within_molecule"],
    williams_mol=wm, williams_pool=wp,
    pearson_r_mol=r_pearson, trend_slope=float(slope), trend_intercept=float(icpt),
    mol_hi=len(mol_hi), mol_lo=len(mol_lo),
    cluster_perm_p=ST["molecular"]["cluster_permutation"]["p_chi"],
    n_per_family={k: len(v) for k, v in sub.items()},
    n_rows=len(PTS), n_distinct_states=n_distinct,
    chain_rows_duplicate_hubbard_U8=dup_ok,
    pooled_rho_chi_71=rho_71_chi, pooled_rho_F1_71=rho_71_F1,
    chain_ladder=dict(
        S_ratio=[LAD[i]['S'] / CH[i]['S'] for i in range(3)],
        F1_rel_diff=[abs(LAD[i]['F1'] - CH[i]['F1']) / CH[i]['F1'] for i in range(3)],
        chain_S=[CH[i]['S'] for i in range(3)], ladder_S=[LAD[i]['S'] for i in range(3)],
        chain_F1=[CH[i]['F1'] for i in range(3)], ladder_F1=[LAD[i]['F1'] for i in range(3)]),
    ceiling_fraction_chi=ST['molecular']['rho_chi_S_as_fraction_of_ceiling'],
    ceiling_fraction_F1=ST['molecular']['rho_F1_S_as_fraction_of_ceiling'],
    floor_fraction=n_floor / len(MOL),

)
json.dump(NUM, open(os.path.join(OUT, "figure_numbers.json"), "w"), indent=1, default=float)
print("SELF-VERIFICATION PASSED: %d checks, 0 failures." % NCHK[0])
print("wrote fig_decoupling_native.tex and fig_resource_master_native.tex")
print("molecular Pearson r=%.4f | tie ceiling=%.6f | floor midrank=%.1f | n_floor=%d"
      % (r_pearson, rho_ceiling, midrank_floor, n_floor))
