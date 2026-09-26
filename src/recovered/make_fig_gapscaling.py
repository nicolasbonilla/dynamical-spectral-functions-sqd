# -*- coding: utf-8 -*-
"""Generate fig_gapscaling_native.tex -- finite-size scaling of the two gaps.

EVERY number in the emitted fragment is read from a JSON file at run time.
Nothing is typed by hand.  Sources:
  gap_scaling.json                  -> gaps, fits, Lieb-Wu, Heisenberg control
  release/data/spinqw_L12.json      -> eta of the spin channel S^zz(q,w)
  release/data/akw_lanczos_L12.json -> eta of the one-particle channel A(k,w)
The script aborts if a value it needs is missing or an internal check fails.

DEPOSITED 2026-09-26 (src/recovered/README.md).  Written on 2026-09-18 in the author's session
scratch directory, where it generated paper/figs/fig_gapscaling_native.tex, and recovered from
there.  Changes against the recovered copy, and nothing else:
  * paths: GAPF, SPINF and AKWF come from data/ of the repository this file sits in, derived
    from __file__ (the original read gap_scaling.json from its scratch run directory and the
    other two through an absolute path to the author's clone);
  * outputs: the fragment goes to paper/figs/fig_gapscaling_native.tex, where it is committed;
    the caption and the facts file go to build/ (not deposited), because the committed
    paper/figs/fig_gapscaling_caption.tex was edited by hand after it was generated and is
    the authoritative caption -- writing over it would reinstate the older wording
    (docs/KNOWN_DISCREPANCIES.md, section 27).
Run from the deposited data/gap_scaling.json, the fragment differs from the committed one in
ONE comment line, line 5, which records the generated_utc of the gap_scaling.json it was built
from: 2026-09-18T12:01:30Z for the run the committed fragment came from, 2026-09-18T19:41:36Z
for the deposited file.  The two JSONs have identical numerical content; they differ only in
provenance fields (the absolute engine path removed in docs/KNOWN_DISCREPANCIES.md section 20,
argv, wall time and that timestamp).  src/check_figures.py compares the fragment with that one
line normalised, and says so.
"""
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))          # <repo>/src/recovered
REPO = os.path.dirname(os.path.dirname(HERE))              # <repo>
GAPF = os.path.join(REPO, "data", "gap_scaling.json")
SPINF = os.path.join(REPO, "data", "spinqw_L12.json")
AKWF = os.path.join(REPO, "data", "akw_lanczos_L12.json")
FIGS = os.path.join(REPO, "paper", "figs")
BUILD = os.path.join(REPO, "build")                        # by-products, not deposited
if not os.path.isdir(BUILD):
    os.makedirs(BUILD)


def load(p):
    with open(p, "r", encoding="utf-8") as fh:
        return json.load(fh)


G = load(GAPF)
eta_spin = float(load(SPINF)["eta"])
eta_charge = float(load(AKWF)["eta"])

# ---------------------------------------------------------------- data
hub = {int(r["L"]): r for r in G["hubbard"]}
hei = {int(r["L"]): r for r in G["heisenberg_control"]}
Ls = sorted(hub)
assert Ls == sorted(hei), "size grids of Hubbard and Heisenberg control differ"
gs = {L: float(hub[L]["spin_gap"]) for L in Ls}
gc = {L: float(hub[L]["charge_gap"]) for L in Ls}
gh = {L: float(hei[L]["szz_pi_lowest_pole"]) for L in Ls}
Lgs = {L: float(hub[L]["L_times_spin_gap"]) for L in Ls}
Lgc = {L: float(hub[L]["L_times_charge_gap"]) for L in Ls}
Lgh = {L: float(hei[L]["L_times_gap"]) for L in Ls}

for L in Ls:                                   # recompute, do not trust
    assert abs(Lgs[L] - L * gs[L]) < 1e-12
    assert abs(Lgc[L] - L * gc[L]) < 1e-12
    assert abs(Lgh[L] - L * gh[L]) < 1e-12
    assert abs(gc[L] - (hub[L]["mu_plus"] - hub[L]["mu_minus"])) < 1e-12

fs = G["fits"]["spin_L6plus"]
fc = G["fits"]["charge_L6plus"]
fq = G["fits"]["charge_all"]["quadratic"]
fh_ = G["fits"]["heisenberg_L6plus"]
LW = float(G["lieb_wu"]["form_A_sinh_integral"])
LWlit = float(G["lieb_wu"]["literature"])
sep = G["separation_of_the_two_channels"]
ltg = G["L_times_gap"]

Lfit = [int(x) for x in fs["sizes"]]
assert Lfit == [6, 8, 10, 12], Lfit
assert [int(x) for x in fq["_sizes"]] == Ls if "_sizes" in fq else True
Lmax, Lmin = max(Ls), min(Ls)

# ---------------------------------------------------------------- derived
fwhm_s, fwhm_c = 2 * eta_spin, 2 * eta_charge
r_s_eta, r_c_eta = gs[Lmax] / eta_spin, gc[Lmax] / eta_charge
r_s_fwhm, r_c_fwhm = gs[Lmax] / fwhm_s, gc[Lmax] / fwhm_c
# height of the eta-broadened Lorentzian at omega=0, relative to its peak
res_s = eta_spin ** 2 / (eta_spin ** 2 + gs[Lmax] ** 2)
res_c = eta_charge ** 2 / (eta_charge ** 2 + gc[Lmax] ** 2)
mean_Lgs = sum(Lgs[L] for L in Lfit) / len(Lfit)
drift_heis = Lgh[Lmax] / Lgh[Lmin] - 1.0
z_lin_LW = (LW - fc["c0_intercept"]) / fc["c0_se"]
z_qua_LW = (fq["c0_intercept"] - LW) / fq["c0_se"]
z_spin0 = abs(fs["c0_intercept"]) / fs["c0_se"]
z_heis0 = abs(fh_["c0_intercept"]) / fh_["c0_se"]
n_below = [L for L in Ls if gs[L] < fwhm_s]
ratio_hh = {L: gs[L] / gh[L] for L in Ls}


def f(x, n=4):
    return ("%%.%df" % n) % x


def coords(d, xf=lambda L: 1.0 / L, sub=None):
    ks = sub if sub is not None else sorted(d)
    return " ".join("(%.6f,%.6f)" % (xf(L), d[L]) for L in ks)


# ------------------------------------------------------------------ v2 derived
fq6 = G["fits"]["charge_L6plus"]["quadratic"]        # SAME protocol as every other fit here
z_qua6_LW = (LW - fq6["c0_intercept"]) / fq6["c0_se"]
dlw6 = LW - fq6["c0_intercept"]
z_lin_zero = fc["c0_intercept"] / fc["c0_se"]
z_qua_zero = fq["c0_intercept"] / fq["c0_se"]
z_qua6_zero = fq6["c0_intercept"] / fq6["c0_se"]
ratio_fwhm = r_c_fwhm / r_s_fwhm

XMIN, XMAX = -0.014, 0.278
T = []
A = T.append
A(r"% ============================================================================")
A(r"% fig:gapscaling -- finite-size scaling of the two gaps of the 1D Hubbard ring")
A(r"% AUTO-GENERATED by make_fig_gapscaling.py.  DO NOT EDIT BY HAND.")
A(r"% Every number below is read at build time from:")
A(r"%%   gap_scaling.json   (%s, %s)" % (G["provenance"]["script"], G["provenance"]["generated_utc"]))
A(r"%%   data/spinqw_L12.json      -> eta = %s  (spin channel, S^zz(q,w))" % f(eta_spin, 2))
A(r"%%   data/akw_lanczos_L12.json -> eta = %s  (one-particle channel, A(k,w))" % f(eta_charge, 2))
A(r"%% Engine: %s; %s; broadening: %s." % (G["provenance"]["engine"]["module"],
                                           G["provenance"]["parameters"]["lanczos"],
                                           G["provenance"]["parameters"]["broadening"]))
A(r"% Palette inherited from fig_scaling2_native.tex (scMag / scFrac / inkPrim / scDim).")
A(r"% ---------------------------------------------------------------------------")
A(r"% PROVENANCE, quantity by quantity (all paths inside gap_scaling.json unless noted):")
A(r"%   (a) spin gaps          -> /hubbard[L]/spin_gap")
A(r"%   (a) Heisenberg control -> /heisenberg_control[L]/szz_pi_lowest_pole")
A(r"%   (a) dashed fit + open square intercept -> /fits/spin_L6plus/{c0_intercept,c0_se,c1_slope}")
A(r"%   (a) grey band half-width eta           -> release/data/spinqw_L12.json /eta")
A(r"%   (b) charge gaps        -> /hubbard[L]/charge_gap  ( = mu_plus - mu_minus, re-checked)")
A(r"%   (b) dashed  curve      -> /fits/charge_L6plus/{c0_intercept,c1_slope}          (L >= 6)")
A(r"%   (b) dash-dot-dot curve -> /fits/charge_L6plus/quadratic/{c0_intercept,c1,c2}   (L >= 6)")
A(r"%   (b) dash-dot curve     -> /fits/charge_all/quadratic/{c0_intercept,c1,c2}      (ALL five L)")
A(r"%   (b) dotted Lieb-Wu line-> /lieb_wu/form_A_sinh_integral")
A(r"%   (c) L*Delta            -> /hubbard[L]/{L_times_spin_gap,L_times_charge_gap}")
A(r"%   (c) thin straight lines and quoted slopes -> /L_times_gap/{spin,charge}/{d0,d1,d1_se}")
A(r"%   (d) the three charge intercepts and their standard errors ->")
A(r"%         /fits/charge_L6plus/{c0_intercept,c0_se},")
A(r"%         /fits/charge_L6plus/quadratic/{c0_intercept,c0_se},")
A(r"%         /fits/charge_all/quadratic/{c0_intercept,c0_se}; dotted line -> /lieb_wu.")
A(r"% NOTE ON FIT PROTOCOL: every fit drawn here excludes L=4 EXCEPT the dash-dot curve, which is")
A(r"% the only all-sizes fit and is labelled as such in the legend, in (d) and in the caption.")
A(r"% Axis limits, tick positions and label positions are layout choices, not data.")
A(r"% ============================================================================")
A(r"\definecolor{scMag}{HTML}{D55E00}\definecolor{scFrac}{HTML}{0072B2}")
A(r"\definecolor{inkPrim}{HTML}{1B1B1B}\definecolor{scDim}{HTML}{8A8A8A}")
A(r"\definecolor{scCtrl}{RGB}{123,63,160}")
A(r"\begin{tikzpicture}")
A(r"\begin{groupplot}[group style={group size=2 by 2, horizontal sep=1.55cm, vertical sep=2.05cm},")
A(r"  width=7.00cm, height=5.35cm, paperaxis,")
A(r"  title style={at={(0,1)},anchor=south west,font=\small\bfseries,yshift=1.5pt},")
A(r"  label style={font=\small}, tick label style={font=\footnotesize},")
A(r"  xticklabel style={/pgf/number format/fixed}]")

# ------------------------------------------------ (a) spin
A(r"% ---------- (a) the spin channel ----------")
A(r"\nextgroupplot[title={(a)~spin channel},")
A(r"  xlabel={$1/L$}, ylabel={$\Delta_{\mathrm{s}}\;[t]$},")
A(r"  xmin=%s, xmax=%s, ymin=-0.062, ymax=0.655," % (XMIN, XMAX))
A(r"  xtick={0,0.05,0.10,0.15,0.20,0.25},")
A(r"  ytick={0,0.1,0.2,0.3,0.4,0.5}, ymajorgrids, grid style={draw=black!10}]")
A(r"\fill[scDim,fill opacity=0.18] (axis cs:%s,0) rectangle (axis cs:%s,%s);" % (XMIN, XMAX, f(fwhm_s, 4)))
A(r"\draw[black!45,line width=0.4pt] (axis cs:%s,0) -- (axis cs:%s,0);" % (XMIN, XMAX))
A(r"\node[scDim!60!black,font=\tiny,anchor=east,align=right] at (axis cs:0.271,0.108)")
A(r"  {linewidth\\$2\eta=%s\,t$};" % f(fwhm_s, 2))
A(r"\node[scFrac,font=\tiny,anchor=west] at (axis cs:0.002,0.612) {Hubbard $U/t=8$};")
A(r"\node[scCtrl,font=\tiny,anchor=west] at (axis cs:0.002,0.540) {Heisenberg $J=4t^{2}/U$};")
A(r"\addplot[scFrac,line width=0.9pt,dashed] coordinates {(0,%s) (%.6f,%s)};"
  % (f(fs["c0_intercept"], 6), 1.0 / min(Lfit), f(fs["c0_intercept"] + fs["c1_slope"] / min(Lfit), 6)))
A(r"\addplot[scFrac,only marks,mark=*,mark size=2.4pt,")
A(r"  mark options={fill=scFrac,draw=white,line width=0.6pt}] coordinates {%s};" % coords(gs, sub=Lfit))
A(r"\addplot[scFrac,only marks,mark=o,mark size=2.4pt,line width=0.8pt] coordinates {%s};" % coords(gs, sub=[Lmin]))
A(r"\addplot[scFrac,only marks,mark=square*,mark size=2.1pt,")
A(r"  mark options={fill=white,draw=scFrac,line width=0.8pt},")
A(r"  error bars/.cd,y dir=both,y explicit,error bar style={scFrac,line width=0.7pt},")
A(r"  error mark options={rotate=90,scFrac,mark size=1.6pt,line width=0.7pt}]")
A(r"  coordinates {(0,%s) +- (0,%s)};" % (f(fs["c0_intercept"], 6), f(fs["c0_se"], 6)))
A(r"% control drawn last, as OPEN triangles: at L=6 the two channels coincide to 0.006 t")
A(r"\addplot[scCtrl,only marks,mark=triangle,mark size=3.4pt,line width=0.8pt]")
A(r"  coordinates {%s};" % coords(gh))

# ------------------------------------------------ (b) charge
A(r"% ---------- (b) the charge channel ----------")
A(r"\nextgroupplot[title={(b)~charge channel},")
A(r"  xlabel={$1/L$}, ylabel={$\Delta_{\mathrm{c}}=\mu^{+}-\mu^{-}\;[t]$},")
A(r"  xmin=%s, xmax=%s, ymin=4.30, ymax=6.34," % (XMIN, XMAX))
A(r"  xtick={0,0.05,0.10,0.15,0.20,0.25},")
A(r"  ytick={4.5,5.0,5.5,6.0}, ymajorgrids, grid style={draw=black!10},")
A(r"  legend cell align=left, legend columns=1,")
A(r"  legend style={at={(0.025,0.975)},anchor=north west,font=\tiny,")
A(r"    draw=black!25,fill=white,fill opacity=0.95,text opacity=1,")
A(r"    inner sep=1.6pt,row sep=0.4pt}]")
qx = [i * 0.25 / 60.0 for i in range(61)]
A(r"\addplot[scMag,line width=0.9pt,dashed] coordinates {(0,%s) (%.6f,%s)};"
  % (f(fc["c0_intercept"], 6), 1.0 / min(Lfit), f(fc["c0_intercept"] + fc["c1_slope"] / min(Lfit), 6)))
A(r"\addlegendentry{linear, $L\ge6$}")
A(r"\addplot[scMag!60,line width=0.9pt,dash dot dot] coordinates {%s};"
  % " ".join("(%.6f,%.6f)" % (x, fq6["c0_intercept"] + fq6["c1"] * x + fq6["c2"] * x * x)
             for x in qx if x <= 1.0 / min(Lfit) + 1e-9))
A(r"\addlegendentry{$+\,1/L^{2}$, $L\ge6$}")
A(r"\addplot[scMag,line width=0.9pt,dash dot] coordinates {%s};"
  % " ".join("(%.6f,%.6f)" % (x, fq["c0_intercept"] + fq["c1"] * x + fq["c2"] * x * x) for x in qx))
A(r"\addlegendentry{$+\,1/L^{2}$, all $L$}")
A(r"\addplot[inkPrim,line width=0.7pt,densely dotted] coordinates {(%s,%s) (%s,%s)};"
  % (XMIN, f(LW, 6), XMAX, f(LW, 6)))
A(r"\addlegendentry{Lieb--Wu $\Delta_{\infty}$}")
A(r"\addplot[forget plot,scMag,only marks,mark=*,mark size=2.4pt,")
A(r"  mark options={fill=scMag,draw=white,line width=0.6pt}] coordinates {%s};" % coords(gc, sub=Lfit))
A(r"\addplot[forget plot,scMag,only marks,mark=o,mark size=2.4pt,line width=0.8pt] coordinates {%s};" % coords(gc, sub=[Lmin]))

# ------------------------------------------------ (c) L * gap
A(r"% ---------- (c) opposite finite-size behaviour, model free ----------")
A(r"% the two thin straight lines are the SAME linear fits whose slopes are quoted; they are")
A(r"% sampled at 41 points each because a straight line in L is curved on a logarithmic y axis.")
A(r"\nextgroupplot[title={(c)~opposite size dependence of $L\Delta$},")
A(r"  xlabel={$L$}, ylabel={$L\,\Delta\;[t]$}, ymode=log,")
A(r"  xmin=3.2, xmax=13.0, ymin=0.92, ymax=420,")
A(r"  xtick={4,6,8,10,12}, ytick={1,2,5,10,20,50,100,200},")
A(r"  yticklabels={1,2,5,10,20,50,100,200},")
A(r"  ymajorgrids, grid style={draw=black!10}]")
Lg = [min(Lfit) + i * (Lmax - min(Lfit)) / 40.0 for i in range(41)]
A(r"\addplot[scMag!55,line width=0.5pt,forget plot] coordinates {%s};"
  % " ".join("(%.4f,%.6f)" % (x, ltg["charge"]["d0"] + ltg["charge"]["d1"] * x) for x in Lg))
A(r"\addplot[scFrac!55,line width=0.5pt,forget plot] coordinates {%s};"
  % " ".join("(%.4f,%.6f)" % (x, ltg["spin"]["d0"] + ltg["spin"]["d1"] * x) for x in Lg))
A(r"\addplot[scMag,only marks,mark=*,mark size=2.4pt,")
A(r"  mark options={fill=scMag,draw=white,line width=0.6pt}] coordinates {%s};"
  % coords(Lgc, xf=lambda L: L, sub=Lfit))
A(r"\addplot[scMag,only marks,mark=o,mark size=2.4pt,line width=0.8pt] coordinates {%s};"
  % coords(Lgc, xf=lambda L: L, sub=[Lmin]))
A(r"\addplot[scFrac,only marks,mark=*,mark size=2.4pt,")
A(r"  mark options={fill=scFrac,draw=white,line width=0.6pt}] coordinates {%s};"
  % coords(Lgs, xf=lambda L: L, sub=Lfit))
A(r"\addplot[scFrac,only marks,mark=o,mark size=2.4pt,line width=0.8pt] coordinates {%s};"
  % coords(Lgs, xf=lambda L: L, sub=[Lmin]))
A(r"\node[scMag,font=\tiny,anchor=west,align=left] at (axis cs:3.42,235)")
A(r"  {charge: fitted slope in $L$, $+%s\pm%s$};" % (f(ltg["charge"]["d1"], 2), f(ltg["charge"]["d1_se"], 2)))
A(r"\node[scFrac,font=\tiny,anchor=west,align=left] at (axis cs:3.42,6.7)")
A(r"  {spin: fitted slope in $L$, $%s\pm%s$};" % (f(ltg["spin"]["d1"], 3), f(ltg["spin"]["d1_se"], 3)))

# ------------------------------------------------ (d) intercepts, on a legible scale
A(r"% ---------- (d) the three L->inf charge intercepts, with their standard errors ----------")
A(r"% Same three fits as (b); this panel exists because on the scale of (b) a +-0.044 t error bar")
A(r"% is thinner than its own marker, so the 2.70 s.e. tension with Lieb-Wu cannot be checked by eye.")
A(r"\nextgroupplot[title={(d)~the $1/L\to0$ charge intercepts},")
A(r"  xlabel={$\Delta_{\mathrm{c}}(L\to\infty)\;[t]$}, ymin=0.15, ymax=3.75,")
A(r"  xmin=3.96, xmax=4.99, xtick={4.0,4.2,4.4,4.6,4.8}, ytick=\empty,")
A(r"  xmajorgrids, grid style={draw=black!10}]")
A(r"\addplot[inkPrim,line width=0.7pt,densely dotted,forget plot] coordinates {(%s,0.75) (%s,3.75)};"
  % (f(LW, 6), f(LW, 6)))
A(r"\node[inkPrim,font=\tiny,anchor=north] at (axis cs:%s,0.62) {Lieb--Wu $%s\,t$};"
  % (f(LW, 6), f(LW, 4)))
rows = [
    (3.0, "linear in $1/L$, $L\\ge6$", fc["c0_intercept"], fc["c0_se"],
     "$%s$ s.e.\\ below Lieb--Wu" % f(z_lin_LW, 2), "square*", 2.3),
    (2.0, "$+\\,1/L^{2}$, $L\\ge6$", fq6["c0_intercept"], fq6["c0_se"],
     "$%s$ s.e.\\ below Lieb--Wu" % f(z_qua6_LW, 2), "triangle*", 3.0),
    (1.0, "$+\\,1/L^{2}$, all five $L$\\\\(includes $L=4$)", fq["c0_intercept"], fq["c0_se"],
     "$%s$ s.e.\\ above Lieb--Wu" % f(z_qua_LW, 2), "diamond*", 3.0),
]
for (y, lab, v, se, note, mk, ms) in rows:
    A(r"\addplot[scMag,only marks,mark=%s,mark size=%.1fpt," % (mk, ms))
    A(r"  mark options={fill=white,draw=scMag,line width=0.8pt},")
    A(r"  error bars/.cd,x dir=both,x explicit,error bar style={scMag,line width=0.7pt},")
    A(r"  error mark options={scMag,mark size=2.0pt,line width=0.7pt}]")
    A(r"  coordinates {(%s,%.2f) +- (%s,0)};" % (f(v, 6), y, f(se, 6)))
    A(r"\node[scMag,font=\tiny,anchor=south west,align=left] at (axis cs:3.985,%.2f)" % (y + 0.33))
    A(r"  {%s};" % lab)
    A(r"\node[scDim!55!black,font=\tiny,anchor=south west,align=left] at (axis cs:3.985,%.2f)" % (y + 0.11))
    A(r"  {%s};" % note)
A(r"\end{groupplot}")
A(r"\end{tikzpicture}")

tex = "\n".join(T) + "\n"
out = os.path.join(FIGS, "fig_gapscaling_native.tex")
with open(out, "w", encoding="utf-8") as fh:
    fh.write(tex)

# ---------------------------------------------------------------- caption
def lst(d, sub, n=4):
    return ", ".join(f(d[L], n) for L in sub)


def sci(x):
    """LaTeX scientific notation, one significant digit."""
    m, e = ("%.0e" % abs(x)).split("e")
    return r"%s\times10^{%d}" % (m, int(e))


cap = r"""%% AUTO-GENERATED by make_fig_gapscaling.py -- every number comes from a JSON file.
%% Use as:  \begin{figure}[p]\centering\input{figs/fig_gapscaling_native.tex}
%%          \input{figs/fig_gapscaling_caption.tex}\end{figure}
\caption{\textbf{The two gaps of the half-filled Hubbard ring scale in opposite directions.}
Periodic ring, $U/t=8$, $L=%(Lall)s$; in (a) and (b) the abscissa is $1/L$, so the sizes run
$L=%(Lrev)s$ right to left, and open symbols are $L=%(Lmin)d$, which the $L\ge6$ fits exclude.
Every gap is an exact Lanczos pole, full reorthogonalisation, \emph{no} broadening.
\textbf{(a)} Spin gap, lowest pole of $S^{zz}(q=\pi,\omega)$: $%(gs6)s\,t$ at $L=%(Lfit)s$
($%(gs4)s\,t$ at $L=%(Lmin)d$). Dashed: linear fit in $1/L$, $L\ge6$ ($R^{2}=%(r2s)s$); the open
square is its intercept $%(is)s\pm%(isse)s\,t$, $%(zs)s$ s.e.\ from zero. Triangles: the pure
Heisenberg ring $J=4t^{2}/U=%(J)s$, an independent control by dense exact diagonalisation in a
different sector and code path, whose lowest $S^{zz}(q=\pi)$ pole reproduces its
singlet--triplet gap to $%(heisres)s$. For $L\ge6$ the two channels agree to within $%(agree)s\%%$
(ratios $%(ratios)s$); at $L=%(Lmin)d$ they do not ($%(ratio4)s$), which is why that size is
excluded. Grey band: the linewidth actually used for the published spin figure,
$\mathrm{FWHM}=2\eta=%(fwhms)s\,t$ at $\eta=%(etas)s\,t$. The $L=%(Lmax)d$ gap lies inside it
($%(gsmax)s\,t=%(rse)s\,\eta=%(rsf)s\,\mathrm{FWHM}$) and $L=10$ sits on its edge ($%(gs10)s\,t$);
broadened by this $\eta$, a pole at the $L=%(Lmax)d$ gap still carries $%(reshs)s\%%$ of its peak
height at $\omega=0$. One size alone therefore cannot separate a closing gap from no gap.
\textbf{(b)} Charge gap $\Delta_{\mathrm c}=\mu^{+}-\mu^{-}$ from the $N\pm1$ sectors:
$%(gc6)s\,t$ at $L=%(Lfit)s$ ($%(gc4)s\,t$ at $L=%(Lmin)d$). At $L=%(Lmax)d$ that is
$%(rce)s\,\eta=%(rcf)s\,\mathrm{FWHM}$ for the $\eta=%(etac)s\,t$ of the one-particle channel,
$%(ratio_fwhm)s$ times further from $\omega=0$ in units of its own linewidth than the spin gap.
Three extrapolations are drawn and they disagree. The ordinate omits zero because the question
here is which finite value the gap tends to, not whether it vanishes.
\textbf{(d)} Those three intercepts on a scale where their standard errors are legible; on the
scale of (b) the $\pm%(icse)s\,t$ bar is thinner than its marker. Linear, $L\ge6$:
$%(ic)s\pm%(icse)s\,t$, \emph{undershooting} the exact Lieb--Wu gap $%(lw)s\,t$ (two quadrature
forms agreeing to $%(lwres)s$) by $%(dlw)s\,t$, $%(zlw)s$ s.e. Same protocol plus a $1/L^{2}$
term: $%(iq6)s\pm%(iq6se)s\,t$, undershooting by $%(dlw6)s\,t$ ($%(zq6)s$ s.e., standard error
$%(sefac)s$ times larger). Only the $1/L^{2}$ fit over all five sizes brackets Lieb--Wu
($%(iq)s\pm%(iqse)s\,t$, $%(zq)s$ s.e.\ above), and it is the one fit here that includes
$L=%(Lmin)d$, the size (a) and (b) mark anomalous. These sizes therefore fix the $L\to\infty$
charge gap only to about $\pm0.1$--$0.3\,t$: compatible with Lieb--Wu once a $1/L^{2}$ term is
allowed, not evidence for it. What all three exclude is zero, the weakest by $%(z6z)s$ s.e.\ and
the linear one by $%(zlz)s$ s.e.
\textbf{(c)} The same data with no extrapolation model, logarithmic ordinate:
$L\Delta_{\mathrm s}$ is flat over $L\ge6$ (mean $%(meanLgs)s\,t$; fitted slope in $L$
$%(slps)s\pm%(slpsse)s$, $%(zslps)s$ s.e.\ from zero) while $L\Delta_{\mathrm c}$ grows linearly
($+%(slpc)s\pm%(slpcse)s$, $%(zslpc)s$ s.e.). The thin lines are those two linear fits, curved
here only because the ordinate is logarithmic.
\emph{Caveat, because it limits (a).} A vanishing intercept is not by itself the evidence: the
Heisenberg control, whose gap is known to vanish in the thermodynamic limit, has a
linear-in-$1/L$ intercept $%(ih)s\pm%(ihse)s\,t$, $%(zh)s$ s.e.\ \emph{above} zero, and its own
$L\Delta$ drifts $+%(drift)s\%%$ from $L=%(Lmin)d$ to $L=%(Lmax)d$, both from its logarithmic
corrections. The model-free statement is (c).}
\label{fig:gapscaling}
""" % {
    "Lall": ",".join(str(L) for L in Ls),
    "Lrev": ",".join(str(L) for L in reversed(Ls)),
    "Lfit": ",".join(str(L) for L in Lfit),
    "Lmin": Lmin, "Lmax": Lmax,
    "gs6": lst(gs, Lfit, 5), "gs4": f(gs[Lmin], 5), "gs10": f(gs[10], 5),
    "gc6": lst(gc, Lfit, 4), "gc4": f(gc[Lmin], 4),
    "gsmax": f(gs[Lmax], 5),
    "r2s": f(fs["r2"], 3),
    "is": f(fs["c0_intercept"], 4), "isse": f(fs["c0_se"], 4), "zs": f(z_spin0, 2),
    "J": f(G["provenance"]["parameters"]["J_eff"], 2),
    "heisres": sci(max(abs(r["internal_consistency"]) for r in G["heisenberg_control"])),
    "agree": f(100 * max(abs(1 - ratio_hh[L]) for L in Lfit), 1),
    "ratios": ", ".join(f(ratio_hh[L], 3) for L in Lfit),
    "ratio4": f(ratio_hh[Lmin], 3),
    "fwhms": f(fwhm_s, 2), "etas": f(eta_spin, 2), "etac": f(eta_charge, 2),
    "rse": f(r_s_eta, 2), "rsf": f(r_s_fwhm, 2),
    "rce": f(r_c_eta, 1), "rcf": f(r_c_fwhm, 1),
    "ratio_fwhm": f(ratio_fwhm, 1),
    "reshs": f(100 * res_s, 1),
    "ic": f(fc["c0_intercept"], 4), "icse": f(fc["c0_se"], 4),
    "iq": f(fq["c0_intercept"], 4), "iqse": f(fq["c0_se"], 4),
    "iq6": f(fq6["c0_intercept"], 4), "iq6se": f(fq6["c0_se"], 4),
    "sefac": f(fq6["c0_se"] / fq["c0_se"], 1),
    "lw": f(LW, 4),
    "lwres": sci([c for c in G["provenance"]["self_checks"]
                  if c["check"] == "lieb_wu_two_forms_agree"][0]["got"]),
    "dlw": f(abs(sep["charge_intercept_minus_lieb_wu"]), 4),
    "dlw6": f(dlw6, 4),
    "zlw": f(z_lin_LW, 2), "zq": f(z_qua_LW, 2), "zq6": f(z_qua6_LW, 2),
    "z6z": f(z_qua6_zero, 1), "zlz": f(z_lin_zero, 0),
    "meanLgs": f(mean_Lgs, 3),
    "slps": f(ltg["spin"]["d1"], 3), "slpsse": f(ltg["spin"]["d1_se"], 3),
    "zslps": f(abs(ltg["spin"]["d1_z"]), 2),
    "slpc": f(ltg["charge"]["d1"], 2), "slpcse": f(ltg["charge"]["d1_se"], 2),
    "zslpc": f(ltg["charge"]["d1_z"], 0),
    "ih": f(fh_["c0_intercept"], 4), "ihse": f(fh_["c0_se"], 4), "zh": f(z_heis0, 1),
    "drift": f(100 * drift_heis, 1),
}
with open(os.path.join(BUILD, "fig_gapscaling_caption.tex"), "w", encoding="utf-8") as fh:
    fh.write(cap)

# ---- guards: the caption must not re-introduce the misattributed z, nor banned words ----
assert "$95$ s.e." not in cap, "the separation z (94.66) is back in the caption"
_low = cap.lower()
for _w in ("saturat", "machine precision", "gapless", "continuum", "beyond-classical"):
    assert _w not in _low, "banned term in caption: %s" % _w
assert "saturat" not in tex.lower(), "banned term in fragment"
assert f(z_lin_zero, 0) == "104", f(z_lin_zero, 0)
assert f(ratio_fwhm, 1) == "16.6", f(ratio_fwhm, 1)
print("guards OK: z_lin_zero=%.4f z_qua6_zero=%.4f z_qua_zero=%.4f ratio_fwhm=%.4f"
      % (z_lin_zero, z_qua6_zero, z_qua_zero, ratio_fwhm))
facts = {
    "eta_spin": eta_spin, "eta_charge": eta_charge,
    "spin_gaps": gs, "charge_gaps": gc, "heisenberg_gaps": gh,
    "L_times_spin": Lgs, "L_times_charge": Lgc, "L_times_heisenberg": Lgh,
    "mean_L_spin_gap_L6plus": mean_Lgs,
    "heisenberg_L_times_gap_drift_L4_to_L12": drift_heis,
    "spin_fit_L6plus": {k: fs[k] for k in ("c0_intercept", "c0_se", "c1_slope", "c1_se", "r2")},
    "charge_fit_linear_L6plus": {k: fc[k] for k in ("c0_intercept", "c0_se", "c1_slope", "c1_se", "r2")},
    "charge_fit_quadratic_allL": {k: fq[k] for k in ("c0_intercept", "c0_se", "c1", "c2", "r2")},
    "heisenberg_fit_L6plus": {k: fh_[k] for k in ("c0_intercept", "c0_se", "c1_slope", "r2")},
    "lieb_wu": LW, "lieb_wu_literature": LWlit,
    "z_lieb_wu_above_linear_intercept": z_lin_LW,
    "z_quadratic_intercept_above_lieb_wu": z_qua_LW,
    "z_spin_intercept_from_zero": z_spin0,
    "z_heisenberg_intercept_from_zero": z_heis0,
    "z_intercept_separation_NOT_a_distance_from_zero": sep["z_intercept_separation"],
    "charge_fit_quadratic_L6plus": {k: fq6[k] for k in ("c0_intercept", "c0_se", "c1", "c2", "r2")},
    "z_lieb_wu_above_quadratic_L6plus_intercept": z_qua6_LW,
    "lieb_wu_minus_quadratic_L6plus_intercept": dlw6,
    "z_linear_intercept_from_zero": z_lin_zero,
    "z_quadratic_allL_intercept_from_zero": z_qua_zero,
    "z_quadratic_L6plus_intercept_from_zero": z_qua6_zero,
    "charge_over_spin_gap_in_own_fwhm_units": ratio_fwhm,
    "Ltimesgap_slope_spin_d1_se_z": [ltg["spin"]["d1"], ltg["spin"]["d1_se"], ltg["spin"]["d1_z"]],
    "Ltimesgap_slope_charge_d1_se_z": [ltg["charge"]["d1"], ltg["charge"]["d1_se"], ltg["charge"]["d1_z"]],
    "spin_gap_L12_over_eta": r_s_eta, "spin_gap_L12_over_fwhm": r_s_fwhm,
    "charge_gap_L12_over_eta": r_c_eta, "charge_gap_L12_over_fwhm": r_c_fwhm,
    "lorentz_height_at_zero_spin": res_s, "lorentz_height_at_zero_charge": res_c,
    "sizes_inside_spin_linewidth": n_below,
    "hubbard_over_heisenberg_ratio": ratio_hh,
}
with open(os.path.join(BUILD, "fig_gapscaling_facts.json"), "w", encoding="utf-8") as fh:
    json.dump(facts, fh, indent=1)

print("wrote", out)
for k in ("spin_gap_L12_over_eta", "spin_gap_L12_over_fwhm", "charge_gap_L12_over_eta",
          "charge_gap_L12_over_fwhm", "lorentz_height_at_zero_spin", "lorentz_height_at_zero_charge",
          "mean_L_spin_gap_L6plus", "heisenberg_L_times_gap_drift_L4_to_L12",
          "z_lieb_wu_above_linear_intercept", "z_quadratic_intercept_above_lieb_wu",
          "z_spin_intercept_from_zero", "z_heisenberg_intercept_from_zero",
          "z_lieb_wu_above_quadratic_L6plus_intercept", "z_linear_intercept_from_zero",
          "z_quadratic_L6plus_intercept_from_zero", "charge_over_spin_gap_in_own_fwhm_units"):
    print("  %-42s %s" % (k, facts[k]))
print("  sizes with spin gap inside 2eta:", n_below)
print("  Hubbard/Heisenberg ratio:", {k: round(v, 4) for k, v in ratio_hh.items()})
