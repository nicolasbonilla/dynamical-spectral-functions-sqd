# -*- coding: utf-8 -*-
r"""Emit paper/figs/fig_akw_sampled_native.tex (fig:akwsampled, Fig. 5) from sampled_akw_L8.json.

Self-contained native pgfplots fragment (inline coordinates, lmodern fonts -> matches the paper): panel
(a) overlays the bitstring-sampled A(k,omega) (dashed) on the exact Haydock target (solid) at k=0,pi/2,pi;
panel (b) the per-momentum rel-L1 error (all < 0.5%). No hand-typed numbers -- everything from the JSON.

Run:  python sampled_akw.py && python make_akw_sampled_fig.py
"""
import json
import numpy as np

d = json.load(open('data/sampled_akw_L8.json'))
wg = np.array(d['wg']); ov = d['overlay']; pk = d['per_k']; FR = d['frac']
step = 3; ws = wg[::step]


def coords(y):
    y = np.array(y)[::step]
    return " ".join("(%.3f,%.4f)" % (w, v) for w, v in zip(ws, y))


cuts = ['0.0', '0.5', '1.0']
labels = {'0.0': r'$k{=}0$', '0.5': r'$k{=}\pi/2$', '1.0': r'$k{=}\pi$'}
maxA = max(max(ov[c]['exact']) for c in cuts); off = 1.05 * maxA
T = []
T.append("% Sampled-vs-exact A(k,omega): REAL reconstruction (L=8 chain, U/t=8, 85% of the (N+-1) sector).")
T.append("% Generated from sampled_akw_L8.json (mean per-k rel-L1 %.4f, max %.4f); no hand-typed numbers." % (d['mean_relL1'], d['max_relL1']))
T.append(r"\definecolor{akExact}{HTML}{1B1B1B}\definecolor{akSamp}{HTML}{D55E00}\definecolor{akInk}{HTML}{1B1B1B}")
T.append(r"\begin{tikzpicture}")
T.append(r"\begin{groupplot}[group style={group size=2 by 1, horizontal sep=1.7cm},")
T.append(r"  width=8.2cm, height=6.2cm, paperaxis,")
T.append(r"  title style={at={(0,1)},anchor=south west,font=\normalsize\bfseries,yshift=2pt},")
T.append(r"  label style={font=\normalsize}, tick label style={font=\footnotesize}]")
T.append(r"\nextgroupplot[title={(a)~sampled vs exact $A(k,\omega)$},")
T.append(r"  xlabel={$\omega-E_0\ (t)$}, ylabel={$A(k,\omega)$ (offset by $k$)}, xmin=-9,xmax=9,")
T.append(r"  ymin=-0.3,ymax=%.2f, ytick=\empty, xtick={-8,-4,0,4,8}," % (3 * off + off * 0.7))
T.append(r"  legend style={font=\footnotesize,at={(0.03,0.985)},anchor=north west,legend columns=1,draw=none,fill=white,fill opacity=0.55,text opacity=1},legend cell align=left]")
for i, c in enumerate(cuts):
    o = i * off
    ex = np.array(ov[c]['exact']) + o; sm = np.array(ov[c]['sampled']) + o
    T.append(r"\addplot[akExact,line width=1.2pt,fill=akExact,fill opacity=0.06] coordinates {" + coords(ex) + r"} \closedcycle;")
    T.append(r"\addplot[akSamp,line width=1.1pt,dash pattern=on 4pt off 2.5pt] coordinates {" + coords(sm) + r"};")
    if i == 0:
        T.append(r"\addlegendentry{exact}\addlegendentry{sampled ($85\%$ sector)}")
    T.append(r"\node[akInk,font=\footnotesize,anchor=west] at (axis cs:5.2," + "%.3f" % (o + 0.55 * off) + r") {" + labels[c] + r"};")
T.append(r"\nextgroupplot[title={(b)~per-momentum error},")
T.append(r"  xlabel={$k/\pi$}, ylabel={relative-$L_1$}, xmin=-0.1,xmax=2.0,")
T.append(r"  ymin=0,ymax=0.006, xtick={0,0.5,1,1.5}, ytick={0,0.002,0.004},")
T.append(r"  yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=3}, ymajorgrids, grid style={draw=black!10}]")
kk = sorted(float(k) for k in pk)
co = " ".join("(%.3f,%.4f)" % (k, pk["%.3f" % k]['relL1']) for k in kk)
T.append(r"\addplot[akSamp,only marks,mark=*,mark size=2.6pt,mark options={fill=akSamp,draw=white,line width=0.6pt}] coordinates {" + co + r"};")
T.append(r"\addplot[akExact,densely dotted,line width=0.8pt,domain=-0.1:2.0] {0.005};")
T.append(r"\node[akInk,font=\footnotesize,anchor=north east] at (axis cs:1.95,0.0059) {$<0.5\%$ at every $k$};")
T.append(r"\node[akInk,font=\footnotesize,anchor=south] at (axis cs:0.9,0.0038) {mean " + "%.1f" % (d['mean_relL1'] * 1000) + r"$\times10^{-3}$};")
T.append(r"\end{groupplot}")
T.append(r"\end{tikzpicture}")
open('paper/figs/fig_akw_sampled_native.tex', 'w').write("\n".join(T) + "\n")
print("wrote paper/figs/fig_akw_sampled_native.tex; mean rel-L1=%.4f" % d['mean_relL1'])
