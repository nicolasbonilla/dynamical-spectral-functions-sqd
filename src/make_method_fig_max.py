# -*- coding: utf-8 -*-
r"""Emit the MAX-LEVEL method-validation figure fragment (fig:method) + aw_method.dat, from method_max.json.
(a) A(w) exact vs bitstring-sampled, with the exact sum rule; (b) rel-L1 vs Krylov order K over 8 seeds
(error band = std), the two physical regimes (discovery / geometric collapse), and the geometric-convergence
line rho^{-K} that Theorem 1 guarantees. Fonts inherit the document body size (11pt) exactly."""
# --- DEPOSIT PATHS (repaired 2026-09-18, second pass) -----------------------
# This figure generator addresses every file it reads and writes by a path relative
# to the repository root ('data/...', 'paper/figs/...'), so it only ever worked when
# launched from that root and raised FileNotFoundError from anywhere else.  Rather
# than rewrite every literal -- which risks changing an output byte -- the process's
# working directory is anchored to the repository THIS FILE lives in.  Run from a
# copy of the tree (as src/check_figures.py does), it anchors to that copy, which is
# the behaviour that check wants.  Paths only; no figure content changed.
import os as _os
_os.chdir(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
# ---------------------------------------------------------------------------
import json, numpy as np
d=json.load(open('data/method_max.json'))
g=np.array(d['grid']); Aex=np.array(d['A_exact']); Asm=np.array(d['A_final'])
# panel (a) data file
with open('paper/figs/aw_method.dat','w') as f:
    f.write('w Aex Asm\n')
    for w,ae,asm in zip(g,Aex,Asm): f.write(f'{w:.6e} {ae:.6e} {asm:.6e}\n')
sr=d['sumrule']; rho=d['rho']; ks=d['kstar']
# panel (a) annotation: rel-L1 of the COMPLETE-SUBSPACE reconstruction actually plotted here
# (exact vs A_final), NOT the finite-shot convergence floor l1[-1] shown in panel (b).
l1_complete=float(np.sum(np.abs(Aex-Asm))/np.sum(np.abs(Aex)))
S=d['sampled']; K=[r['K'] for r in S]; l1=[r['l1'] for r in S]; sd=[r['l1_std'] for r in S]
Smean=[r['S'] for r in S]
FL=0.010
lo=" ".join(f"({r['K']},{max(r['l1']-r['l1_std'],FL):.4f})" for r in S)
hi=" ".join(f"({r['K']},{r['l1']+r['l1_std']:.4f})" for r in S)
pts=" ".join(f"({r['K']},{r['l1']:.4f})" for r in S)
fit=" ".join(f"({f['K']},{f['l1']:.4f})" for f in d['fit'])

L=r"""% native \input fragment (fig:method) — fonts = document body (11pt) exactly
\definecolor{inkPrim}{HTML}{1B1B1B}\definecolor{inkMute}{HTML}{6E6E6E}\definecolor{clExact}{HTML}{1F3A5F}
\definecolor{coral}{HTML}{F0653F}\definecolor{gridClr}{HTML}{ECEAE5}\definecolor{thm}{HTML}{2E7D5B}
\begin{tikzpicture}[font=\rmfamily]
\begin{groupplot}[
  group style={group size=2 by 1, horizontal sep=1.9cm},
  width=7.7cm, height=5.9cm,
  axis line style={inkPrim, line width=0.6pt}, tick style={inkPrim, line width=0.6pt},
  title style={yshift=1pt},
]
% ---- (a) A(w): exact vs bitstring-sampled ----
\nextgroupplot[
  xlabel={frequency \; $\omega/t$}, ylabel={$A(\omega)$},
  title={(a)\quad bitstring reconstruction of $A(\omega)$},
  xmin=3, xmax=14, ymin=0, ymax=0.50, ymajorgrids, y grid style={gridClr},
  legend style={draw=none, fill=white, fill opacity=0.75, text opacity=1, at={(0.985,0.985)}, anchor=north east}, legend cell align=left,
]
  \addplot[draw=none, fill=clExact, fill opacity=0.16] table[x=w,y=Aex]{aw_method.dat} \closedcycle;
  \addplot[clExact, line width=1.3pt] table[x=w,y=Aex]{aw_method.dat}; \addlegendentry{exact}
  \addplot[coral, line width=1.1pt, dash pattern=on 3pt off 2.2pt] table[x=w,y=Asm]{aw_method.dat}; \addlegendentry{bitstring-sampled}
  \node[anchor=north east, align=right, inkPrim] at (rel axis cs:0.985,0.56)
     {rel-$L_1=__L1F__$\\[-1pt] $\int\! A=__SR__$};
% ---- (b) convergence: two regimes + geometric law (Thm 1) ----
\nextgroupplot[
  xlabel={Krylov order \; $K$ \; (time-evolution steps)}, ylabel={rel-$L_1$ error},
  title={(b)\quad geometric convergence, verified},
  xmin=-0.4, xmax=12.4, ymode=log, ymin=0.013, ymax=2.6, ymajorgrids, y grid style={gridClr},
  ytick={0.02,0.05,0.1,0.2,0.5}, yticklabels={$2\%$,$5\%$,$10\%$,$20\%$,$50\%$},
  legend style={draw=none, fill=white, fill opacity=0.85, text opacity=1, at={(0.985,0.985)}, anchor=north east}, legend cell align=left,
]
  \fill[inkMute, fill opacity=0.06] (axis cs:-0.4,0.013) rectangle (axis cs:__KS__,2.6);
  \node[inkPrim, anchor=center, align=center] at (axis cs:1.75,0.19) {discovery\\ regime};
  \addplot[name path=lo, draw=none, forget plot] coordinates {__LO__};
  \addplot[name path=hi, draw=none, forget plot] coordinates {__HI__};
  \addplot[coral, fill opacity=0.16, forget plot] fill between[of=lo and hi];
  \addplot[thm, line width=1pt, dashed] coordinates {__FIT__}; \addlegendentry{geometric (Thm 1)}
  \addplot[coral, line width=1.3pt, mark=*, mark size=2.2pt,
           mark options={fill=coral, draw=white, line width=0.4pt}] coordinates {__PTS__}; \addlegendentry{sampled (8 seeds)}
  \node[inkPrim, anchor=south east] at (axis cs:12.2,0.033) {$\rho\approx__RHO__$};
\end{groupplot}
\end{tikzpicture}
"""
L=(L.replace('__SR__',f'{sr:.3f}').replace('__L1F__',f'{l1_complete:.3f}').replace('__RHO__',f'{rho:.1f}')
    .replace('__KS__',f'{ks}').replace('__LO__',lo).replace('__HI__',hi).replace('__FIT__',fit).replace('__PTS__',pts))
open('paper/figs/fig_method_native_frag.tex','w',encoding='utf-8').write(L)
print(f'wrote fragment + aw_method.dat | sumrule={sr:.4f} rho={rho:.2f} '
      f'complete-subspace rel-L1={l1_complete:.3f} (panel a) | finite-shot floor={l1[-1]:.3f} (panel b)')
