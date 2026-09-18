# -*- coding: utf-8 -*-
r"""Native fig_gflownet (honest dequantization guard) from VERIFIED 5-seed results
(P1_review_ML_for_SQD/data/GFlowNet_SQD_resultados_verificados.md sec 3.4; N2 CAS(10e,12o), dim=120,
1000 shots, real Docker pyscf+torch). Big improvement = cheap CLASSICAL prior (ibm -> ibm+cheap);
the generative GFlowNet does NOT beat the fair classical control -> a learned generator cannot
manufacture the dominant determinants."""
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
# verified mHa, shots=1000, dim=120, 5 seeds: (label-lines, mean, std, class)
METH=[(r'raw',43.22,4.47,'base'),(r'ibm',39.24,3.95,'base'),
      (r'ibm\\ $+$cheap',26.98,0.66,'ctrl'),(r'gfn\\ (no cheap)',25.96,3.40,'gen'),
      (r'gfn\\ (fused)',22.94,3.30,'gen')]
COL={'base':'inkMute','ctrl':'blue','gen':'coral'}
def y(v): return v/10.0
bars=[]
for i,(n,v,s,c) in enumerate(METH):
    bars.append(
      f"  \\draw[fill={COL[c]},draw=none] ({i}-0.30,0) rectangle ({i}+0.30,{y(v):.3f});"
      f" \\draw[inkPrim,line width=0.8pt] ({i},{y(v-s):.3f})--({i},{y(v+s):.3f});"
      f" \\draw[inkPrim,line width=0.8pt] ({i}-0.08,{y(v-s):.3f})--({i}+0.08,{y(v-s):.3f});"
      f" \\draw[inkPrim,line width=0.8pt] ({i}-0.08,{y(v+s):.3f})--({i}+0.08,{y(v+s):.3f});"
      f" \\node[font=\\normalsize,inkPrim,anchor=south] at ({i},{y(v+s)+0.04:.3f}) {{${v:.1f}$}};"
      f" \\node[font=\\normalsize,inkPrim,anchor=north,align=center] at ({i},-0.08) {{{n}}};")
bars="\n".join(bars)

tex=r"""\documentclass[11pt,border=8pt]{standalone}
\usepackage[T1]{fontenc}\usepackage{lmodern}\usepackage{amsmath}
\usepackage{tikz}\usetikzlibrary{arrows.meta,decorations.pathreplacing}
\definecolor{inkPrim}{HTML}{1B1B1B}\definecolor{inkMute}{HTML}{9A968E}
\definecolor{coral}{HTML}{D55E00}\definecolor{blue}{HTML}{0072B2}\definecolor{gridClr}{HTML}{ECEAE5}
\definecolor{chem}{HTML}{2E7D5B}
\begin{document}\begin{tikzpicture}[font=\rmfamily,x=1.95cm,y=1.15cm]
\foreach \yv in {0,10,20,30,40,50,60}{
  \draw[gridClr,line width=0.5pt] (-0.55,\yv/10) -- (4.6,\yv/10);
  \node[font=\normalsize,inkPrim,anchor=east] at (-0.6,\yv/10) {\yv};}
\draw[inkPrim,line width=0.6pt] (-0.55,0) -- (-0.55,6.15);
\draw[inkPrim,line width=0.6pt] (-0.55,0) -- (4.6,0);
\node[font=\normalsize,inkPrim,rotate=90,anchor=south] at (-1.05,2.5) {energy error vs FCI \; (mHa)};
""" + bars + r"""
% chemical accuracy
\draw[chem,dashed,line width=0.7pt] (-0.55,0.16) -- (4.6,0.16);
\node[font=\normalsize,inkPrim,anchor=east] at (4.55,0.33) {chemical accuracy $1.6$};
% (i) the big gain is the cheap classical prior: label above ALL bars, arrow down to ibm+cheap
\node[font=\normalsize,align=left,anchor=south west,inkPrim] at (0.15,5.62)
  {big drop $=$ cheap \emph{classical} prior};
\draw[-{Stealth[length=2mm]},blue,line width=0.9pt] (1.00,5.42) to[bend left=16] (1.82,2.88);
% (ii) generative barely differs from the fair control: 2-line label above the 3 short right bars
\draw[decorate,decoration={brace,amplitude=4pt,mirror},coral,line width=0.7pt] (1.70,3.55) -- (4.30,3.55);
\node[font=\normalsize,align=center,anchor=south,inkPrim] at (3.05,3.68)
  {generative $\approx$ classical control\\ \emph{no clean win}};
\end{tikzpicture}\end{document}
"""
open('fig_gflow_native.tex','w',encoding='utf-8').write(tex)
print("wrote fig_gflow_native.tex")
