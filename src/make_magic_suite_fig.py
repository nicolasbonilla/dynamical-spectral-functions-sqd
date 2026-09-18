# -*- coding: utf-8 -*-
r"""Generate fig_magic_suite.tex (native TikZ/pgfplots dumbbell plot) directly from the REAL,
FCI-verified n19_suite.json. No hand-typed numbers. Normalized fermionic magic FAF/FAF_max at
equilibrium (open) -> dissociation (filled) for all 19 molecules, sorted by dissociation magic.
Design: vibrant regime colours in the DATA, neutral dark text, polished dots, subtle row bands."""
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
import json
d=json.load(open('data/n19_suite.json'))
M=[m for m in d['mols'] if 'error' not in m]
for m in M:
    m['fe']=m['FAF_eq']/m['fafmax']; m['fd']=m['FAF_diss']/m['fafmax']
M.sort(key=lambda m:m['fd'])           # lowest at bottom (y=0), highest at top

COV={'H2','HF','F2','CO','N2','H2O','NH3','BeH2','H4','H6'}
MRF={'C2','O2','OH','CN','NO'}          # inherently multireference / open-shell
ION={'LiF','LiH','BH','BeH'}            # ionic / weakly-correlated
def grp(x): return 'cov' if x in COV else ('mrf' if x in MRF else 'ion')
COL={'cov':'covOI','mrf':'mrfOI','ion':'ionOI'}
def texlab(s):
    for a,b in [('2','$_2$'),('3','$_3$'),('4','$_4$'),('6','$_6$')]: s=s.replace(a,b)
    return s

rows=[]
for i,m in enumerate(M):
    g=grp(m['mol']); c=COL[g]
    rows.append((i,m['mol'],g,c,m['fe'],m['fd']))

body=[]
# subtle alternating row bands (elegance + readability), behind everything
for (i,mol,g,c,fe,fd) in rows:
    if i%2==0:
        body.append(f"  \\fill[band] (axis cs:-0.02,{i-0.5}) rectangle (axis cs:1.08,{i+0.5});")
# connecting arrows + polished dots (open eq -> filled diss)
for (i,mol,g,c,fe,fd) in rows:
    body.append(f"  \\draw[{c},line width=2.4pt,-{{Stealth[length=2.8mm]}}] (axis cs:{fe:.4f},{i}) -- (axis cs:{fd:.4f},{i});")
    body.append(f"  \\addplot[only marks,mark=*,mark size=2.6pt,draw={c},fill=white,line width=1.2pt] coordinates {{({fe:.4f},{i})}};")
    body.append(f"  \\addplot[only marks,mark=*,mark size=3.8pt,draw=white,fill={c},line width=0.8pt] coordinates {{({fd:.4f},{i})}};")

ylabels=", ".join(texlab(mol) for (_,mol,*_ ) in rows)
yticks=", ".join(str(i) for i in range(len(rows)))
def rowof(name): return next(r[0] for r in rows if r[1]==name)
c2i=rowof('C2'); c2fe=next(r[4] for r in rows if r[1]=='C2')

tmpl=r"""\documentclass[border=8pt]{standalone}
\usepackage[T1]{fontenc}\usepackage{lmodern}\usepackage{amsmath,amssymb}
\usepackage{tikz}\usetikzlibrary{arrows.meta}
\usepackage{pgfplots}\pgfplotsset{compat=1.18}
\definecolor{inkPrim}{HTML}{1B1B1B}\definecolor{inkMute}{HTML}{6E6E6E}
\definecolor{covOI}{HTML}{E8543A}   % vibrant vermilion - covalent bond-breaking
\definecolor{mrfOI}{HTML}{10A87B}   % vibrant emerald   - multireference / open-shell
\definecolor{ionOI}{HTML}{2E6FE0}   % vibrant blue      - ionic / weakly correlated
\definecolor{gridClr}{HTML}{E7E4DE}\definecolor{band}{HTML}{F5F2EC}
\begin{document}
\begin{tikzpicture}[font=\rmfamily]
\begin{axis}[
  width=13.6cm, height=13.8cm, clip=false,
  xmin=-0.02, xmax=1.08, ymin=-0.7, ymax=__YMAX__,
  xtick={0,0.25,0.5,0.75,1.0}, xticklabels={$0$,$0.25$,$0.5$,$0.75$,$1$},
  ytick={__YTICKS__}, yticklabels={__YLABELS__},
  ytick style={draw=none}, y tick label style={font=\normalsize\bfseries, inkPrim},
  xmajorgrids, major grid style={gridClr, line width=0.5pt},
  axis line style={inkPrim, line width=0.6pt},
  xtick style={inkPrim}, x tick label style={font=\small, inkPrim},
  xlabel={\textbf{fermionic magic}\ \ $\mathcal F_1/\mathcal F_1^{\max}$ \quad\small\textcolor{inkPrim}{$0=$ Gaussian\, $\to$\, $1=$ maximal non-Gaussianity}},
  xlabel style={font=\normalsize, inkPrim}, ylabel={}, axis on top=false,
]
% vertical reference at 0
\draw[inkMute,line width=0.5pt,dotted] (axis cs:0,-0.7) -- (axis cs:0,__YMAX__);
__BODY__
% --- legend: neutral text + colour swatch (bottom-right empty triangle) ---
\node[anchor=west,font=\small,inkPrim] at (axis cs:0.47,3.35)
   {\textcolor{covOI}{$\blacksquare$}\ \ \textbf{covalent}: magic switches \emph{on}};
\node[anchor=west,font=\small,inkPrim] at (axis cs:0.47,2.45)
   {\textcolor{mrfOI}{$\blacksquare$}\ \ \textbf{multiref.\,/\,open-shell}: large already};
\node[anchor=west,font=\small,inkPrim] at (axis cs:0.47,1.55)
   {\textcolor{ionOI}{$\blacksquare$}\ \ \textbf{ionic\,/\,weak}: stays near-Gaussian};
% eq/diss marker key
\node[anchor=west,font=\footnotesize,inkPrim] at (axis cs:0.47,0.55)
   {\tikz{\draw[covOI,line width=1.2pt] (0,0) circle (2.2pt);}\,equilibrium \ \ $\rightarrow$ \ \ \tikz{\fill[covOI] (0,0) circle (3.4pt);}\,dissociation};
% --- callout: C2 already multireference at equilibrium ---
\node[anchor=south west,font=\small\itshape,inkPrim,align=left] at (axis cs:__C2X__,__C2Y__)
   {already multireference\\ at equilibrium};
\draw[inkMute,line width=0.6pt,-{Stealth[length=1.6mm]}] (axis cs:__C2AX__,__C2AY__) -- (axis cs:__C2FE__,__C2I__);
\end{axis}
\end{tikzpicture}
\end{document}
"""
N=len(rows)
out=(tmpl.replace('__YMAX__',f'{N-0.3:.1f}')
        .replace('__YTICKS__',yticks).replace('__YLABELS__',ylabels)
        .replace('__BODY__',"\n".join(body))
        .replace('__C2X__',f'{c2fe+0.03:.3f}').replace('__C2Y__',f'{c2i+0.22:.2f}')
        .replace('__C2AX__',f'{c2fe+0.02:.3f}').replace('__C2AY__',f'{c2i+0.22:.2f}')
        .replace('__C2FE__',f'{c2fe:.4f}').replace('__C2I__',f'{c2i}'))
open('fig_magic_suite.tex','w',encoding='utf-8').write(out)
print("wrote fig_magic_suite.tex  rows=",N)
