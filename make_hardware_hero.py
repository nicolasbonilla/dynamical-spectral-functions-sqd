# -*- coding: utf-8 -*-
r"""HARDWARE figure (fig:heron), TWO panels, emitted as an \input FRAGMENT so \normalsize == the paper
body size EXACTLY (no include-scaling font games). Rules enforced: (1) BOTH real hardware runs shown;
(2) every font = body \normalsize; (3) all text neutral (inkPrim/inkMute) -- colour only in data;
(4) all text inside the axes, never over a curve. Panels:
(a) ibm_fez -- Hubbard L=6 spectral function A(omega) (full-sector coverage -> exact).
(b) ibm_marrakesh -- N2 ground-state energy vs recovery step: noisy device reaches 0.59 mHa (below
    chemical accuracy) and beats the noiseless statevector simulation of the same circuit (29.5 mHa).
Data: heron_spectral.json, hw_lucj_n2_result.json."""
import json, numpy as np
h=json.load(open('heron_spectral.json')); g=np.array(h['grid'])
A=np.array(h['A_hw'])   # the REAL ibm_fez hardware-sampled reconstruction (==A_exact to 8e-15 by full-sector coverage)
with open('paper/figs/heron_hot.dat','w') as f:
    f.write('w A\n')
    for w,a in zip(g,A): f.write(f'{w:.4f} {a:.5f}\n')
Amax=float(A.max())
n=json.load(open('hw_lucj_n2_result.json')); fci=n['e_fci']
hw=[abs(e-fci)*1000 for e in n['hist_hw']]; sim=[abs(e-fci)*1000 for e in n['hist_sim']]
hwc=" ".join(f"({i},{v:.3f})" for i,v in enumerate(hw))
simc=" ".join(f"({i},{v:.3f})" for i,v in enumerate(sim))

frag=r"""% fragment: \input into the figure environment -> \normalsize equals the body font exactly.
\definecolor{inkPrim}{HTML}{1B1B1B}\definecolor{inkMute}{HTML}{6E6E6E}
\definecolor{warmFill}{HTML}{E8770C}\definecolor{warmLine}{HTML}{7A1E10}
\definecolor{hwWarm}{HTML}{C4360C}\definecolor{simSlate}{HTML}{9AA1AC}\definecolor{band}{HTML}{ECE7DF}
\begin{tikzpicture}[font=\rmfamily]
\pgfplotsset{hpanel/.style={width=8.5cm,height=6.0cm,axis line style={inkPrim,line width=0.6pt},
  tick style={inkPrim,line width=0.6pt},tick label style={font=\normalsize},label style={font=\normalsize},
  title style={font=\normalsize,yshift=1pt},ymajorgrids,major grid style={inkMute,opacity=0.10,line width=0.3pt}}}
% ===================== (a) ibm_fez : spectral function A(omega) =====================
% CLEAN: curve + title only; all technical detail (shots, mitigation, coverage) lives once, in the caption.
\begin{axis}[name=axa,hpanel,
  title={(a)\ \ \texttt{ibm\_fez} $\cdot$ 12 qubits: \ $A(\omega)$},
  xlabel={frequency \ $\omega-E_0$}, ylabel={$A(\omega)$},
  xmin=__XMIN__,xmax=__XMAX__,ymin=0,ymax=__YMAXA__]
  \addplot[draw=none,fill=warmFill,fill opacity=0.28] table[x=w,y=A]{figs/heron_hot.dat} \closedcycle;
  \addplot[warmLine,line width=1.3pt] table[x=w,y=A]{figs/heron_hot.dat};
\end{axis}
% ===================== (b) ibm_marrakesh : N2 ground-state energy =====================
% CLEAN: two curves + which-is-which + the two headline numbers; specs/jobs live once, in the caption.
\begin{axis}[name=axb,hpanel,at={(axa.south east)},anchor=south west,xshift=1.4cm,
  title={(b)\ \ \texttt{ibm\_marrakesh} $\cdot$ 24 qubits: \ $\mathrm{N_2}$ energy},
  xlabel={configuration-recovery step}, ylabel={error to FCI \ (mHa)},
  ymode=log,xmin=-0.35,xmax=4.55,ymin=0.7,ymax=90,
  xtick={0,1,2,3,4},ytick={1,3,10,30},yticklabels={$1$,$3$,$10$,$30$}]
  \fill[band] (axis cs:-0.35,0.7) rectangle (axis cs:4.55,1.6);
  \node[anchor=south west,font=\normalsize,inkPrim] at (axis cs:-0.28,0.76) {chem.\ accuracy};
  \addplot[simSlate,line width=1.4pt,dash pattern=on 5pt off 3pt,mark=square,mark size=2.6pt,
           mark options={draw=simSlate,fill=white,line width=1pt}] coordinates {__SIMC__};
  \addplot[hwWarm,line width=2.2pt,mark=*,mark size=3.2pt,mark options={fill=hwWarm,draw=white,line width=0.7pt}]
     coordinates {__HWC__};
  % which curve is which: a leader from each label lands ON its own curve; plus the two headline numbers
  \node[anchor=south,font=\normalsize,inkPrim] (nl) at (axis cs:1.85,44)
    {noiseless (same circuit), $29.5$~mHa};
  \draw[-{Stealth[length=1.7mm]},inkPrim,line width=0.5pt] (nl.south) -- (axis cs:1.55,28.9);
  \node[anchor=south,font=\normalsize,inkPrim] (nh) at (axis cs:3.05,4.0) {noisy hardware};
  \draw[-{Stealth[length=1.7mm]},inkPrim,line width=0.5pt] (nh.south) -- (axis cs:3.05,1.72);
  \node[anchor=south east,font=\normalsize,inkPrim] at (axis cs:4.5,1.72) {$0.59$~mHa};
\end{axis}
\end{tikzpicture}
"""
out=(frag.replace('__XMIN__',f'{g[0]:.2f}').replace('__XMAX__',f'{g[-1]:.2f}')
        .replace('__YMAXA__',f'{Amax*1.10:.3f}').replace('__SIMC__',simc).replace('__HWC__',hwc))
open('paper/figs/fig_hardware_hero_frag.tex','w',encoding='utf-8').write(out)
print('wrote fig_hardware_hero_frag.tex (fragment, 2 panels)')
