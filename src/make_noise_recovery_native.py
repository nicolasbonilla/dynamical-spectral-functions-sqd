# -*- coding: utf-8 -*-
r"""Emit paper/figs/fig_noise_recovery_native.tex: SELF-CONTAINED native pgfplots fragment for the
molecular noise-robustness figure. (a) energy error vs bit-flip rate for 6 representative molecules
(S-CoRe solid, naive dashed), log-y, chemical-accuracy band; (b) full 19-molecule suite at eps=2%
(S-CoRe bars, naive ticks). Real data: molecular_noise_sweep.json + molecular_noise.json."""
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
sweep={m['mol']:m for m in json.load(open('data/molecular_noise_sweep.json'))['mols']}
suite=sorted(json.load(open('data/molecular_noise.json'))['mols'], key=lambda x:x['score'])
CHEM=1.6; FL=3e-3
show=[('H6','#6A4C93'),('NH3','#E8880C'),('N2','#D7263D'),('CO','#12876F'),('C2','#2E6F95')]

def coords(m,key):
    return " ".join(f"({r['eps']*100:.0f},{max(r[key],FL):.4f})" for r in sweep[m]['rows'])
def band(m):
    # log-scale-friendly +-1sigma band: floor the lower edge at 45% of the mean so it never collapses
    lo=" ".join(f"({r['eps']*100:.0f},{max(r['score']-r['score_std'],0.45*r['score'],FL):.4f})" for r in sweep[m]['rows'])
    hi=" ".join(f"({r['eps']*100:.0f},{max(r['score']+r['score_std'],FL):.4f})" for r in sweep[m]['rows'])
    return lo,hi

# panel (b): sorted molecules
blabels=[m['mol'] for m in suite]
bpos=list(range(len(suite)))
FL2=2.2e-3
cl=lambda v: max(v,FL2)
dumb=" ".join(f"\\draw[black!32,line width=0.8pt] (axis cs:{cl(m['naive']):.4f},{i}) -- (axis cs:{cl(m['score']):.4f},{i});"
              for i,m in enumerate(suite))
naive_pts=" ".join(f"({cl(m['naive']):.4f},{i})" for i,m in enumerate(suite))
score_pts=" ".join(f"({cl(m['score']):.4f},{i})" for i,m in enumerate(suite))
ytk=",".join(str(i) for i in bpos)
ylb=",".join(bl.replace('2','$_2$').replace('3','$_3$').replace('4','$_4$').replace('6','$_6$') for bl in blabels)

lines=[]
lines.append(r"% self-contained native pgfplots fragment (fig:noiserec)")
for lab,c in show: lines.append(f"\\definecolor{{nr{lab}}}{{HTML}}{{{c[1:]}}}")
lines.append(r"\definecolor{nrOK}{HTML}{0E7C7B}\definecolor{nrBad}{HTML}{5E4B8B}\definecolor{nrNv}{HTML}{D1495B}\definecolor{nrChem}{HTML}{2E7D32}")
lines.append(r"\begin{tikzpicture}")
# `group name=nrg` is load-bearing, not cosmetic.  pgfplots names the panels of a
# groupplot <group name> c<col>r<row>, and panel (b) positions itself against panel
# (a) through that node.  Giving the first \\nextgroupplot an explicit name=
# REPLACES that automatic node: panel (b) then finds no anchor and is composed on top
# of panel (a) (that is exactly what happened, and it reached a rendered page).  An
# explicit GROUP name leaves the automatic panel nodes intact and still gives the
# out-of-axis annotation below something stable to hang from.
lines.append(r"\begin{groupplot}[group style={group size=2 by 1, group name=nrg, horizontal sep=1.45cm}, height=9.5cm, paperaxis,")
lines.append(r"  title style={at={(0,1)},anchor=south west,font=\bfseries,yshift=1pt}]")
# ---- (a) ----
# 0.505 + 0.425 (was 0.54 + 0.46): inside a figure* the linewidth is 510pt, and the
# widened horizontal sep has to come out of the panels, not out of the margin.
lines.append(r"\nextgroupplot[width=0.505\linewidth, ymode=log, xmin=-1, xmax=31, ymin=3e-3, ymax=80,")
lines.append(r"  xlabel={per-qubit bit-flip rate $\varepsilon$ (\%)}, ylabel={energy error vs FCI (mHa)}, title={(a)},")
# The legend used to sit at the top left (0.03,0.97).  Measured, not eyeballed: an A/B
# render at 600 dpi with the legend fill made transparent recovered 590 ink pixels of
# nrH6 (#6A4C93) under the box -- the H6 naive curve, hidden by the legend at its own
# bottom-right corner.  The floor of the panel is the C2 band lower edge, which never
# comes below 0.024 mHa for eps in [6,24] while the axis starts at 3e-3, so the bottom
# centre has ~49pt of clear height for a 2-row box.  The same A/B probe recovers 0 px
# there.  Panel (b) was probed the same way and already hid nothing.
lines.append(r"  legend style={at={(0.5,0.03)},anchor=south,font=\normalsize,draw=black!20},legend columns=3,legend cell align=left]")
lines.append(rf"\draw[nrChem,dashed,line width=0.6pt] (axis cs:-1,{CHEM}) -- (axis cs:31,{CHEM});")
# The label used to sit INSIDE the plot area with a white fill at 70% opacity, which
# washed out the lower edge of the N2 band and the upper edge of the CO band.  Two
# in-plot placements were then tried and both rejected by measurement (600 dpi render
# differenced against a no-label render): near y=1.6 the panel has no ink-free window
# anywhere, and the nearest one at all is 0.7 decades below the rule, reachable only
# by a leader that crosses the CO band.  It is emitted outside the axis instead, at
# the end of this file; the rule itself runs the full width again.
lines.append(r"% chemical-accuracy label: emitted after the groupplot, see below.")
for lab,c in show:
    lo,hi=band(lab)
    lines.append(f"\\addplot[name path={lab}lo,draw=none,forget plot] coordinates {{{lo}}};")
    lines.append(f"\\addplot[name path={lab}hi,draw=none,forget plot] coordinates {{{hi}}};")
    lines.append(f"\\addplot[nr{lab},opacity=0.15,forget plot] fill between[of={lab}lo and {lab}hi];")
for lab,c in show:
    lines.append(f"\\addplot[nr{lab},line width=1.4pt,mark=*,mark size=1.3pt,mark options={{fill=nr{lab},draw=white,line width=0.4pt}}] coordinates {{{coords(lab,'score')}}}; \\addlegendentry{{{lab.replace('2','$_2$').replace('3','$_3$').replace('6','$_6$')}}}")
    lines.append(f"\\addplot[nr{lab},line width=0.9pt,densely dashed,forget plot] coordinates {{{coords(lab,'naive')}}};")
# ---- (b) ----
lines.append(r"\nextgroupplot[width=0.425\linewidth, xmode=log, xmin=2.2e-3, xmax=14, enlarge y limits=0.03,")
lines.append(rf"  ytick={{{ytk}}}, yticklabels={{{ylb}}}, y tick label style={{font=\normalsize}},")
lines.append(r"  xtick={0.01,0.1,1,10}, xticklabels={$0.01$,$0.1$,$1$,$10$}, xminorticks=false,")
lines.append(r"  xlabel={energy error at $\varepsilon=2\%$ (mHa)}, title={(b)},")
lines.append(r"  legend style={at={(0.03,0.03)},anchor=south west,font=\normalsize,draw=black!20},legend cell align=left]")
lines.append(rf"\fill[nrChem,fill opacity=0.06] (axis cs:2.2e-3,-0.6) rectangle (axis cs:{CHEM},{len(suite)-0.4});")
lines.append(rf"\draw[nrChem,dashed,line width=0.7pt] (axis cs:{CHEM},-0.6) -- (axis cs:{CHEM},{len(suite)-0.4});")
lines.append(rf"\node[nrChem,font=\tiny,anchor=south,rotate=90] at (axis cs:{CHEM},{len(suite)/2}) {{}};")
lines.append(" ".join(dumb.split()) if False else dumb)
lines.append(rf"\addplot[only marks,mark=o,mark size=1.9pt,nrNv,line width=0.8pt] coordinates {{{naive_pts}}}; \addlegendentry{{naive}}")
lines.append(rf"\addplot[only marks,mark=*,mark size=2.1pt,nrOK] coordinates {{{score_pts}}}; \addlegendentry{{S-CoRe}}")
lines.append(rf"\node[black,font=\normalsize,anchor=north east] at (axis cs:{CHEM*0.92},{len(suite)-0.7}) {{$1.6$ mHa}};")
lines.append(r"\end{groupplot}")
# OUTSIDE the axis: pgfplots clips nodes to the plot box.  It hangs off the automatic
# groupplot node, which is why the group carries an explicit `group name` above.
lines.append(r"\node[nrChem,font=\footnotesize,anchor=south east,inner sep=1pt]")
lines.append(r"      at (nrg c1r1.north east) {chemical accuracy ($1.6$\,mHa)};")
lines.append(r"\end{tikzpicture}")
open('paper/figs/fig_noise_recovery_native.tex','w',encoding='utf-8').write("\n".join(lines)+"\n")
print("wrote fig_noise_recovery_native.tex")
