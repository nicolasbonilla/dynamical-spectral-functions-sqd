# -*- coding: utf-8 -*-
r"""Trim the 19 PyMOL orbital PNGs to content and emit a NATIVE TikZ gallery (fig_gallery_native.tex)
with the paper's lmodern font and shared identity (coral accent on the fermionic-magic axis, colourblind-safe
class colours matching Fig.~molsuite, readable type). Magic shown as a per-molecule bar (length = F1/F1max,
colour = bonding regime). Data from n19_suite.json (FCI-verified). Runs natively (relative paths)."""
import json, os
from PIL import Image

os.makedirs('gtrim', exist_ok=True)
d={m['mol']:m for m in json.load(open('data/n19_suite.json'))['mols'] if 'FAF_diss' in m}
order=sorted(d, key=lambda m:d[m]['FAF_diss']/d[m]['fafmax'])
COV={'H2','HF','F2','CO','N2','H2O','NH3','BeH2','H4','H6'}; MRF={'C2','O2','OH','CN','NO'}
def reg(m): return 'cov' if m in COV else ('mrf' if m in MRF else 'ion')
def tl(s):
    import re; return re.sub(r'(\d+)', r'$_{\1}$', s)

# --- trim each render to its non-transparent bounding box (+ padding) ---
for lab in order:
    im=Image.open(f'cubes/{lab}_r.png').convert('RGBA')
    bb=im.getbbox()
    if bb:
        p=10; bb=(max(0,bb[0]-p),max(0,bb[1]-p),min(im.width,bb[2]+p),min(im.height,bb[3]+p))
        im=im.crop(bb)
    im.save(f'gtrim/{lab}.png')

# --- layout ---
colx=[0.0,3.35,6.70,10.05,13.40]; rowy=[0.0,-3.95,-7.90,-11.85]
IW,IH=2.95,2.18                      # image bounding box (keepaspectratio)
cells=[]
for k,lab in enumerate(order):
    r,c=divmod(k,5); cx,cy=colx[c],rowy[r]; fd=d[lab]['FAF_diss']/d[lab]['fafmax']; g=reg(lab)
    s=[]
    s.append(f"  \\node[inner sep=0] (im{k}) at ({cx:.2f},{cy:.2f}) "
             f"{{\\includegraphics[width={IW}cm,height={IH}cm,keepaspectratio]{{gtrim/{lab}.png}}}};")
    s.append(f"  \\node[font=\\large\\bfseries, inkprim] at ({cx:.2f},{cy-1.48:.2f}) {{{tl(lab)}}};")
    # magic bar: track (light) + fill (regime colour), length proportional to fd
    x0,x1=cx-1.02,cx+1.02; yb=cy-1.90
    s.append(f"  \\fill[bartrack,rounded corners=0.6pt] ({x0:.2f},{yb-0.080:.2f}) rectangle ({x1:.2f},{yb+0.080:.2f});")
    s.append(f"  \\fill[{g},rounded corners=0.6pt] ({x0:.2f},{yb-0.080:.2f}) rectangle ({x0+2.04*fd:.3f},{yb+0.080:.2f});")
    s.append(f"  \\node[font=\\normalsize, inkprim] at ({cx:.2f},{cy-2.32:.2f}) {{$\\mathcal F_1/\\mathcal F_1^{{\\max}}={fd:.2f}$}};")
    cells.append("\n".join(s))

# legend in the empty 20th cell (row 3, col 4)
lx,ly=colx[4],rowy[3]
legend=rf"""
  \node[anchor=north west,align=left,font=\normalsize] at ({lx-1.95:.2f},{ly+1.15:.2f}) {{%
    \tikz{{\fill[phasep] (0,0) rectangle (0.34,0.24);}}~orbital phase $+$\\[3pt]
    \tikz{{\fill[phasen] (0,0) rectangle (0.34,0.24);}}~orbital phase $-$\\[6pt]
    \textcolor{{cov}}{{$\blacksquare$}}~\textbf{{covalent}}\\[3pt]
    \textcolor{{mrf}}{{$\blacksquare$}}~\textbf{{multiref. / open-shell}}\\[3pt]
    \textcolor{{ion}}{{$\blacksquare$}}~\textbf{{ionic / weakly corr.}}\\[6pt]
    \tikz{{\fill[bartrack](0,0)rectangle(0.75,0.18);\fill[inkprim](0,0)rectangle(0.48,0.18);}}~bar $=\mathcal F_1/\mathcal F_1^{{\max}}$
  }};"""

# bottom "single-reference -> multireference" axis: the paper's coral accent, the KEY physics direction
ax0,ax1,ay=colx[0]-1.4,colx[4]+1.4,rowy[3]-3.55
arrow=rf"""
  \draw[-{{Stealth[length=3.4mm]}},line width=1.7pt,coral] ({ax0:.2f},{ay:.2f}) -- ({ax1:.2f},{ay:.2f});
  \node[anchor=north west,font=\normalsize\itshape,inkprim] at ({ax0:.2f},{ay-0.16:.2f}) {{single-reference}};
  \node[anchor=north east,font=\normalsize\itshape,inkprim] at ({ax1:.2f},{ay-0.16:.2f}) {{multireference}};
  \node[anchor=south,font=\large\bfseries,inkprim] at ({(ax0+ax1)/2:.2f},{ay+0.16:.2f}) {{increasing fermionic magic $\mathcal F_1$}};"""

tex=r"""\documentclass[11pt,border=8pt]{standalone}
\usepackage[T1]{fontenc}\usepackage{lmodern}\usepackage{amsmath,amssymb}
\usepackage{graphicx}\usepackage{tikz}\usetikzlibrary{arrows.meta}
\definecolor{inkprim}{HTML}{1B1B1B}\definecolor{inkmute}{HTML}{6E6E6E}\definecolor{coral}{HTML}{F0653F}
\definecolor{cov}{HTML}{E8543A}\definecolor{mrf}{HTML}{10A87B}\definecolor{ion}{HTML}{2E6FE0}
\definecolor{bartrack}{HTML}{E7E4DE}
\definecolor{phasep}{HTML}{E8770C}\definecolor{phasen}{HTML}{1F6FB2}
\begin{document}\begin{tikzpicture}[font=\rmfamily]
""" + "\n".join(cells) + "\n" + legend + "\n" + arrow + r"""
\end{tikzpicture}\end{document}
"""
open('fig_gallery_native.tex','w',encoding='utf-8').write(tex)
print("wrote fig_gallery_native.tex ; trimmed",len(order),"pngs")
