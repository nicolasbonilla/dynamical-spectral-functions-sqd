# -*- coding: utf-8 -*-
r"""Emit paper/figs/fig_decoupling_native.tex: a SELF-CONTAINED native pgfplots fragment (inline
coordinates, no external .dat) for the two-sided decoupling figure. Data: n19_suite.json (molecules,
FAF vs log2|S|), cost_vs_ent.json (Hubbard U-sweep), verified free-fermion |S|. Fonts inherit lmodern."""
import json, numpy as np
mols=json.load(open('n19_suite.json'))['mols']
COV={'H2','HF','F2','CO','N2','H2O','NH3','BeH2','H4','H6'}; MRF={'C2','O2','OH','CN','NO'}
reg=lambda mol:'cov' if mol in COV else ('mrf' if mol in MRF else 'ion')
mF=[]; yS=[]; byreg={'cov':[],'mrf':[],'ion':[]}
for m in mols:
    for f,s in [(m['FAF_eq'],max(m['n99_eq'],1)),(m['FAF_diss'],max(m['n99_diss'],1))]:
        y=np.log2(s); mF.append(f); yS.append(y); byreg[reg(m['mol'])].append((f,y))
mF=np.array(mF); yS=np.array(yS); r=np.corrcoef(mF,yS)[0,1]
a,b=np.polyfit(mF,yS,1)                                  # trend line
regpts={k:" ".join(f"({f:.3f},{y:.3f})" for f,y in v) for k,v in byreg.items()}
xlo,xhi=mF.min(),mF.max()

rows=json.load(open('cost_vs_ent.json'))['rows']
def series(L,N): return sorted([(x['FAF'],np.log2(x['S'])) for x in rows if x['L']==L and x['N']==N])
hub={(6,6):'#D1495B',(8,8):'#009E73',(8,6):'#E0A21E'}
hubcoord={k:" ".join(f"({f:.3f},{s:.3f})" for f,s in series(*k)) for k in hub}

freeL=[4,6,8,10]; freeS=[35,336,3496,37361]
freecoord=" ".join(f"({L},{np.log2(S):.3f})" for L,S in zip(freeL,freeS))

tex=r"""% self-contained native pgfplots fragment (fig:decoupling)
\definecolor{dcH1}{HTML}{D1495B}\definecolor{dcH2}{HTML}{009E73}\definecolor{dcH3}{HTML}{E0A21E}
\definecolor{dcFree}{HTML}{3B7A57}\definecolor{dcCov}{HTML}{E8543A}\definecolor{dcMrf}{HTML}{10A87B}
\definecolor{dcIon}{HTML}{2E6FE0}\definecolor{dcTrend}{HTML}{9A9A9A}
\begin{tikzpicture}
\begin{groupplot}[group style={group size=3 by 1, horizontal sep=0.9cm},
  width=6.0cm, height=5.4cm, paperaxis,
  title style={at={(0.5,1.0)},anchor=south,font=\bfseries,yshift=2pt},
  label style={font=\normalsize}, tick label style={font=\normalsize},]
% ---- (a) molecules (coloured by bonding regime) ----
\nextgroupplot[xlabel={$\mathcal{F}_1=4\,\mathrm{tr}[\gamma(1{-}\gamma)]$},ylabel={$\log_2|\mathcal{S}|$},
  title={(a)~molecules},xmin=-0.6,ymin=-0.4,enlarge x limits=0.05]
\addplot[dcTrend,line width=1pt,dashed,domain=__XLO__:__XHI__,forget plot] {__A__*x+__B__};
\addplot[only marks,mark=*,mark size=2pt,draw=white,line width=0.4pt,fill=dcCov] coordinates {__COV__};
\addplot[only marks,mark=*,mark size=2pt,draw=white,line width=0.4pt,fill=dcMrf] coordinates {__MRF__};
\addplot[only marks,mark=*,mark size=2pt,draw=white,line width=0.4pt,fill=dcIon] coordinates {__ION__};
\node[anchor=south east] at (rel axis cs:0.98,0.02) {$r=__R__$};
% ---- (b) Hubbard chains ----
\nextgroupplot[xlabel={$\mathcal{F}_1$\; ($U\!\uparrow$)},ylabel={$\log_2|\mathcal{S}|$},
  title={(b)~Hubbard chains},legend style={at={(0.97,0.97)},anchor=north east,draw=black!20},
  legend cell align=left]
\addplot[dcH1,mark=*,mark size=2pt,line width=1.1pt] coordinates {__H66__}; \addlegendentry{$(6,6)$}
\addplot[dcH2,mark=square*,mark size=2pt,line width=1.1pt] coordinates {__H88__}; \addlegendentry{$(8,8)$}
\addplot[dcH3,mark=triangle*,mark size=2.4pt,line width=1.1pt] coordinates {__H86__}; \addlegendentry{$(8,6)$}
% ---- (c) free fermions (y-label dropped; shares the |S| meaning, own scale) ----
\nextgroupplot[xlabel={chain length $L$},ylabel={},
  title={(c)~free fermions},xtick={4,6,8,10},enlarge x limits=0.16,ymin=4,ymax=17.5]
\addplot[dcFree,mark=square*,mark size=2.4pt,line width=1.2pt] coordinates {__FREE__};
\node[anchor=south east,inner sep=1.5pt] at (axis cs:4,5.13){$35$};
\node[anchor=south east,inner sep=1.5pt] at (axis cs:6,8.39){$336$};
\node[anchor=south east,inner sep=1.5pt] at (axis cs:8,11.77){$3496$};
\node[anchor=north east,inner sep=1.5pt] at (axis cs:10,15.19){$37361$};
\end{groupplot}
\end{tikzpicture}
"""
tex=(tex.replace('__R__',f'{r:.2f}')
        .replace('__COV__',regpts['cov']).replace('__MRF__',regpts['mrf']).replace('__ION__',regpts['ion'])
        .replace('__A__',f'{a:.4f}').replace('__B__',f'{b:.4f}')
        .replace('__XLO__',f'{xlo:.3f}').replace('__XHI__',f'{xhi:.3f}')
        .replace('__H66__',hubcoord[(6,6)]).replace('__H88__',hubcoord[(8,8)]).replace('__H86__',hubcoord[(8,6)])
        .replace('__FREE__',freecoord))
open('paper/figs/fig_decoupling_native.tex','w',encoding='utf-8').write(tex)
print(f'wrote fig_decoupling_native.tex | molecular r={r:.3f} | free |S| {list(zip(freeL,freeS))}')
