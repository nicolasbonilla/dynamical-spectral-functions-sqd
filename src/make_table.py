# -*- coding: utf-8 -*-
"""Emit rebuild/table_molecules.tex (scratch build dir, created on demand) from the REAL FCI-verified n19_suite.json + n19_spectral.json.
No hand-typed numbers."""
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
import json, os
d=json.load(open('data/n19_suite.json'))
sp={m['mol']:m for m in json.load(open('data/n19_spectral.json'))['mols'] if 'relL1_final' in m}
M=[m for m in d['mols'] if 'error' not in m]
def l1tex(x):
    if x<1e-13: return r'$<\!10^{-13}$'
    e=int(f'{x:.0e}'.split('e')[1]); a=x/10**e
    return rf'${a:.0f}\times10^{{{e}}}$'
COV={'H2','HF','F2','CO','N2','H2O','NH3','BeH2','H4','H6'}
MRF={'C2','O2','OH','CN','NO'}; ION={'LiF','LiH','BH','BeH'}
def grp(x): return 0 if x in COV else (1 if x in MRF else 2)
def tl(s):
    for a,b in [('2','$_2$'),('3','$_3$'),('4','$_4$'),('6','$_6$')]: s=s.replace(a,b)
    return s
for m in M: m['fd']=m['FAF_diss']/m['fafmax']; m['grp']=grp(m['mol'])
M.sort(key=lambda m:(m['grp'],-m['fd']))
GN=[r'\emph{Covalent bond-breaking}',r'\emph{Multireference / open-shell}',
    r'\emph{Ionic / weakly correlated}']
NL=r'\\'
L=[]
L.append(r'\begin{table}[t]\centering\small')
L.append(r'\caption{\textbf{The nineteen-molecule fermionic-magic suite}, computed and \textbf{verified')
L.append(r'against exact active-space FCI to $<10^{-5}$~Ha at every geometry} (cc-pVDZ). The fermionic')
L.append(r'AntiFlatness $F_1$ ($k{=}1$ Majorana antiflatness; $F_1{=}0$ iff Gaussian/mean-field, maximum')
L.append(r'$F_1^{\max}{=}2n_{\rm o}$) is listed near equilibrium and at dissociation, in a valence')
L.append(r'bond-breaking active space CAS$(n_e,n_{\rm o})$. Fermionic magic switches on for covalent')
L.append(r'bond-breaking, is already large at equilibrium for inherently multireference species')
L.append(r'($\mathrm{C_2}$, $\mathrm{O_2}$), and stays near-Gaussian for ionic and weakly-correlated')
L.append(r'bonds---so $F_1$ marks the multireference, strongly-correlated regime rather than bond length')
L.append(r'\emph{per se}. The last two columns validate the \emph{method}: at the dissociation geometry the')
L.append(r'bitstring-sampled reconstruction of the (orbital-projected) spectral function $A(\omega)$ matches')
L.append(r'the exact Lehmann result to relative-$L_1$ error, using a sampled subspace that is a fraction')
L.append(r'$|\mathcal{S}|/D$ of the full $(N{+}1)$ sector. $^{\dagger}$: near-degenerate ground state.}')
L.append(r'\label{tab:molecules}')
L.append(r'\begin{tabular}{llcccccc}\hline\hline')
L.append(r'Molecule & Bonding & CAS$(n_e,n_{\rm o})$ & $F_1^{\rm eq}$ & $F_1^{\rm diss}$ & $F_1^{\rm diss}/F_1^{\max}$ & rel-$L_1$ & $|\mathcal{S}|/D$ '+NL+r' \hline')
cur=-1
for m in M:
    if m['grp']!=cur:
        cur=m['grp']; L.append(r'\multicolumn{8}{l}{'+GN[cur]+r'}'+NL)
    dag=r'$^{\dagger}$' if m['gap_eq']<0.005 else ''
    s=sp.get(m['mol'],{}); rl=l1tex(s['relL1_final']) if s else '--'; fr=f"{s['frac_final']:.2f}" if s else '--'
    L.append(f"{tl(m['mol'])}{dag} & {m['bonding']} & ({m['nelec'][0]}+{m['nelec'][1]},{m['ncas']}) & "
             f"{m['FAF_eq']:.3f} & {m['FAF_diss']:.3f} & {m['fd']:.2f} & {rl} & {fr} {NL}")
L.append(r'\hline\hline\end{tabular}')
L.append(r'\end{table}')
OUT = os.path.join('rebuild', 'table_molecules.tex')   # scratch build directory, created on demand; not tracked in the deposit
os.makedirs(os.path.dirname(OUT), exist_ok=True)
open(OUT, 'w', encoding='utf-8').write('\n'.join(L) + '\n')
print('wrote %s  [NOTE: the AUTHORITATIVE table shipped with the paper is '
      'paper/table_molecules.tex, which carries hand-added caption text; '
      'only the data rows are regenerated here]' % OUT)
print('order:', [(m['mol'],round(m['fd'],2)) for m in M])
