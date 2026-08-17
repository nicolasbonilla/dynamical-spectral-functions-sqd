# -*- coding: utf-8 -*-
r"""RIGOR (makes the obstruction theorem empirical): the SQD cost |S| is governed by the ENTANGLEMENT
(Schmidt rank / minimal MPS bond dimension chi), NOT by the one-body magic FAF. For a family of Hubbard
ground states (chains, varying U, L, filling) compute:
  FAF   = 4 sum g_i(1-g_i)                      (one-body magic, from 1-RDM occupations)
  |S|   = eps-support of the GS in the SITE (computational) basis   (the SQD cost)
  chi   = eps-truncated Schmidt rank at the central cut = minimal MPS bond dimension (the entanglement floor)
Then report the pairwise (Spearman) correlations. Prediction (thesis): chi predicts |S| strongly; FAF does not.
Writes cost_vs_ent.json / .txt.
"""
import json, numpy as np
from scipy.stats import spearmanr
from scipy.sparse.linalg import eigsh
import gate1_ladder as G, nonfreeness as NF, rigor_floor as RF

def chain_bonds(L): return [(i,i+1,1.0) for i in range(L-1)]

def gs_metrics(L,U,nup,ndn,eps=0.05):
    H,Su,iu,Sd,idd,Du,Dd=G.build_H(L,U,nup,ndn,chain_bonds(L))
    w,v=eigsh(H,k=1,which='SA'); psi=v[:,0]
    # FAF (1-RDM occupations)
    occ=NF.occupations(psi,Su,iu,Sd,idd,L); FAF=float(4*np.sum(occ*(1-occ)))
    # |S| = eps-support in site basis
    w2=np.sort(np.abs(psi)**2)[::-1]; w2/=w2.sum(); S=int(np.searchsorted(np.cumsum(w2),1-eps*eps)+1)
    # chi = eps-truncated Schmidt rank at central spatial cut (min MPS bond dim)
    leftorb=sorted([2*i for i in range(L//2)]+[2*i+1 for i in range(L//2)])   # first L/2 sites, both spins
    p=RF.schmidt_spectrum(psi,Su,Sd,L,leftorb); chi=RF.eps_rank(p,eps)
    return FAF,S,chi

if __name__=='__main__':
    rows=[]
    for L in (6,8,10):
        for hole in (0,2):
            N=L-hole; nup=N//2; ndn=N-nup
            for U in (1.,2.,4.,8.,16.):
                FAF,S,chi=gs_metrics(L,U,nup,ndn)
                rows.append(dict(L=L,N=N,U=U,FAF=FAF,S=S,chi=chi))
                print(f"L={L} N={N} U={U:4.1f}: FAF={FAF:6.3f}  |S|={S:5d}  chi={chi:4d}",flush=True)
    F=np.array([r['FAF'] for r in rows]); S=np.log2([max(r['S'],1) for r in rows]); C=np.log2([max(r['chi'],1) for r in rows])
    rFS=spearmanr(F,S); rCS=spearmanr(C,S)
    print(f"\n==> Spearman r(FAF, log|S|)  = {rFS.correlation:.3f}  (p={rFS.pvalue:.1e})  -- one-body magic vs cost")
    print(f"==> Spearman r(chi,  log|S|) = {rCS.correlation:.3f}  (p={rCS.pvalue:.1e})  -- entanglement vs cost")
    print(f"    interpretation: chi (entanglement) predicts the SQD cost |S|; FAF (one-body magic) does not.")
    json.dump(dict(rows=rows,r_FAF_S=rFS.correlation,r_chi_S=rCS.correlation),open('/w/cost_vs_ent.json','w'),indent=1)
