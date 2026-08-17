# -*- coding: utf-8 -*-
r"""Dynamical structure factor S(q,omega) of the 1D Hubbard model (PBC), exact, by the
Haydock continued-fraction (Lanczos) method -- the density-density response, a NEUTRAL
(particle-conserving) excitation, computed in the half-filled N-particle sector.

  n_q = sum_j e^{i q j} (n_{j,up}+n_{j,dn})              (q != 0, so <n_q>=0)
  S(q,w) = -(1/pi) Im <psi0| n_{-q} [w-(H-E0)+i eta]^{-1} n_q |psi0>
         = sum_m |<m|n_q|psi0>|^2 delta(w-(E_m-E0))
Same sector engine / Haydock as akw_lanczos (exact spectral functions).
"""
import json, time, numpy as np
from scipy.sparse.linalg import eigsh
from akw_lanczos import build_H_explicit, strings, haydock
import sys
t0=time.time(); log=lambda *a: print(f"[{time.time()-t0:6.1f}s]",*a,flush=True)

def occ_matrix(L,n):
    S,_=strings(L,n)
    return np.array([[ (m>>i)&1 for i in range(L)] for m in S],dtype=float)   # (D,L)

def run(L=12,U=8.0,t=1.0,eta=0.20,nw=600,wwin=(-0.5,13.5),nl=200):
    nup=nd=L//2
    H,Du,Dd=build_H_explicit(L,U,nup,nd,t)
    E0,V0=eigsh(H,k=1,which='SA'); E0=float(E0[0]); Psi=V0[:,0].reshape(Du,Dd)
    log(f"L={L} gs dim={Du*Dd} E0={E0:.4f}")
    upocc=occ_matrix(L,nup); dnocc=occ_matrix(L,nd)
    wg=np.linspace(wwin[0],wwin[1],nw); zt=wg+1j*eta
    Happly=lambda x: H@x - E0*x                 # M = H - E0 ; poles at w = E_m-E0
    out={}
    for n in range(1,L):                        # q = 2 pi n / L, skip q=0
        q=2*np.pi*n/L
        ph=np.exp(1j*q*np.arange(L))
        nq_up=upocc@ph; nq_dn=dnocc@ph
        seed=((nq_up[:,None]+nq_dn[None,:])*Psi).ravel()
        G=haydock(Happly, seed, nl, zt)
        S=-(1.0/np.pi)*np.imag(G)
        out[round(q/np.pi,5)]=S.tolist()
        log(f"  q/pi={2*n/L:.3f}  sum S dw={np.trapezoid(S,wg):.3f}")
    return dict(L=L,U=U,eta=eta,E0=E0,wg=wg.tolist(),
                S={f"{qq:.5f}":v for qq,v in out.items()})

if __name__=="__main__":
    L=int(sys.argv[1]) if len(sys.argv)>1 else 12
    d=run(L=L)
    json.dump(d,open(f'/w/sqw_L{L}.json','w'))
    log(f"WROTE sqw_L{L}.json ({len(d['S'])} q-points)")
