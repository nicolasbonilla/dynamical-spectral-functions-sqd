# -*- coding: utf-8 -*-
r"""Dynamical SPIN structure factor S^{zz}(q,omega) of the 1D Hubbard model (PBC), exact, by the
Haydock continued-fraction (Lanczos) method -- the spin-spin response, a NEUTRAL, spin-conserving
excitation in the half-filled N-particle sector. Companion of the charge response S(q,omega).

  Sz_q = sum_j e^{i q j} (n_{j,up} - n_{j,dn})/2          (q != 0, so <Sz_q>=0)
  Szz(q,w) = -(1/pi) Im <psi0| Sz_{-q} [w-(H-E0)+i eta]^{-1} Sz_q |psi0>
The spin sector -> Heisenberg with J=4t^2/U at strong coupling: gapless 2-spinon continuum,
lower edge (des Cloizeaux-Pearson) (pi J/2)|sin q|, upper edge pi J |sin(q/2)|.
"""
import json, time, numpy as np
from scipy.sparse.linalg import eigsh
from akw_lanczos import build_H_explicit, strings, haydock
t0=time.time(); log=lambda *a: print(f"[{time.time()-t0:6.1f}s]",*a,flush=True)
def occ_matrix(L,n):
    S,_=strings(L,n); return np.array([[(m>>i)&1 for i in range(L)] for m in S],dtype=float)
L=12;U=8.0;t=1.0;eta=0.10;nw=600;nl=200; nup=nd=L//2
H,Du,Dd=build_H_explicit(L,U,nup,nd,t)
E0,V0=eigsh(H,k=1,which='SA'); E0=float(E0[0]); Psi=V0[:,0].reshape(Du,Dd)
log(f"E0={E0:.6f}")
upocc=occ_matrix(L,nup); dnocc=occ_matrix(L,nd)
wg=np.linspace(-0.1,3.5,nw); zt=wg+1j*eta
Happly=lambda x: H@x - E0*x
out={}
for n in range(1,L):
    q=2*np.pi*n/L; ph=np.exp(1j*q*np.arange(L))
    nq_up=upocc@ph; nq_dn=dnocc@ph
    seed=(((nq_up[:,None]-nq_dn[None,:])/2.0)*Psi).ravel()   # Sz_q |psi0>
    Szz=-(1.0/np.pi)*np.imag(haydock(Happly,seed,nl,zt))
    out[round(q/np.pi,5)]=Szz.tolist()
    log(f"  q/pi={2*n/L:.3f} sum={np.trapz(Szz,wg):.3f}")
d=dict(L=L,U=U,eta=eta,E0=E0,J=4*t*t/U,wg=wg.tolist(),S={f"{qq:.5f}":v for qq,v in out.items()})
json.dump(d,open('spinqw_L12.json','w'))
log("WROTE spinqw_L12.json (eta=0.10)")
# lower/upper 2-spinon edges from the data (3% threshold) for the figure
qs=np.array(sorted(out.keys())); M=np.array([out[q] for q in qs]); wgn=np.array(wg)
lines=["qpi lo hi peak"]
for q in qs:
    A=np.array(out[q]); thr=0.03*M.max(); sig=wgn[A>thr]
    lo=sig.min() if len(sig) else 0.0; hi=sig.max() if len(sig) else 0.0; pk=wgn[np.argmax(A)]
    lines.append(f"{q:.4f} {lo:.3f} {hi:.3f} {pk:.3f}")
open('paper/figs/spinqw_edges.dat','w').write("\n".join(lines)+"\n")
log("WROTE spinqw_edges.dat"); print("EDGES:\n"+"\n".join(lines))
