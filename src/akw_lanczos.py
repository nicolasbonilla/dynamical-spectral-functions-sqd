# -*- coding: utf-8 -*-
r"""Momentum-resolved spectral function A(k,omega) of the 1D Hubbard model (PBC)
via the Haydock continued-fraction (Lanczos) method in a SECTOR basis
(up-string (x) dn-string), so it scales to L=10-12 without full diagonalization.

  G^+_k(z) = <psi0| c_{k,up} [z-(H-E0)]^{-1} c^dag_{k,up} |psi0> ,  z=omega+i*eta
  G^-_k(z) = <psi0| c^dag_{k,up} [z-(E0-H)]^{-1} c_{k,up} |psi0>
  A(k,w) = -(1/pi) Im[G^+_k(w) + G^-_k(w)]
Ground state by sparse Lanczos; each G by ~n_lanc Haydock steps (sparse matvecs).
Validated against exact sector diagonalization at L=6.
"""
# --- DEPOSIT PATHS (repaired 2026-09-18) -----------------------------------
# This script used to hard-code its output under '/w/'.  /w was the working
# directory of the Docker container the published runs were made in; outside that
# container the documented pipeline wrote nothing a reader could find, and data/
# was in fact repopulated BY HAND.  That made `make data` and the README recipe
# untrue.  Repaired: every path is now an ARGUMENT with a default RELATIVE TO THIS
# REPOSITORY, so a clean clone reproduces into its own tree.
#     read   <repo>/data/<name>      override with  --in  PATH
#     write  <repo>/data/<name>      override with  --out PATH
#     write  <repo>/build/<name>     for by-products that are NOT part of the deposit
# Paths only -- no physics and no computational default was changed here.
import os as _os, sys as _sys
_REPO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))

def _flag(name):
    """Value of a `--name VALUE` (or `-n VALUE`) command-line flag, else None."""
    a = _sys.argv[1:]
    for f in ('--' + name, '-' + name[0]):
        if f in a and a.index(f) + 1 < len(a):
            return a[a.index(f) + 1]
    return None

def _argv_positional():
    """argv[1:] with the --in/--out flags and their values removed."""
    a, keep, i = _sys.argv[1:], [], 0
    while i < len(a):
        if a[i] in ('--out', '-o', '--in', '-i'):
            i += 2
            continue
        keep.append(a[i]); i += 1
    return keep

def _outpath(name, sub='data'):
    """Absolute path to write `name` to: --out if given, else <repo>/<sub>/<name>."""
    p = _os.path.abspath(_flag('out') or _os.path.join(_REPO, sub, name))
    _os.makedirs(_os.path.dirname(p), exist_ok=True)
    return p

def _inpath(name, sub='data'):
    """Absolute path to read `name` from: --in if given, else <repo>/<sub>/<name>."""
    return _os.path.abspath(_flag('in') or _os.path.join(_REPO, sub, name))
# ---------------------------------------------------------------------------
import json, time, numpy as np, scipy.sparse as sp
# --- NUMPY_TRAPEZOID_BRIDGE ------------------------------------------------
# numpy 2.0 ADDED np.trapezoid and REMOVED np.trapz.  Files in this repository use
# both names, so without this bridge no single numpy version runs the whole deposit:
# numpy 1.x breaks the files that call trapezoid, numpy 2.x breaks the files that call
# trapz (this guardian included).  requirements.txt asks for numpy>=1.24; with the
# bridge that is true again.
if not hasattr(np, "trapezoid"):
    np.trapezoid = np.trapz          # numpy < 2.0
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid          # numpy >= 2.0
# ---------------------------------------------------------------------------
from itertools import combinations
from scipy.sparse.linalg import eigsh, LinearOperator
# --- DETERMINISTIC ARPACK START VECTOR (added 2026-09-18, second pass) ------
# eigsh() with no v0= lets ARPACK draw its own start vector from an UNSEEDED
# generator, so E0 converges to a slightly different point on every run (the last
# few digits move) and every quantity derived from it moves with it.  Measured, not
# hypothetical: three consecutive calls on the same matrix gave
# -2.0481308860914536 / ...504 / ...522, and in the leakage certificate a small
# subspace selection moved by 3.2%.  The eigenpair is the same to ARPACK's
# tolerance -- no physics changes -- but a deposit must be bit-reproducible.
# Deliberately NOT a numpy global seed: this touches only the ARPACK start vector.
# The seed 20260918 is the one used by scaling_lanczos.py and the certificate suite.
def _v0(n):
    """Fixed, dimension-dependent ARPACK start vector (never orthogonal to the GS)."""
    return np.random.default_rng(20260918).standard_normal(n)
# ---------------------------------------------------------------------------
t0=time.time(); log=lambda *a: print(f"[{time.time()-t0:6.1f}s]",*a,flush=True)

def strings(L,n):
    S=[0]*0; S=[]
    for c in combinations(range(L),n):
        m=0
        for b in c: m|=(1<<b)
        S.append(m)
    S.sort(); idx={m:i for i,m in enumerate(S)}
    return S,idx

def hop(L,n,t=1.0):
    """single-spin hopping matrix -t sum_<ij>(c^dag_i c_j + h.c.), PBC, in n-string basis."""
    S,idx=strings(L,n); D=len(S); rows=[];cols=[];val=[]
    for a,m in enumerate(S):
        for i in range(L):
            j=(i+1)%L
            for (p,q) in ((i,j),(j,i)):   # c^dag_p c_q
                if (m>>q)&1 and not (m>>p)&1:
                    m2=(m & ~(1<<q))|(1<<p)
                    lo,hi=min(p,q),max(p,q)
                    mask=m & (((1<<hi)-1) ^ ((1<<(lo+1))-1))  # occupied strictly between q and p
                    sign=-1.0 if (bin(mask).count('1')&1) else 1.0
                    rows.append(idx[m2]); cols.append(a); val.append(-t*sign)
    return sp.csr_matrix((val,(rows,cols)),shape=(D,D)),S,idx

def cdag_map(L,n):
    """list over site j of sparse matrix c^dag_j : n-string -> (n+1)-string."""
    Sa,ia=strings(L,n); Sb,ib=strings(L,n+1)
    ops=[]
    for j in range(L):
        rows=[];cols=[];val=[]
        for a,m in enumerate(Sa):
            if not (m>>j)&1:
                m2=m|(1<<j)
                sign=-1.0 if (bin(m & ((1<<j)-1)).count('1')&1) else 1.0
                rows.append(ib[m2]); cols.append(a); val.append(sign)
        ops.append(sp.csr_matrix((val,(rows,cols)),shape=(len(Sb),len(Sa))))
    return ops  # c_j = ops[j].T

def sector_H(L,U,nup,ndn,t=1.0):
    Tu,Su,iu=hop(L,nup,t); Td,Sd,idd=hop(L,ndn,t)
    Du=len(Su); Dd=len(Sd)
    # diagonal U: number of doubly-occupied sites (fast matrix product, no big broadcast)
    upocc=np.array([[ (m>>i)&1 for i in range(L)] for m in Su],dtype=float)  # (Du,L)
    dnocc=np.array([[ (m>>i)&1 for i in range(L)] for m in Sd],dtype=float)  # (Dd,L)
    diagU=U*(upocc@dnocc.T).ravel()                                          # (Du*Dd,)
    def matvec(x):
        X=x.reshape(Du,Dd)
        Y=Tu@X + (Td@X.T).T           # X @ Td (Td symmetric): fast sparse-times-dense both terms
        return Y.ravel()+diagU*x
    H=LinearOperator((Du*Dd,Du*Dd),matvec=matvec,dtype=float)
    return H,Su,iu,Sd,idd,Du,Dd

def build_H_explicit(L,U,nup,nd,t=1.0):
    """explicit sparse H = Tu(x)I + I(x)Td + diag(U*doublons), for a fast, reliable ground state."""
    Tu,Su,iu=hop(L,nup,t); Td,Sd,idd=hop(L,nd,t)
    Du=len(Su); Dd=len(Sd)
    upocc=np.array([[ (m>>i)&1 for i in range(L)] for m in Su],dtype=float)
    dnocc=np.array([[ (m>>i)&1 for i in range(L)] for m in Sd],dtype=float)
    diagU=U*(upocc@dnocc.T).ravel()
    H=(sp.kron(Tu,sp.identity(Dd),format='csr')
       +sp.kron(sp.identity(Du),Td,format='csr')
       +sp.diags(diagU))
    return H.tocsr(),Du,Dd

def haydock(H,v0,nl,z):
    """continued fraction <v0|(z-H)^{-1}|v0> via Lanczos (v0 need not be normalized)."""
    nrm2=np.vdot(v0,v0).real
    if nrm2<1e-14: return np.zeros_like(z)
    v=v0/np.sqrt(nrm2); vp=np.zeros_like(v); b=0.0; a=[];bs=[]
    for _ in range(nl):
        w=H(v); an=np.vdot(v,w).real; a.append(an)
        w=w-an*v-b*vp
        bn=np.sqrt(np.vdot(w,w).real); bs.append(bn)
        if bn<1e-10: break
        vp=v; v=w/bn; b=bn
    a=np.array(a); bs=np.array(bs)
    G=z-a[-1]
    for m in range(len(a)-2,-1,-1):
        G=(z-a[m]) - bs[m]**2 / G
    G=nrm2/G
    return G

def run(L=6,U=8.0,t=1.0,eta=0.15,nl=120,nw=600,ntw=1,wwin=None):
    nup=nd=L//2
    Hexpl,Du,Dd=build_H_explicit(L,U,nup,nd,t)
    # ground state (explicit sparse -> fast, reliable ARPACK)
    E0,V0=eigsh(Hexpl,k=1,which='SA',v0=_v0(Hexpl.shape[0])); E0=float(E0[0]); psi0=V0[:,0]
    log(f"L={L} U={U}: gs dim={Du*Dd}  E0={E0:.6f}")
    Psi=psi0.reshape(Du,Dd)
    Ha,_,_,_,_,Dua,_=sector_H(L,U,nup+1,nd,t)
    Hr,_,_,_,_,Dur,_=sector_H(L,U,nup-1,nd,t)
    cdU=cdag_map(L,nup)                 # c^dag_j : nup -> nup+1  (Dua x Du)
    cU =cdag_map(L,nup-1)               # c^dag_j : nup-1 -> nup  ; c_j = cU[j].T : nup -> nup-1
    if wwin is None: wwin=(-1.0*U-2,1.0*U+2)
    mu=U/2.0                                  # half-filling Fermi level (particle-hole symmetric)
    wg=np.linspace(wwin[0],wwin[1],nw)         # omega measured from E_F
    zt=(wg+mu)+1j*eta                          # evaluate G at raw energy = omega + mu
    ks=[2*np.pi*n/L for n in range(L)]
    out={}
    for n,k in enumerate(ks):
        ph=np.exp(1j*k*np.arange(L))/np.sqrt(L)
        # add: |i> = c^dag_{k,up}|psi0> in (nup+1,nd)
        add=np.zeros((Dua,Dd),dtype=complex)
        for j in range(L): add+=ph[j]*(cdU[j]@Psi)
        Gp=haydock(lambda x:_apply(Ha,x)-E0*x, add.ravel(), nl, zt)   # add: poles at (E_m^{N+1}-E0)-mu
        # remove: |r> = c_{k,up}|psi0> in (nup-1,nd)
        rem=np.zeros((Dur,Dd),dtype=complex)
        for j in range(L): rem+=np.conj(ph[j])*(cU[j].T@Psi)
        Gm=haydock(lambda x:E0*x-_apply(Hr,x), rem.ravel(), nl, zt)   # remove: poles at (E0-E_m^{N-1})-mu
        A=-(1.0/np.pi)*np.imag(Gp+Gm)
        out[round((k/np.pi)%2,5)]=A
        log(f"  k/pi={2*n/L:.3f}  sum A dw={np.trapezoid(A,wg):.3f}")
    return dict(L=L,U=U,eta=eta,E0=E0,wg=wg.tolist(),
                A={f"{kk:.5f}":v.tolist() for kk,v in out.items()})

def _apply(Hlin,x):
    return Hlin.matvec(x)

if __name__=="__main__":
    import sys
    _a=_argv_positional(); L=int(_a[0]) if _a else 12
    d=run(L=L,U=8.0,eta=0.18,nl=150,nw=640,wwin=(-9,9))
    json.dump(d,open(_outpath(f'akw_lanczos_L{L}.json'),'w'))
    log(f"WROTE akw_lanczos_L{L}.json ({len(d['A'])} native k-points)")
