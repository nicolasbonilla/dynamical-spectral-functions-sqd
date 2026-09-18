# -*- coding: utf-8 -*-
r"""Multi-seed head-to-head benchmark (real ED, L=6 Hubbard): at matched subspace size |S|, which
determinant selector best reconstructs A(omega)? Q = quantum time-evolution sampling (STOCHASTIC ->
seed-averaged over 12 seeds with spread bands); K = classical polynomial Krylov; SCI = CIPSI-GF
(both deterministic). Reuses the validated setup/spectrum of headtohead.py. Writes headtohead_ms.json.
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
import json, time, numpy as np
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
from collections import Counter
import headtohead as HT
t0=time.time(); log=lambda *a: print(f"[{time.time()-t0:6.1f}s]",*a,flush=True)
SEEDS=list(range(12)); Ks=[10,20,30,40,60,80,120,160,200,260]

def run_U(L,U,eta=0.15,K=16,dt=0.5,shots=3000):
    M,dim,C,Cd,Nocc,Szocc,H=HT.setup(L,U)
    gi=np.where((Nocc==L)&(Szocc==0))[0]
    w,v=np.linalg.eigh(H[gi][:,gi].toarray()); E0=w[0]; psi0=np.zeros(dim,complex); psi0[gi]=v[:,0]
    phi=Cd[0]@psi0
    si=np.where((Nocc==L+1)&(Szocc==1))[0]; nS=len(si)
    H7=H[si][:,si].toarray(); En,Vn=np.linalg.eigh(H7); phi7=phi[si]; wnorm=float(np.vdot(phi7,phi7).real)
    coef=Vn.conj().T@phi7; poles=En-E0; grid=np.linspace(poles.min()-1,poles.max()+1,600)
    A_ex=HT.spectrum(poles,np.abs(coef)**2,grid,eta); norm=np.trapezoid(A_ex,grid)
    def relL1(S):
        HS=H7[np.ix_(S,S)]; Em,Um=np.linalg.eigh(HS); a=Um.conj().T@phi7[S]
        return float(np.trapezoid(np.abs(HT.spectrum(Em-E0,np.abs(a)**2,grid,eta)-A_ex),grid)/norm)
    kk=[k for k in Ks if k<=nS]
    # --- deterministic Krylov (K) ---
    Qb=[]; q=phi7/np.linalg.norm(phi7)
    for _ in range(min(60,nS)):
        for u in Qb: q=q-np.vdot(u,q)*u
        n=np.linalg.norm(q)
        if n<1e-8: break
        q=q/n; Qb.append(q.copy()); q=H7@Qb[-1]
    wcfg=(np.abs(np.array(Qb))**2).sum(0); order_K=list(np.argsort(wcfg)[::-1])
    curveK=[relL1(np.array(sorted(order_K[:k]))) for k in kk]
    # --- deterministic CIPSI-GF (SCI) ---
    Eref=float((phi7.conj()@(H7@phi7)).real/wnorm); Hdiag=np.diag(H7)
    Sset=[int(np.argmax(np.abs(phi7)**2))]; order_SCI=[Sset[0]]
    while len(order_SCI)<min(260,nS):
        rep=np.zeros(nS,complex); rep[Sset]=phi7[Sset]; Hrep=H7@rep
        score=np.abs(Hrep)**2/np.maximum(np.abs(Eref-Hdiag),0.1)
        for i in order_SCI: score[i]=-1
        add=int(np.argmax(score))
        if score[add]<=0: break
        order_SCI.append(add); Sset.append(add)
    curveSCI=[relL1(np.array(sorted(order_SCI[:k]))) for k in kk if k<=len(order_SCI)]
    # --- quantum Q: seed-averaged ---
    Qcurves=[]
    for sd in SEEDS:
        rng=np.random.default_rng(sd); freq=Counter()
        for k in range(K+1):
            vk=Vn@(np.exp(-1j*En*k*dt)*coef); pr=np.abs(vk)**2; pr/=pr.sum()
            for dd in rng.choice(nS,size=shots,p=pr): freq[int(dd)]+=1
        order_Q=[i for i,_ in freq.most_common()]
        Qcurves.append([relL1(np.array(sorted(order_Q[:k]))) for k in kk])
    Qc=np.array(Qcurves); Qm=Qc.mean(0).tolist(); Qs=Qc.std(0).tolist()
    log(f"U={U}: sector={nS}  |S|=200 -> Q={Qm[kk.index(200)]:.3f}+/-{Qs[kk.index(200)]:.3f} "
        f"K={curveK[kk.index(200)]:.3f} SCI={curveSCI[kk.index(200)]:.3f}")
    return {'U':U,'nSector':nS,'Ks':kk,'Q_mean':Qm,'Q_std':Qs,'K':curveK,'SCI':curveSCI}

if __name__=='__main__':
    out={f"U{U}":run_U(6,U) for U in (4.0,8.0,12.0)}
    json.dump(out,open(_outpath('headtohead_ms.json'),'w'),indent=1)
    log("WROTE headtohead_ms.json")
