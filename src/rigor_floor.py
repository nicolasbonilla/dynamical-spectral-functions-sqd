# -*- coding: utf-8 -*-
r"""RIGOR: the shared-floor / CP-overhead structure, on the SAME state, rigorously.
For the saturated time-evolved response seed Psi = e^{-iHt*}phi (t* = argmax spatial entanglement) on doped
2-leg Hubbard ladders (U=8), at the central rung cut, we compute on the IDENTICAL vector:
  - chi_eps  = eps-truncated Schmidt rank (min r : sum_{k>r} p_k <= eps^2), p_k = normalized squared singular
               values of the Fock-embedded matricization. By Eckart-Young this LOWER-BOUNDS any determinant
               support: a state supported on |S| determinants has matricization rank <= |S|, so |S| >= chi_eps.
               chi_eps is ALSO the minimal MPS bond dimension at that cut -> the SHARED FLOOR.
  - Seps     = eps-support size = # of largest-|coeff|^2 determinants capturing 1-eps^2 of the 2-norm
               (the determinant count actually needed to represent Psi to 2-norm eps in the fixed basis).
  - 2^{S_ent}= entropy floor (<= chi_eps).
Reports chi_eps, Seps, their ratio (the CP-overhead above the entanglement floor), and 2^{S_ent}, for several
eps, across sizes. Writes rigor_floor.json.
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
from scipy.sparse.linalg import eigsh, expm_multiply
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
import gate1_ladder as G
t0=time.time(); log=lambda *a: print(f"[{time.time()-t0:7.1f}s]",*a,flush=True)

def left_mask(L,Lr):
    half=Lr//2; left=[]
    for leg in (0,1):
        for r in range(Lr):
            site=leg*Lr+r
            if r<half: left+=[2*site,2*site+1]
    return sorted(left)

def schmidt_spectrum(vec, Su1, Sd1, L, leftorb):
    """normalized squared singular values p_k of the Fock-embedded matricization across leftorb|rest."""
    leftset=set(leftorb); nl=len(leftorb); M=2*L; nr=M-nl
    lpos={o:i for i,o in enumerate(sorted(leftorb))}
    rpos={o:i for i,o in enumerate(o for o in range(M) if o not in leftset)}
    vm=vec.reshape(len(Su1),len(Sd1)); A=np.zeros((1<<nl,1<<nr),complex)
    for a,su in enumerate(Su1):
        row=vm[a]
        for b,sd in enumerate(Sd1):
            amp=row[b]
            if abs(amp)<1e-12: continue
            li=0; ri=0
            for i in range(L):
                if (su>>i)&1:
                    o=2*i; (li:=li|(1<<lpos[o])) if o in leftset else (ri:=ri|(1<<rpos[o]))
                if (sd>>i)&1:
                    o=2*i+1; (li:=li|(1<<lpos[o])) if o in leftset else (ri:=ri|(1<<rpos[o]))
            A[li,ri]+=amp
    s=np.linalg.svd(A,compute_uv=False); p=(s[s>1e-12]**2); p/=p.sum()
    return np.sort(p)[::-1]

def eps_rank(p_desc, eps):
    """min r : sum_{k>r} p_k <= eps^2, p_desc sorted descending and normalized."""
    tail=np.cumsum(p_desc[::-1])[::-1]              # tail[r] = sum_{k>=r} p
    # sum_{k>r} = tail[r+1]; find smallest r with tail[r+1..] <= eps^2
    e2=eps*eps
    for r in range(len(p_desc)+1):
        rem = tail[r] if r<len(tail) else 0.0       # sum_{k>=r}
        if rem<=e2: return r                        # keeping k<r (i.e. r terms 0..r-1) leaves <=eps^2
    return len(p_desc)

def run(Lr, U=8.0, hole=2, tmax=8.0, nt=11, epslist=(0.1,0.05)):
    L=2*Lr; bonds=G.ladder_bonds(Lr); N=L-hole; nup=N//2; ndn=N-nup
    H,Su,iu,Sd,idd,Du,Dd=G.build_H(L,U,nup,ndn,bonds)
    w,v=eigsh(H,k=1,which='SA',v0=_v0(H.shape[0])); psi=v[:,0].reshape(Du,Dd)
    Su1,iu1=G.strings(L,nup+1); r=[];c=[];val=[]
    for a,m in enumerate(Su):
        if not (m>>0)&1: r.append(iu1[m|1]); c.append(a); val.append(1.0)
    Cd=sp.csr_matrix((val,(r,c)),shape=(len(Su1),Du)); phi=(Cd@psi).reshape(-1)
    H1,Su1b,iu1b,Sd1,idd1,Du1,Dd1=G.build_H(L,U,nup+1,ndn,bonds); nS=H1.shape[0]
    phi=phi/np.linalg.norm(phi); leftorb=left_mask(L,Lr)
    traj=expm_multiply(-1j*H1, phi, start=0.0, stop=tmax, num=nt)
    Sents=[float(-(pp*np.log2(pp)).sum()) for pp in
           (schmidt_spectrum(traj[k],Su1b,Sd1,L,leftorb) for k in range(nt))]
    kstar=int(np.argmax(Sents)); Psi=traj[kstar]/np.linalg.norm(traj[kstar])
    p=schmidt_spectrum(Psi,Su1b,Sd1,L,leftorb)          # Schmidt spectrum at t*
    Sent=float(-(p*np.log2(p)).sum())
    coeff2=np.sort((np.abs(Psi)**2))[::-1]; coeff2/=coeff2.sum()   # determinant weight spectrum
    row={'Lr':Lr,'L':L,'N':N,'nSector':nS,'S_ent':Sent,'two_pow_S':2**Sent,'kstar':kstar}
    for eps in epslist:
        chi=eps_rank(p,eps); Se=eps_rank(coeff2,eps)
        row[f'chi_eps@{eps}']=chi; row[f'Seps@{eps}']=Se; row[f'ratio@{eps}']=Se/max(chi,1)
    log(f"2x{Lr} (sector={nS}) S_ent={Sent:.2f} 2^S={2**Sent:.0f} | "
        +" | ".join(f"eps={eps}: chi_eps={row[f'chi_eps@{eps}']} <= Seps={row[f'Seps@{eps}']} "
                     f"(x{row[f'ratio@{eps}']:.1f})" for eps in epslist))
    return row

if __name__=='__main__':
    out=[run(Lr) for Lr in (3,4,5)]
    json.dump(out,open(_outpath('rigor_floor.json'),'w'),indent=1); log("WROTE rigor_floor.json")
