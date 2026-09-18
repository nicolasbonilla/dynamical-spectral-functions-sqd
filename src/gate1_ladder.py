# -*- coding: utf-8 -*-
r"""GATE 1 (Plan B make-or-break): does RESPONSE-targeted determinant selection give a smaller /
better-scaling subspace for the spectral function A(omega) than ENERGY-biased selection, on a doped
2-leg Hubbard ladder? Tests the LR-SCI insight (arXiv:2510.02949) in the sampling setting.
Sector-basis ED (up-string (x) down-string). For a frontier orbital seed we compare, at matched
spectral accuracy rel-L1<thr:
  (E)  energy-biased  : CIPSI grown from the HF determinant by Epstein-Nesbet energy score (ground-
                        state-relevant subspace), then A(omega) reconstructed in it.
  (R)  response-target: configs ranked by the time-evolved-seed weight of phi=c^dag_p|psi0> (the
                        spectral/response-relevant subspace).
Reports |S|/D at rel-L1<thr for each, per ladder size -> the curve that decides Plan B.
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
import json, time, itertools, numpy as np
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
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
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

def ladder_bonds(Lr, tleg=1.0, trung=1.0):
    """2-leg ladder: site index = leg*Lr + r, leg in {0,1}, r in 0..Lr-1 (open legs, rungs)."""
    B=[]
    for leg in (0,1):
        for r in range(Lr-1): B.append((leg*Lr+r, leg*Lr+r+1, tleg))
    for r in range(Lr): B.append((r, Lr+r, trung))
    return B

def strings(L,n):
    S=[]
    for c in itertools.combinations(range(L),n):
        m=0
        for b in c: m|=(1<<b)
        S.append(m)
    S.sort(); return S,{m:i for i,m in enumerate(S)}

def hop_bonds(L,n,bonds):
    """single-spin hopping -t sum_bonds(c^dag_i c_j + h.c.) in the n-string basis, for a bond list."""
    S,idx=strings(L,n); D=len(S); r=[];c=[];v=[]
    for a,m in enumerate(S):
        for (i,j,t) in bonds:
            for (p,q) in ((i,j),(j,i)):
                if (m>>q)&1 and not (m>>p)&1:
                    m2=(m&~(1<<q))|(1<<p)
                    lo,hi=min(p,q),max(p,q)
                    mask=m & (((1<<hi)-1) ^ ((1<<(lo+1))-1))
                    sg=-1.0 if (bin(mask).count('1')&1) else 1.0
                    r.append(idx[m2]); c.append(a); v.append(-t*sg)
    return sp.csr_matrix((v,(r,c)),shape=(D,D)),S,idx

def build_H(L,U,nup,ndn,bonds):
    Tu,Su,iu=hop_bonds(L,nup,bonds); Td,Sd,idd=hop_bonds(L,ndn,bonds)
    Du,Dd=len(Su),len(Sd)
    upo=np.array([[(m>>i)&1 for i in range(L)] for m in Su],float)
    dno=np.array([[(m>>i)&1 for i in range(L)] for m in Sd],float)
    diagU=U*(upo@dno.T).ravel()
    H=(sp.kron(Tu,sp.identity(Dd),format='csr')+sp.kron(sp.identity(Du),Td,format='csr')+sp.diags(diagU)).tocsr()
    return H,Su,iu,Sd,idd,Du,Dd

def run(Lr, U=8.0, hole=2, eta=0.15, thr=0.05, K=20, dt=0.5, p_site=0):
    L=2*Lr; bonds=ladder_bonds(Lr)
    N=L-hole                      # doped: half filling (=L) minus 'hole' electrons
    nup=N//2; ndn=N-nup
    H,Su,iu,Sd,idd,Du,Dd=build_H(L,U,nup,ndn,bonds)
    w,v=eigsh(H,k=1,which='SA',v0=_v0(H.shape[0])); E0=float(w[0]); psi=v[:,0].reshape(Du,Dd)
    # response seed phi = c^dag_{p,up} |psi0>  ->  (nup+1, ndn) sector
    cdU=[]; Su1,iu1=strings(L,nup+1)
    # build c^dag_p (up) map nup->nup+1
    r=[];c=[];val=[]
    for a,m in enumerate(Su):
        if not (m>>p_site)&1:
            m2=m|(1<<p_site); sg=-1.0 if (bin(m&((1<<p_site)-1)).count('1')&1) else 1.0
            r.append(iu1[m2]); c.append(a); val.append(sg)
    Cd=sp.csr_matrix((val,(r,c)),shape=(len(Su1),Du))
    phi_mat=(Cd@psi)                                   # (Du1, Dd)
    H1,Su1b,iu1b,Sd1,idd1,Du1,Dd1=build_H(L,U,nup+1,ndn,bonds)
    phi=phi_mat.reshape(-1); nS=H1.shape[0]; wn=float(phi@phi)
    Harr=H1.toarray(); Em,Um=np.linalg.eigh(Harr); coef=Um.conj().T@phi
    poles=Em-E0; grid=np.linspace(poles.min()-1,poles.max()+1,600)
    def spec(pw,ww):
        A=np.zeros_like(grid)
        for a,b in zip(pw,ww): A+=b*(eta/np.pi)/((grid-a)**2+eta**2)
        return A
    A_ex=spec(poles,np.abs(coef)**2); nrm=np.trapezoid(A_ex,grid)
    def relL1(idxs):
        Ss=np.sort(np.array(idxs)); HS=Harr[np.ix_(Ss,Ss)]
        e,u=np.linalg.eigh(HS); a=u.conj().T@phi[Ss]
        return float(np.trapezoid(np.abs(spec(e-E0,np.abs(a)**2)-A_ex),grid)/nrm)
    def frac_for(order):
        for k in range(2,nS+1,max(1,nS//80)):
            if relL1(order[:k])<thr: return k/nS, k
        return 1.0,nS
    # (R) response-targeted: rank by time-evolved-seed weight
    wc=np.zeros(nS)
    for k in range(K+1):
        vk=Um@(np.exp(-1j*Em*k*dt)*coef); wc+=np.abs(vk)**2
    orderR=list(np.argsort(wc)[::-1]); fracR,kR=frac_for(orderR)
    # (E) energy-biased: CIPSI from HF determinant (of the N+1 sector), Epstein-Nesbet score
    hf=int(np.argmax([Harr[i,i] for i in range(nS)]))  # placeholder; use lowest-diagonal as HF-like
    hf=int(np.argmin(np.diag(Harr)))
    Sset=[hf]; Hdiag=np.diag(Harr); Eref=Harr[hf,hf]
    while len(Sset)<nS:
        rep=np.zeros(nS);
        for i in Sset: rep[i]=1.0
        Hrep=Harr@rep; sc=(Hrep**2)/np.maximum(np.abs(Eref-Hdiag),0.1)
        for i in Sset: sc[i]=-1
        add=int(np.argmax(sc));
        if sc[add]<=0: break
        Sset.append(add)
    orderE=Sset+[i for i in range(nS) if i not in set(Sset)]; fracE,kE=frac_for(orderE)
    log(f"ladder 2x{Lr} (L={L},N={N},hole={hole}) sector(N+1)={nS}: "
        f"RESPONSE |S|={kR} frac={fracR:.3f} | ENERGY |S|={kE} frac={fracE:.3f} | ratio E/R={kE/max(kR,1):.2f}")
    return dict(Lr=Lr,L=L,N=N,nSector=nS,fracR=fracR,kR=kR,fracE=fracE,kE=kE)

if __name__=='__main__':
    out=[]
    for Lr in (3,4,5):
        out.append(run(Lr))
    json.dump(out,open(_outpath('gate1_ladder.json'),'w'),indent=1)
    log("WROTE gate1_ladder.json")
