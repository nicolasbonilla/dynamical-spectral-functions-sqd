# -*- coding: utf-8 -*-
r"""G1 (landmark): put the resource story in the geometry where the advantage is claimed.
Chain vs 2-leg LADDER, matched site count, half filling, same U. For each we compute the SAME
quantities the paper uses:
  F1   = 4 sum g(1-g)      one-body fermionic magic (1-RDM, orbital-rotation invariant)
  |S|  = eps-support of the GS in the site (computational) basis   (the SQD cost)
  chi  = eps-truncated Schmidt rank at the central SPATIAL cut = minimal MPS bond dimension
  frac = (N+1)-sector fraction that a response-targeted subspace needs for spectral rel-L1<thr
Thesis, demonstrated in 2D: at MATCHED F1 the ladder (2-bond perpendicular cut, area-law boundary 2)
needs a HIGHER chi and a LARGER |S|/frac than the chain (1-bond cut) -- the cost is set by chi, not F1,
exactly where tensor networks lose their edge. Reuses gate1_ladder / nonfreeness / rigor_floor.
"""
# --- DEPOSIT PATHS (repaired 2026-09-18, second pass) -----------------------
# Eighteen scripts were repaired earlier today because they hard-coded their output
# under '/w/', the working directory of the Docker container the published runs were
# made in.  THIS FILE WAS NOT AMONG THEM, and it was broken in a quieter way: it wrote
# to the RELATIVE path 'data/...', which lands in whatever directory the reader happens
# to be standing in, and raises FileNotFoundError from anywhere except the repository
# root.  The Data Availability Statement claims that every deposited computation script
# writes into the repository's own data/ directory; that sentence was FALSE for this
# file until now.  Repaired exactly like the other eighteen: the default is computed
# from THIS FILE's location, and --out overrides it.
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


def _outpath(name, sub='data'):
    """Absolute path to write `name` to: --out if given, else <repo>/<sub>/<name>."""
    p = _os.path.abspath(_flag('out') or _os.path.join(_REPO, sub, name))
    _os.makedirs(_os.path.dirname(p), exist_ok=True)
    return p


def _repopath(name, sub='data'):
    """Absolute path <repo>/<sub>/<name>.  Never overridden: for a script's SECOND
    output, which --out (a single flag) cannot address unambiguously."""
    p = _os.path.abspath(_os.path.join(_REPO, sub, name))
    _os.makedirs(_os.path.dirname(p), exist_ok=True)
    return p
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
import gate1_ladder as G, nonfreeness as NF, rigor_floor as RF
t0=time.time(); log=lambda *a: print(f"[{time.time()-t0:6.1f}s]",*a,flush=True)

def chain_bonds(L):  return [(i,i+1,1.0) for i in range(L-1)]                       # open chain
def chain_left(L):   return sorted([2*i for i in range(L//2)]+[2*i+1 for i in range(L//2)])
def ladder_left(Lr):                                                               # perpendicular cut: left half-columns of BOTH legs
    sites=[r for r in range(Lr//2)]+[Lr+r for r in range(Lr//2)]
    return sorted([2*s for s in sites]+[2*s+1 for s in sites])

def spec(grid,pw,ww,eta):
    A=np.zeros_like(grid)
    for a,b in zip(pw,ww): A+=b*(eta/np.pi)/((grid-a)**2+eta**2)
    return A

def metrics(L,U,bonds,leftorb,hole=0,eps=0.05,eta=0.15,thr=0.05,K=20,dt=0.5,p_site=0,do_frac=True):
    N=L-hole; nup=N//2; ndn=N-nup
    H,Su,iu,Sd,idd,Du,Dd=G.build_H(L,U,nup,ndn,bonds)
    w,v=eigsh(H,k=1,which='SA',v0=_v0(H.shape[0])); E0=float(w[0]); psi=v[:,0]
    occ=NF.occupations(psi,Su,iu,Sd,idd,L); F1=float(4*np.sum(occ*(1-occ)))
    w2=np.sort(np.abs(psi)**2)[::-1]; w2/=w2.sum(); S=int(np.searchsorted(np.cumsum(w2),1-eps*eps)+1)
    p=RF.schmidt_spectrum(psi,Su,Sd,L,leftorb); chi=RF.eps_rank(p,eps)
    D=Du*Dd
    res=dict(L=L,U=U,N=N,F1=round(F1,4),Sdet=int(S),D=int(D),chi=int(chi))
    if do_frac:
        # (N+1,Sz+1) response-targeted fraction for the frontier-orbital seed phi=c^dag_p|psi0>
        psi_mat=psi.reshape(Du,Dd); Su1,iu1=G.strings(L,nup+1)
        r=[];c=[];val=[]
        for a,m in enumerate(Su):
            if not (m>>p_site)&1:
                m2=m|(1<<p_site); sg=-1.0 if (bin(m&((1<<p_site)-1)).count('1')&1) else 1.0
                r.append(iu1[m2]); c.append(a); val.append(sg)
        import scipy.sparse as sp
        Cd=sp.csr_matrix((val,(r,c)),shape=(len(Su1),Du)); phi=(Cd@psi_mat).reshape(-1)
        H1,_,_,_,_,Du1,Dd1=G.build_H(L,U,nup+1,ndn,bonds); nS=H1.shape[0]
        Harr=H1.toarray(); Em,Um=np.linalg.eigh(Harr); coef=Um.conj().T@phi
        poles=Em-E0; grid=np.linspace(poles.min()-1,poles.max()+1,600)
        A_ex=spec(grid,poles,np.abs(coef)**2,eta); nrm=np.trapz(A_ex,grid)
        def relL1(idxs):
            Ss=np.sort(np.array(idxs)); e,u=np.linalg.eigh(Harr[np.ix_(Ss,Ss)]); a=u.conj().T@phi[Ss]
            return float(np.trapz(np.abs(spec(grid,e-E0,np.abs(a)**2,eta)-A_ex),grid)/nrm)
        wc=np.zeros(nS)
        for k in range(K+1):
            vk=Um@(np.exp(-1j*Em*k*dt)*coef); wc+=np.abs(vk)**2
        order=list(np.argsort(wc)[::-1])
        frac=1.0; kk=nS
        for k in range(2,nS+1,max(1,nS//80)):
            if relL1(order[:k])<thr: frac=k/nS; kk=k; break
        res.update(nSector=int(nS),kR=int(kk),frac=round(float(frac),4))
    return res

if __name__=='__main__':
    U=8.0; out={'U':U,'eps':0.05,'note':'chain(OBC) vs 2-leg ladder, half filling, U=8; F1,|S|,chi,(N+1) response frac','pairs':[]}
    cases=[('chain-L6',   lambda: metrics(6, U, chain_bonds(6),  chain_left(6))),
           ('ladder-2x3', lambda: metrics(6, U, G.ladder_bonds(3), ladder_left(3))),
           ('chain-L8',   lambda: metrics(8, U, chain_bonds(8),  chain_left(8))),
           ('ladder-2x4', lambda: metrics(8, U, G.ladder_bonds(4), ladder_left(4))),
           ('chain-L10',  lambda: metrics(10,U, chain_bonds(10), chain_left(10), do_frac=False)),
           ('ladder-2x5', lambda: metrics(10,U, G.ladder_bonds(5), ladder_left(5), do_frac=False))]
    res={}
    for name,fn in cases:
        r=fn(); res[name]=r
        log(f"{name:12s} L={r['L']} F1={r['F1']:.3f} |S|={r['Sdet']:5d} chi={r['chi']:3d} "
            f"frac={r.get('frac','--')} (N+1 |S|={r.get('kR','--')}/{r.get('nSector','--')})")
    out['results']=res
    # matched-site comparisons
    for cl,ld,sites in [('chain-L6','ladder-2x3',6),('chain-L8','ladder-2x4',8),('chain-L10','ladder-2x5',10)]:
        c,l=res[cl],res[ld]
        log(f"== {sites} sites ==  F1: chain {c['F1']:.2f} vs ladder {l['F1']:.2f} | "
            f"chi: {c['chi']} -> {l['chi']} (x{l['chi']/max(c['chi'],1):.1f}) | "
            f"|S|: {c['Sdet']} -> {l['Sdet']} (x{l['Sdet']/max(c['Sdet'],1):.1f}) | "
            f"frac: {c.get('frac','--')} -> {l.get('frac','--')}")
    json.dump(out,open(_outpath('ladder_vs_chain.json'),'w'),indent=1)
    log("WROTE ladder_vs_chain.json")
