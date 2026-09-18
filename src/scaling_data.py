# -*- coding: utf-8 -*-
r"""SCALING toward the beyond-classical regime (real data for fig:scaling). For L=4,6,8 (8,12,16 qubits)
half-filled Hubbard U=8: (a) fermionic magic FAF of the ground state = 4*sum_i g_i(1-g_i) over spin-orbital
natural occupations g_i (exact identity for number-conserving states; validated vs full-space at L=4);
(b) subspace FRACTION of the (N+1) sector needed to reconstruct A(omega) to rel-L1<0.05 from the
time-evolved seed. Sector engine reused from akw_lanczos.py. Writes scaling_data.json."""
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
import akw_lanczos as AK
t0=time.time(); log=lambda *a: print(f"[{time.time()-t0:6.1f}s]",*a,flush=True)

def onerdm_up(psi_mat,S,idx,L):
    """up-spin 1-RDM Gamma_up[i,j]=<psi|c^dag_i c_j|psi> in the up-string basis (down traced out)."""
    rho=psi_mat@psi_mat.conj().T            # (Du,Du) reduced density matrix of up sector
    G=np.zeros((L,L))
    for a,m in enumerate(S):
        for j in range(L):
            if not (m>>j)&1: continue
            for i in range(L):
                if i==j: m2=m; sign=1.0
                elif (m>>i)&1: continue
                else:
                    m2=(m&~(1<<j))|(1<<i)
                    lo,hi=min(i,j),max(i,j)
                    mask=m & (((1<<hi)-1) ^ ((1<<(lo+1))-1))
                    sign=-1.0 if (bin(mask).count('1')&1) else 1.0
                G[i,j]+=float(rho[idx[m2],a].real)*sign
    return G

def faf_from_rdm(g):   # g: spin-orbital occupations in [0,1]; FAF=4 sum g(1-g) over all spin-orbitals
    return float(8.0*np.sum(g*(1.0-g)))     # x2 for up+down (equal at Sz=0)

def spectral_fraction(L,U,eta=0.15,K=18,dt=0.5,thr=0.05):
    nup=L//2; nd=L//2
    H,Du,Dd=AK.build_H_explicit(L,U,nup,nd)
    w,v=eigsh(H,k=1,which='SA',v0=_v0(H.shape[0])); E0=float(w[0]); psi=v[:,0]; psi_mat=psi.reshape(Du,Dd)
    Su,iu=AK.strings(L,nup); g=np.linalg.eigvalsh(onerdm_up(psi_mat,Su,iu,L)); g=np.clip(g,0,1)
    FAF=faf_from_rdm(g)
    # (N+1,Sz=+1) sector: nup+1 up, nd down. seed = c^dag_{0up}|GS>
    cdU=AK.cdag_map(L,nup)                      # up-string c^dag_j : nup -> nup+1
    Su1,iu1=AK.strings(L,nup+1)
    seed_mat=(cdU[0]@psi_mat)                   # (Du1, Dd)
    H1,Du1,Dd1=AK.build_H_explicit(L,U,nup+1,nd)
    seed=seed_mat.reshape(-1); nS=H1.shape[0]
    Em,Um=np.linalg.eigh(H1.toarray()); coef=Um.conj().T@seed
    poles=Em-E0; grid=np.linspace(poles.min()-1,poles.max()+1,600)
    def spec(pw,ww):
        A=np.zeros_like(grid)
        for a,b in zip(pw,ww): A+=b*(eta/np.pi)/((grid-a)**2+eta**2)
        return A
    A_ex=spec(poles,np.abs(coef)**2); nrm=np.trapezoid(A_ex,grid)
    # time-evolution union weight over configs, rank, grow subspace
    wc=np.zeros(nS)
    for k in range(K+1):
        vk=Um@(np.exp(-1j*Em*k*dt)*coef); wc+=np.abs(vk)**2
    order=np.argsort(wc)[::-1]; Harr=H1.toarray()
    def relL1_at(k):
        Sset=np.sort(order[:k]); HS=Harr[np.ix_(Sset,Sset)]
        Emk,Umk=np.linalg.eigh(HS); a=Umk.conj().T@seed[Sset]
        return float(np.trapezoid(np.abs(spec(Emk-E0,np.abs(a)**2)-A_ex),grid)/nrm)
    # coarse fractions, then bisection between the last-fail and first-pass k
    fracs=[0.3,0.45,0.6,0.7,0.78,0.85,0.92,1.0]; ks=sorted(set(max(2,int(f*nS)) for f in fracs))
    frac=1.0; klo=ks[0]
    for k in ks:
        if relL1_at(k)<thr:
            khi=k
            while khi-klo>max(2,nS//60):     # bisect for a tighter fraction
                km=(klo+khi)//2
                if relL1_at(km)<thr: khi=km
                else: klo=km
            frac=khi/nS; break
        klo=k
    return FAF, int(nS), frac, len(g)

out={'U':8.0,'points':[]}
for L in (4,6,8):
    FAF,nS,frac,nso=spectral_fraction(L,8.0)
    out['points'].append({'L':L,'qubits':2*L,'FAF':FAF,'FAF_per_site':FAF/L,'Np1_sector':nS,'frac':frac})
    log(f"L={L} ({2*L}q): FAF={FAF:.3f} (per-site {FAF/L:.3f})  (N+1)sector={nS}  frac@relL1<0.05={frac:.3f}")
json.dump(out,open(_outpath('scaling_data.json'),'w'),indent=1)
log("WROTE scaling_data.json")
