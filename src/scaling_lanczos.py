# -*- coding: utf-8 -*-
r"""SCALING toward beyond-classical (fig:scaling) — LANCZOS/HAYDOCK version that scales past the dense-ED
frontier (L=8) to L=10 (20q) and L=12 (24q). NEVER densifies: ground state by eigsh, exact and subspace
A(omega) by Haydock continued fraction with a MATRIX-FREE RESTRICTED matvec (embed -> full H.dot -> extract;
no submatrix is ever materialized), config ranking by Krylov time evolution (expm_multiply).
Computes ONE or more L (argv), checkpoints per L into scaling_data.json. Validate on L=4,6,8 first:
must reproduce FAF 6.61/9.39/12.81 and frac 0.92/0.83/0.56 of the dense reference before trusting L=10,12.
Reuses akw_lanczos primitives (build_H_explicit sparse, haydock, strings, cdag_map)."""
import json, sys, time, os, numpy as np, scipy.sparse as sp
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
from scipy.sparse.linalg import eigsh, expm_multiply
# --- DETERMINISTIC ARPACK START VECTOR (added 2026-09-18) -------------------
# eigsh() with no v0= lets ARPACK draw its own random start vector from an
# unseeded generator.  E0 then converges to a slightly different point every run
# (the last few digits move), and every rel-L1 downstream moves with it: two runs
# of this script on the same machine did NOT agree digit-for-digit.  Nothing about
# the physics changes -- the eigenpair is the same to ARPACK's tolerance -- but the
# deposit must be bit-reproducible, so the start vector is now fixed.
# Deliberately NOT a numpy global seed: this touches only the ARPACK start vector.
def _v0(n):
    """Fixed, dimension-dependent ARPACK start vector (never orthogonal to the GS)."""
    return np.random.default_rng(20260918).standard_normal(n)
# ---------------------------------------------------------------------------
import akw_lanczos as AK
t0=time.time(); log=lambda *a: print(f"[{time.time()-t0:7.1f}s]",*a,flush=True)

def onerdm_up(psi_mat,S,idx,L):
    """up-spin 1-RDM Gamma_up[i,j]=<psi|c^dag_i c_j|psi> in the up-string basis (down traced out)."""
    rho=psi_mat@psi_mat.conj().T
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

def faf_from_rdm(g):   # FAF = 4 sum_{spin-orbitals} g(1-g); x2 for up+down equal at Sz=0
    return float(8.0*np.sum(g*(1.0-g)))

def spec(Hmv,seed,E0,grid,eta,nl):
    """A(omega) = -Im<seed|(w+i.eta-(H-E0))^{-1}|seed>/pi via Haydock; seed need not be normalized."""
    z=grid+1j*eta
    G=AK.haydock(lambda x: Hmv(x)-E0*x, seed, nl, z)
    return -G.imag/np.pi

def compute(L,U=8.0,eta=0.15,K=18,dt=0.5,thr=0.05,nl=260,ngrid=600):
    nup=L//2; nd=L//2
    log(f"L={L}: building GS sector ({nup},{nd}) ...")
    H0,Du0,Dd0=AK.build_H_explicit(L,U,nup,nd)
    w,v=eigsh(H0,k=1,which='SA',v0=_v0(H0.shape[0])); E0=float(w[0]); psi=v[:,0]; psi_mat=psi.reshape(Du0,Dd0)
    Su0,iu0=AK.strings(L,nup)
    g=np.linalg.eigvalsh(onerdm_up(psi_mat,Su0,iu0,L)); g=np.clip(g,0,1); FAF=faf_from_rdm(g)
    log(f"L={L}: E0={E0:.5f}  FAF={FAF:.4f}  (natural occ sum check n_up={g.sum():.3f})")
    # (N+1, Sz=+1) sector: nup+1 up, nd down ; seed = c^dag_{0,up}|GS>
    cdU=AK.cdag_map(L,nup)
    seed_mat=(cdU[0]@psi_mat)
    H1,Du1,Dd1=AK.build_H_explicit(L,U,nup+1,nd)
    nS=H1.shape[0]; seed=seed_mat.reshape(-1).astype(complex)
    log(f"L={L}: (N+1) sector dim nS={nS}  seed_norm2={np.vdot(seed,seed).real:.4f}")
    # spectral window from extremal eigenvalues of H1
    emin=float(eigsh(H1,k=1,which='SA',return_eigenvectors=False,v0=_v0(H1.shape[0]))[0])
    emax=float(eigsh(H1,k=1,which='LA',return_eigenvectors=False,v0=_v0(H1.shape[0]))[0])
    grid=np.linspace(emin-E0-1.0,emax-E0+1.0,ngrid)
    H1mv=lambda x: H1.dot(x)
    A_ex=spec(H1mv,seed,E0,grid,eta,min(nl,nS)); nrm=np.trapz(np.abs(A_ex),grid)
    log(f"L={L}: exact A(w) done (window [{grid[0]:.2f},{grid[-1]:.2f}], integral={np.trapz(A_ex,grid):.4f})")
    # config ranking by Krylov time-evolution union weight
    vt=seed/np.sqrt(np.vdot(seed,seed).real); wc=np.abs(vt)**2
    for k in range(1,K+1):
        vt=expm_multiply(-1j*dt*H1,vt); wc+=np.abs(vt)**2
    order=np.argsort(wc)[::-1]
    log(f"L={L}: config ranking done (top weight {wc[order[0]]:.3e})")
    # restricted-matvec Haydock: no submatrix ever materialized
    def relL1_at(k):
        Sset=np.sort(order[:k])
        def rmv(xs):
            xf=np.zeros(nS,dtype=complex); xf[Sset]=xs
            return H1.dot(xf)[Sset]
        A_s=spec(rmv,seed[Sset],E0,grid,eta,min(nl,k))
        return float(np.trapz(np.abs(A_s-A_ex),grid)/nrm)
    # coarse grid spans LOW fractions too (large systems need only a small fraction); klo starts at 2 so
    # the search can descend below the first passing point (fixes the upper-bound bug when frac<0.3).
    fracs=[0.02,0.05,0.08,0.12,0.16,0.20,0.25,0.30,0.40,0.50,0.60,0.70,0.78,0.85,0.92,1.0]
    ks=sorted(set(max(2,int(f*nS)) for f in fracs))
    frac=1.0; klo=2
    for k in ks:
        r=relL1_at(k); log(f"L={L}:   frac={k/nS:.3f} (k={k})  relL1={r:.4f}")
        if r<thr:
            khi=k
            while khi-klo>max(2,nS//200):
                km=(klo+khi)//2
                if relL1_at(km)<thr: khi=km
                else: klo=km
            frac=khi/nS; break
        klo=k
    log(f"L={L}: DONE  FAF={FAF:.4f}  frac@relL1<{thr}={frac:.4f}  nS={nS}")
    return {'L':L,'qubits':2*L,'FAF':FAF,'FAF_per_site':FAF/L,'Np1_sector':int(nS),'frac':float(frac),
            'E0':E0,'method':'lanczos-haydock'}

if __name__=='__main__':
    Ls=[int(x) for x in sys.argv[1:]] or [4,6,8]
    # Resolved from this file's own location, not from the caller's cwd: running it
    # from anywhere else used to start a FRESH scaling_data.json there (line below
    # falls back to an empty dict) and leave data/scaling_data.json untouched.
    fn=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    'data','scaling_data.json')
    ref={4:(6.607,0.9167),6:(9.385,0.83),8:(12.814,0.5625)}   # dense-ED reference for validation
    data=json.load(open(fn)) if os.path.exists(fn) else {'U':8.0,'points':[]}
    data.setdefault('points',[])
    for L in Ls:
        p=compute(L)
        if L in ref:
            fe,fr=ref[L]; ok=abs(p['FAF']-fe)<0.05 and abs(p['frac']-fr)<0.06
            log(f"VALIDATE L={L}: FAF {p['FAF']:.3f} vs {fe} | frac {p['frac']:.3f} vs {fr}  -> {'OK' if ok else 'MISMATCH!!'}")
        data['points']=[q for q in data['points'] if q['L']!=L]+[p]
        data['points'].sort(key=lambda q:q['L'])
        json.dump(data,open(fn,'w'),indent=1)
        log(f"checkpoint saved: L={L} -> {fn}")
    log("ALL DONE")
