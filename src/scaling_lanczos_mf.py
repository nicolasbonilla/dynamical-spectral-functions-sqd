# -*- coding: utf-8 -*-
r"""MATRIX-FREE, STAGED, RESUMABLE scaling computation to push to L=14 (28 qubits) on a RAM-limited box.
NEVER builds an explicit Hamiltonian (uses AK.sector_H LinearOperator) and NEVER materializes a submatrix.
Runs ONE stage per process so the OS reclaims all memory between stages; each stage checkpoints to ckpt/.
Stages:  gs (ground state + FAF)  ->  exact (exact A(w))  ->  rank (config ranking)  ->  frac (bisection).
Usage:  python scaling_lanczos_mf.py <L> <stage>       e.g.  python scaling_lanczos_mf.py 14 gs
Memory economy: small eigsh ncv; memory-light Krylov exponential propagator (few vectors, substeps).
Validated on L=10 (must reproduce frac=0.3625 of the dense/expm reference) before trusting L=14."""
import json, sys, os, time, gc, numpy as np
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
from scipy.sparse.linalg import eigsh, LinearOperator
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
from scipy.linalg import expm as dense_expm
import akw_lanczos as AK
t0=time.time(); log=lambda *a: print(f"[{time.time()-t0:7.1f}s]",*a,flush=True)
CK='ckpt'; os.makedirs(CK,exist_ok=True)
eta=0.15; ngrid=600; nl=260

def onerdm_up(psi_mat,S,idx,L):
    rho=psi_mat@psi_mat.conj().T; G=np.zeros((L,L))
    for a,m in enumerate(S):
        for j in range(L):
            if not (m>>j)&1: continue
            for i in range(L):
                if i==j: m2=m; sign=1.0
                elif (m>>i)&1: continue
                else:
                    m2=(m&~(1<<j))|(1<<i); lo,hi=min(i,j),max(i,j)
                    mask=m & (((1<<hi)-1) ^ ((1<<(lo+1))-1)); sign=-1.0 if (bin(mask).count('1')&1) else 1.0
                G[i,j]+=float(rho[idx[m2],a].real)*sign
    return G

def sectorHmv(L,U,nup,ndn):
    """matrix-free matvec for the (nup,ndn) Hubbard sector; returns (matvec, dim, Du, Dd, diag)."""
    Tu,Su,iu=AK.hop(L,nup); Td,Sd,idd=AK.hop(L,ndn)
    Du=len(Su); Dd=len(Sd)
    upocc=np.array([[ (m>>i)&1 for i in range(L)] for m in Su],dtype=np.float64)
    dnocc=np.array([[ (m>>i)&1 for i in range(L)] for m in Sd],dtype=np.float64)
    diagU=U*(upocc@dnocc.T).ravel()
    def mv(x):
        X=x.reshape(Du,Dd); Y=Tu@X + (Td@X.T).T
        return Y.ravel()+diagU*x
    return mv,Du*Dd,Du,Dd,diagU,(Su,iu)

def haydock_grid(mv,v0,E0,grid,nlanc=None):
    z=grid+1j*eta
    return AK.haydock(lambda x: mv(x)-E0*x, v0, min(nlanc or nl,len(v0)), z)

def krylov_expm(mv,v,dt,m):
    """e^{-i dt H} v via a small Lanczos basis (m vectors); memory = m complex vectors."""
    n=v.shape[0]; V=np.empty((m,n),dtype=complex); beta=np.linalg.norm(v)
    if beta<1e-14: return v.copy()
    V[0]=v/beta; al=np.zeros(m); be=np.zeros(m)
    w=mv(V[0]); a=np.vdot(V[0],w).real; al[0]=a; w=w-a*V[0]; mm=m
    for j in range(1,m):
        b=np.linalg.norm(w); be[j-1]=b
        if b<1e-12: mm=j; break
        V[j]=w/b; w=mv(V[j]); a=np.vdot(V[j],w).real; al[j]=a; w=w-a*V[j]-b*V[j-1]
    T=np.diag(al[:mm])+np.diag(be[:mm-1],1)+np.diag(be[:mm-1],-1)
    E=dense_expm(-1j*dt*T)[:,0]
    out=(beta*(V[:mm].T@E)); del V; return out

def stage_gs(L,U=8.0):
    nup=nd=L//2
    log(f"L={L} GS: sector ({nup},{nd}) matrix-free eigsh ...")
    mv,dim,Du,Dd,diagU,(Su,iu)=sectorHmv(L,U,nup,nd)
    H=LinearOperator((dim,dim),matvec=mv,dtype=float)
    ncv=max(8,2*1+2)   # tiny Krylov basis for the k=1 ground state -> low memory
    w,v=eigsh(H,k=1,which='SA',ncv=12,maxiter=5000,tol=1e-9,v0=_v0(dim))
    E0=float(w[0]); psi=v[:,0].copy(); del v; gc.collect()
    psi_mat=psi.reshape(Du,Dd)
    g=np.linalg.eigvalsh(onerdm_up(psi_mat,Su,iu,L)); g=np.clip(g,0,1); FAF=float(8.0*np.sum(g*(1-g)))
    np.savez(f'{CK}/L{L}_gs.npz',E0=E0,FAF=FAF,psi=psi,nup=nup,nd=nd)
    log(f"L={L} GS DONE: E0={E0:.5f} FAF={FAF:.4f} n_up={g.sum():.3f} -> {CK}/L{L}_gs.npz")

def stage_exact(L,U=8.0):
    d=np.load(f'{CK}/L{L}_gs.npz'); E0=float(d['E0']); psi=d['psi']; nup=int(d['nup']); nd=int(d['nd'])
    _,_,Du0,Dd0,_,_=sectorHmv(L,U,nup,nd)   # only to get Du0/Dd0 shape (cheap)
    psi_mat=psi.reshape(Du0,Dd0); del psi
    cdU=AK.cdag_map(L,nup); seed=(cdU[0]@psi_mat).reshape(-1).astype(complex); del psi_mat; gc.collect()
    mv1,dim1,Du1,Dd1,diagU1,_=sectorHmv(L,U,nup+1,nd); nS=dim1
    log(f"L={L} EXACT: (N+1) dim={nS}  seed_norm2={np.vdot(seed,seed).real:.4f}  Haydock ...")
    H1=LinearOperator((nS,nS),matvec=mv1,dtype=float)
    emin=float(eigsh(H1,k=1,which='SA',ncv=10,return_eigenvectors=False,tol=1e-6,v0=_v0(H1.shape[0]))[0])
    emax=float(eigsh(H1,k=1,which='LA',ncv=10,return_eigenvectors=False,tol=1e-6,v0=_v0(H1.shape[0]))[0])
    grid=np.linspace(emin-E0-1.0,emax-E0+1.0,ngrid)
    A_ex=(-haydock_grid(mv1,seed,E0,grid).imag/np.pi); nrm=float(np.trapz(np.abs(A_ex),grid))
    np.savez(f'{CK}/L{L}_exact.npz',E0=E0,seed=seed,grid=grid,A_ex=A_ex,nrm=nrm,nS=nS,nup=nup,nd=nd)
    log(f"L={L} EXACT DONE: window[{grid[0]:.2f},{grid[-1]:.2f}] integral={np.trapz(A_ex,grid):.4f} -> saved")

def stage_rank(L,U=8.0,K=18,dt=0.5,sub=16):
    d=np.load(f'{CK}/L{L}_exact.npz'); seed=d['seed']; nup=int(d['nup']); nd=int(d['nd']); nS=int(d['nS'])
    mv1,_,_,_,_,_=sectorHmv(L,U,nup+1,nd)
    ddt=dt/sub; m=6                       # finer substep + smaller Krylov basis -> lower RAM (fits ~1.4GB @ L=14)
    v=seed/np.sqrt(np.vdot(seed,seed).real); del seed; gc.collect()
    wc=(np.abs(v)**2).astype(np.float32); nsteps=K*sub   # float32 accumulator halves its footprint
    log(f"L={L} RANK: {nsteps} substeps (ddt={ddt:.4f}, Krylov m={m}) ...")
    for s in range(1,nsteps+1):
        v=krylov_expm(mv1,v,ddt,m)
        if s%sub==0: wc+=(np.abs(v)**2).astype(np.float32)
        if s%(sub*3)==0: log(f"L={L} RANK: t={s*ddt:.1f}/{K*dt}")
    order=np.argsort(wc)[::-1].astype(np.int32); del v,wc; gc.collect()
    np.save(f'{CK}/L{L}_order.npy',order)
    log(f"L={L} RANK DONE: saved order ({order.shape[0]} configs) -> {CK}/L{L}_order.npy")

def stage_frac(L,U=8.0,thr=0.05,nl_sub=150):
    de=np.load(f'{CK}/L{L}_exact.npz'); E0=float(de['E0']); seed=de['seed']; grid=de['grid']
    nS=int(de['nS']); nup=int(de['nup']); nd=int(de['nd'])
    order=np.load(f'{CK}/L{L}_order.npy')
    mv1,_,_,_,_,_=sectorHmv(L,U,nup+1,nd)
    # recompute the exact A with the SAME reduced Lanczos depth as the subspaces (consistent comparison, fast)
    A_ex=(-haydock_grid(mv1,seed,E0,grid,nl_sub).imag/np.pi); nrm=float(np.trapz(np.abs(A_ex),grid))
    log(f"L={L} FRAC: exact A recomputed at nl_sub={nl_sub} (integral={np.trapz(A_ex,grid):.4f})")
    def relL1_at(k):
        Sset=np.sort(order[:k])
        def rmv(xs):
            xf=np.zeros(nS,dtype=complex); xf[Sset]=xs; return mv1(xf)[Sset]
        A_s=(-haydock_grid(rmv,seed[Sset],E0,grid,nl_sub).imag/np.pi)
        return float(np.trapz(np.abs(A_s-A_ex),grid)/nrm)
    fracs=[0.02,0.05,0.08,0.12,0.16,0.20,0.25,0.30,0.40,0.50,0.60,0.70,0.78,0.85,0.92,1.0]
    ks=sorted(set(max(2,int(f*nS)) for f in fracs)); frac=1.0; klo=2
    for k in ks:
        r=relL1_at(k); log(f"L={L} FRAC: frac={k/nS:.3f} (k={k}) relL1={r:.4f}")
        if r<thr:
            khi=k
            while khi-klo>max(2,nS//200):
                km=(klo+khi)//2
                if relL1_at(km)<thr: khi=km
                else: klo=km
            frac=khi/nS; break
        klo=k
    dg=np.load(f'{CK}/L{L}_gs.npz'); FAF=float(dg['FAF'])
    p={'L':L,'qubits':2*L,'FAF':FAF,'FAF_per_site':FAF/L,'Np1_sector':nS,'frac':float(frac),
       'E0':E0,'method':'lanczos-haydock-matrixfree'}
    # Resolved from this file's own location, not from the caller's cwd.
    fn=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    'data','scaling_data.json')
    data=json.load(open(fn)); data['points']=[q for q in data['points'] if q['L']!=L]+[p]
    data['points'].sort(key=lambda q:q['L']); json.dump(data,open(fn,'w'),indent=1)
    log(f"L={L} FRAC DONE: FAF={FAF:.4f} frac={frac:.4f} -> appended to {fn}")

if __name__=='__main__':
    L=int(sys.argv[1]); stage=sys.argv[2]
    {'gs':stage_gs,'exact':stage_exact,'rank':stage_rank,'frac':stage_frac}[stage](L)
