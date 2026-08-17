# -*- coding: utf-8 -*-
r"""CONVENTION-FREE complexity: nonfreeness of the correction vector across frequency (the Z-criterion test).
Nonfreeness of a pure state = entropy of the Gaussian state with the same 1-RDM = -sum_i H2(g_i), g_i = natural
occupations (eigenvalues of the 1-RDM). It is ORBITAL-ROTATION INVARIANT (Gottlieb-Mauser; = min-over-bases
orbital correlation, Ding-Schilling) -> immune to the JW/basis-ordering problem that broke our entanglement
numbers. Test: does nonfreeness(|chi(w)>) DIP in the coherent quasiparticle window and RISE in the incoherent
Hubbard satellite? That is the physical easy/hard criterion. Compares to A(w). Sizes 2x3,4,5. Writes nonfree.txt.
"""
import time, numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh, gmres
import gate1_ladder as G, akw_lanczos as AK
t0=time.time()
def emit(s): print(f"[{time.time()-t0:6.1f}s] {s}",flush=True); open('/w/nonfree.txt','a').write(s+"\n")
def H2(x):  # binary entropy in bits, safe
    x=np.clip(x,1e-14,1-1e-14); return -(x*np.log2(x)+(1-x)*np.log2(1-x))

def onerdm_up(chi_md, Su1, iu1, L):
    """gamma^up[i,j] = <chi| a^dag_i a_j |chi>, chi_md shape (len(Su1), Ddn). Down is a spectator."""
    rho=chi_md@chi_md.conj().T                       # (Du1,Du1) up reduced (trace down)
    g=np.zeros((L,L),complex)
    for a,m in enumerate(Su1):
        for j in range(L):
            if not (m>>j)&1: continue
            m1=m^(1<<j); sj=(-1)**bin(m&((1<<j)-1)).count('1')
            for i in range(L):
                if (m1>>i)&1: continue
                m2=m1|(1<<i); si=(-1)**bin(m1&((1<<i)-1)).count('1')
                g[i,j]+=si*sj*rho[iu1[m2],a]
    return g

def occupations(chi, Su1,iu1,Sd1,idd1,L):
    Du1,Dd=len(Su1),len(Sd1); md=chi.reshape(Du1,Dd)
    gu=onerdm_up(md,Su1,iu1,L)
    gd=onerdm_up(md.T.copy(),Sd1,idd1,L)             # symmetric role for down (spectator=up)
    ou=np.clip(np.linalg.eigvalsh((gu+gu.conj().T)/2).real,0,1)
    od=np.clip(np.linalg.eigvalsh((gd+gd.conj().T)/2).real,0,1)
    return np.concatenate([ou,od])

def run(Lr,U=8.0,hole=2,eta=0.15,nw=9,K=140):
    L=2*Lr; bonds=G.ladder_bonds(Lr); N=L-hole; nup=N//2; ndn=N-nup
    H,Su,iu,Sd,idd,Du,Dd=G.build_H(L,U,nup,ndn,bonds)
    w,v=eigsh(H,k=1,which='SA'); E0=float(w[0]); psi=v[:,0].reshape(Du,Dd)
    Su1,iu1=G.strings(L,nup+1); r=[];c=[];val=[]
    for a,m in enumerate(Su):
        if not (m>>0)&1: r.append(iu1[m|1]); c.append(a); val.append(1.0)
    Cd=sp.csr_matrix((val,(r,c)),shape=(len(Su1),Du)); phi=(Cd@psi).reshape(-1)
    H1,Su1b,iu1b,Sd1,idd1,Du1,Dd1=G.build_H(L,U,nup+1,ndn,bonds); nS=H1.shape[0]
    phi=phi/np.linalg.norm(phi); I=sp.identity(nS,format='csc'); Hc=H1.tocsr()
    grid=np.linspace(-6,24,1400); A=-(1/np.pi)*np.imag(AK.haydock(lambda x:H1@x,phi,K,grid+E0+1j*eta))
    supp=grid[A>0.02*A.max()]; wpts=np.linspace(supp.min(),supp.max(),nw)
    prof=[]
    for om in wpts:
        x,info=gmres(Hc-((E0+om)+1j*eta)*I,phi,rtol=1e-7,atol=1e-10,restart=200,maxiter=600); x=x/np.linalg.norm(x)
        g=occupations(x,Su1b,iu1b,Sd1,idd1,L); nf=float(H2(g).sum()); Aom=float(np.interp(om,grid,A))
        prof.append((float(om),Aom,nf))
    lo=min(prof,key=lambda t:t[2]); hi=max(prof,key=lambda t:t[2])
    emit(f"2x{Lr}: nonfreeness(w) [bits] "+" ".join(f"w{om:.1f}:{nf:.2f}(A{Aom:.2f})" for om,Aom,nf in prof))
    emit(f"==> 2x{Lr}: nonfreeness MIN={lo[2]:.2f}@w{lo[0]:.1f}(coherent)  MAX={hi[2]:.2f}@w{hi[0]:.1f}(incoherent)  ratio={hi[2]/max(lo[2],1e-9):.2f}")

if __name__=='__main__':
    open('/w/nonfree.txt','w').write("")
    for Lr in (3,4,5): run(Lr)
    emit("DONE nonfreeness")
