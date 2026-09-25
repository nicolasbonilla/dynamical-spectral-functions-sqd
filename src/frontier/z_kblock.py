# -*- coding: utf-8 -*-
"""The momentum-block route of Sec. V.E, measured.

A momentum-resolved probe phi_k = c^dag_{k,up}|Psi_0> is an eigenvector of lattice
translation, so K(H, phi_k) lies inside one total-momentum block of the (N+1) sector, of
exactly dim/L determinants.  This script rotates the sector to the plane-wave determinant
basis (Slater-determinant overlaps, unitarity checked), confirms H is block diagonal there,
measures the coordinate support of K(H, phi_k) inside the block from the exact projector,
ranks determinants by the accumulated Born weight of the time-evolved probe (K = 16,
dt = 0.5, the sampled_akw protocol of c0_lib.born_order), and reports the smallest |S| at
which the leakage branch of Theorem 1 drops below the trivial bound (eta = 0.18 t, U/t = 8).
Uses the deposited builder (akw_lanczos) and certificate (leakage_certificate.lambda_leak).

    python src/frontier/z_kblock.py 6      (seconds)
    python src/frontier/z_kblock.py 8      (minutes; dense sector of 3920)
writes data/c3_frontier/adversarial/z_kblock_L<L>.json

One count is threshold-sensitive and is NOT quoted in the paper: |suppK| is read from the
projector onto the eigenvectors whose coefficient exceeds 1e-10 of the probe norm, each
degenerate group normalized to one, so a roundoff-level coefficient in another block counts
as a full vector there.  At L = 8, k = 2(2pi/8) this gives |suppK| = 980, two blocks,
although the half-filled ground state is non-degenerate (gap 0.24 t) and phi_k lies in one
block.  The certificate does not share the artefact -- it uses the amplitudes themselves --
and its first non-vacuous |S| (483 of 490 there) is what the paper quotes.
"""
import os, json
import sys, numpy as np, scipy.sparse as sp
from itertools import combinations
from scipy.sparse.linalg import eigsh
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import akw_lanczos as AK
from leakage_certificate import lambda_leak
def strings(L,n): return sorted(sum(1<<b for b in c) for c in combinations(range(L),n))
def occ(m,L): return [i for i in range(L) if (m>>i)&1]
def Mrot(L,n,U):
    S=strings(L,n); M=np.zeros((len(S),len(S)),complex)
    for a,J in enumerate(S):
        oj=occ(J,L)
        for b,K in enumerate(S):
            M[a,b]=np.linalg.det(U[np.ix_(oj,occ(K,L))])
    return M,S
L=int(sys.argv[1]); Uint=8.0; eta=0.18
nup=nd=L//2
H0,Du,Dd=AK.build_H_explicit(L,Uint,nup,nd)
E0,V=eigsh(H0,k=1,which='SA',v0=np.random.default_rng(3).standard_normal(H0.shape[0])); E0=E0[0]
Psi=V[:,0].reshape(Du,Dd)
cd=AK.cdag_map(L,nup)
H1,Du1,Dd1=AK.build_H_explicit(L,Uint,nup+1,nd); H1=H1.toarray(); D=H1.shape[0]
Uo=np.array([[np.exp(1j*2*np.pi*n*j/L)/np.sqrt(L) for n in range(L)] for j in range(L)])
Mu,Su=Mrot(L,nup+1,Uo); Md,Sd=Mrot(L,nd,Uo)
M=np.kron(Mu,Md)
print('unitary defect',np.abs(M.conj().T@M-np.eye(D)).max())
Hk=M.conj().T@H1@M
# total momentum label of each k-determinant
P=np.array([ (sum(occ(a,L))+sum(occ(b,L)))%L for a in Su for b in Sd])
off=np.abs(Hk[P[:,None]!=P[None,:]]).max(); print('L',L,'D',D,'max offblock |H_k|',off, 'block sizes',np.bincount(P))
Em,Um=np.linalg.eigh(Hk)
RES=[]
for kn in range(L//2+1):
    phi_r=sum(np.exp(1j*2*np.pi*kn*j/L)/np.sqrt(L)*(cd[j]@Psi).ravel() for j in range(L))
    phik=M.conj().T@phi_r; n2=np.vdot(phik,phik).real
    blk=np.unique(P[np.abs(phik)>1e-10])
    # Krylov support via projector diag
    c=Um.conj().T@phik; live=np.abs(c)>1e-10*np.sqrt(n2)
    E=Em[live];Uv=Um[:,live];cv=c[live];cols=[];i=0
    while i<len(E):
        j=i
        while j+1<len(E) and E[j+1]-E[i]<1e-9: j+=1
        v=Uv[:,i:j+1]@cv[i:j+1]; nv=np.linalg.norm(v)
        if nv>1e-10: cols.append(v/nv)
        i=j+1
        
    Kb=np.array(cols).T; dg=(np.abs(Kb)**2).sum(1)
    suppK=int((dg>1e-12).sum()); bsz=int((P==blk[0]).sum())
    # Born ranking of time-evolved phi_k in k-basis, K=16, dt=0.5
    v0=phik/np.sqrt(n2); co=Um.conj().T@v0; wc=np.zeros(D)
    for s in range(17): wc+=np.abs(Um@(np.exp(-1j*Em*s*0.5)*co))**2
    order=np.argsort(-wc, kind='stable')
    first=None; rows=[]
    for k in range(5, D+1):
        S=np.sort(order[:k]); Q=np.setdiff1d(np.arange(D),S)
        w=np.vdot(phik[S],phik[S]).real/n2
        th,Ur=np.linalg.eigh(Hk[np.ix_(S,S)]); cc=Ur.conj().T@phik[S]
        lam,_=lambda_leak(Hk[np.ix_(Q,S)],th,cc,Ur,eta)
        B=(1+np.sqrt(w))*(np.sqrt(max(1-w,0))+lam/(np.sqrt(n2)*eta)); triv=1+w
        if B<triv: first=(k,w,B,triv); break
        if k>bsz+10: break
    print('k=%d*2pi/%d block=%s size %d  |suppK|=%d (%.3f of block)  dimK=%d  first non-vacuous |S|=%s'%(kn,L,blk,bsz,suppK,suppK/bsz,Kb.shape[1],first))
    RES.append(dict(k_index=kn, block=int(blk[0]), block_size=bsz, suppK=suppK, dimK=int(Kb.shape[1]),
                    first_nonvacuous=(None if first is None else int(first[0])),
                    w_S=(None if first is None else float(first[1])), bound=(None if first is None else float(first[2]))))
OUT=os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),'data','c3_frontier','adversarial')
json.dump(dict(L=L,U=Uint,eta=eta,D=D,rows=RES),open(os.path.join(OUT,'z_kblock_L%d.json'%L),'w'),indent=1)
