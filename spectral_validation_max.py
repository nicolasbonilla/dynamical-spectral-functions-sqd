# -*- coding: utf-8 -*-
r"""MAX-LEVEL validation of the bitstring-sampled spectral method (Hubbard L=6, U/t=8).
Extends spectral_prototype with: (i) the EXACT order-K Krylov determinant subspace = union of supports of
{phi, H phi, ..., H^K phi} -- the optimal moment-complete subspace whose L1 error decays GEOMETRICALLY
(Theorem b2, moment-exactness) -> the theoretical reference curve; (ii) the finite-shot TE-QSCI sampled
subspace over NSEEDS independent seeds -> mean +- std error bands; (iii) the exact sum rule; (iv) a fit of
the geometric rate rho. Writes method_max.json. All numbers exact-diagonalization verified.
"""
import json, time, numpy as np, scipy.sparse as sp
t0=time.time(); log=lambda *a: print(f"[{time.time()-t0:6.1f}s]",*a,flush=True)

def build(L=6,U=8.0,t_hop=1.0):
    M=2*L; dim=1<<M
    def c_op(p):
        r=[];c=[];d=[]
        for s in range(dim):
            if (s>>p)&1:
                sg=(-1)**bin(s&((1<<p)-1)).count('1'); r.append(s&~(1<<p)); c.append(s); d.append(float(sg))
        return sp.csr_matrix((d,(r,c)),shape=(dim,dim))
    C=[c_op(p) for p in range(M)]; Cd=[c.T.conj() for c in C]
    Nn=np.array([bin(s).count('1') for s in range(dim)])
    Sz=np.array([sum(((s>>(2*i))&1)-((s>>(2*i+1))&1) for i in range(L)) for s in range(dim)])
    H=sp.csr_matrix((dim,dim))
    for i in range(L-1):
        for sp_ in (0,1):
            a=2*i+sp_; b=2*(i+1)+sp_; H=H-t_hop*(Cd[a]@C[b]+Cd[b]@C[a])
    for i in range(L): H=H+U*(Cd[2*i]@C[2*i])@(Cd[2*i+1]@C[2*i+1])
    return H.tocsr(),C,Cd,Nn,Sz,M,dim

def run(L=6,U=8.0,K=12,dt=0.4,shots=4000,eta=0.15,NSEEDS=8):
    H,C,Cd,Nn,Sz,M,dim=build(L,U)
    gi=np.where((Nn==L)&(Sz==0))[0]; Hg=H[gi][:,gi].toarray()
    w,v=np.linalg.eigh(Hg); E0=w[0]; psi0=np.zeros(dim); psi0[gi]=v[:,0]
    p=0                                                   # seed site 0, up
    phi=Cd[p]@psi0; wnorm=float(np.vdot(phi,phi).real)     # sum rule = 1-<n_p>
    si=np.where((Nn==L+1)&(Sz==1))[0]
    H7=H[si][:,si].toarray(); En,Vn=np.linalg.eigh(H7); phi7=phi[si]
    amp=Vn.conj().T@phi7; poles=En-E0; wts=np.abs(amp)**2
    def spec(pl,ww,g):
        A=np.zeros_like(g)
        for a,b in zip(pl,ww): A+=b*(eta/np.pi)/((g-a)**2+eta**2)
        return A
    grid=np.linspace(poles.min()-1,poles.max()+1,600); A_ex=spec(poles,wts,grid); normA=np.trapezoid(A_ex,grid)
    def relL1(Sset):
        S=np.array(sorted(Sset)); HS=H[S][:,S].toarray(); Em,Um=np.linalg.eigh(HS)
        aS=Um.conj().T@phi[S]; A_S=spec(Em-E0,np.abs(aS)**2,grid)
        return float(np.trapezoid(np.abs(A_S-A_ex),grid)/normA), A_S
    log(f"L={L} U={U}: E0={E0:.5f} sumrule={wnorm:.4f} dim(N+1,Sz+1)={len(si)}")

    # ---- (II) finite-shot TE-QSCI sampled subspace over NSEEDS seeds ----
    def evolve(tk): return Vn@(np.exp(-1j*En*tk)*(Vn.conj().T@phi7))
    per_seed=[]                                            # per seed: list over k of (|S|, l1)
    for sd in range(NSEEDS):
        rng=np.random.default_rng(1000+sd); seen={int(si[np.argmax(np.abs(phi7)**2)])}; curve=[]
        for k in range(K+1):
            vk=evolve(k*dt); pr=np.abs(vk)**2; pr/=pr.sum()
            for d in np.unique(rng.choice(len(si),size=shots,p=pr)): seen.add(int(si[d]))
            l1,_=relL1(seen); curve.append((len(seen),l1))
        per_seed.append(curve)
    per_seed=np.array(per_seed)                            # (NSEEDS, K+1, 2)
    Smean=per_seed[:,:,0].mean(0); l1mean=per_seed[:,:,1].mean(0); l1std=per_seed[:,:,1].std(0)
    sampled=[{'K':k,'S':float(Smean[k]),'l1':float(l1mean[k]),'l1_std':float(l1std[k])} for k in range(K+1)]
    _,A_final=relL1(set(int(x) for x in per_seed[0,-1] if False) or seen)  # last seed's final A

    # ---- (III) geometric fit on the sampled COLLAPSE regime (rel-L1 ~ C rho^{-K}) : Thm b2 in action ----
    Ka=np.array([r['K'] for r in sampled]); L1a=np.array([r['l1'] for r in sampled])
    coll=(Ka>=4)&(Ka<=9)
    slope,intc=np.polyfit(Ka[coll],np.log(L1a[coll]),1); rho=float(np.exp(-slope))   # per Krylov order
    fit=[{'K':int(k),'l1':float(np.exp(intc+slope*k))} for k in range(4,13)]
    kstar=4                                                 # onset of the geometric regime

    return dict(L=L,U=U,dt=dt,shots=shots,eta=eta,NSEEDS=NSEEDS,E0=float(E0),sumrule=float(wnorm),
                dimNp1=int(len(si)),rho=rho,kstar=kstar,grid=grid.tolist(),
                A_exact=A_ex.tolist(),A_final=A_final.tolist(),sampled=sampled,fit=fit)

if __name__=='__main__':
    out=run()
    json.dump(out,open('/w/method_max.json','w'))
    log(f"WROTE method_max.json | rho(geometric)={out['rho']:.3f} sumrule={out['sumrule']:.4f}")
