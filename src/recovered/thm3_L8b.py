# -*- coding: utf-8 -*-
"""The paper's flagship sampled configuration: L=8, U/t=8, eta=0.18t, top-TE-Born fraction 0.85
   (exactly the protocol of release/src/sampled_akw.py, which orders determinants by the
   time-averaged Born weight and keeps the top FR=0.85).  Site seed instead of momentum seed."""
import numpy as np, sys, os, time, json
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from thm3_core import hubbard_problem
t0=time.time(); log=lambda *a: print("[%6.1fs]"%(time.time()-t0),*a,flush=True)

L,U,eta=8,8.0,0.18
Hs,phi,E0=hubbard_problem(L=L,U=U); n=Hs.shape[0]
log(f"sector dim = {n}, ||phi||^2 = {phi@phi:.6f}")
E,V=np.linalg.eigh(Hs); coef=V.T@phi; poles=E-E0; wts=coef**2
lo,hi=poles.min()-1,poles.max()+1
grid=np.unique(np.concatenate([np.linspace(lo-2,hi+2,30001),
     -np.logspace(0,5,2000)[::-1], np.logspace(0,5,2000)]))
def spec(pl,ww,CH=256):
    out=np.zeros_like(grid)
    for i in range(0,len(pl),CH):
        p=pl[i:i+CH]; w=ww[i:i+CH]
        out+=((w[None,:]*(eta/np.pi))/((grid[:,None]-p[None,:])**2+eta**2)).sum(1)
    return out
A=spec(poles,wts); normA=np.trapz(A,grid); nrm2=float(phi@phi)
log(f"normA={normA:.6f} (should be {nrm2:.6f})")

K,dt=16,0.5; wc=np.zeros(n)
for k in range(K+1):
    vk=V@(np.exp(-1j*E*k*dt)*coef); wc+=np.abs(vk)**2
om=np.argsort(wc)[::-1]
out=[]
for FR in (0.50,0.70,0.85,0.95):
    S=np.sort(om[:int(round(FR*n))])
    mask=np.zeros(n,bool); mask[S]=True
    HS=Hs[np.ix_(S,S)]; Em,Um=np.linalg.eigh(HS); c=Um.T@phi[S]
    A_S=spec(Em-E0,c**2)
    LHS=float(np.trapz(np.abs(A-A_S),grid))
    phiS=np.where(mask,phi,0.0); wS=float(phiS@phiS)/nrm2; nphiQ=np.linalg.norm(phi-phiS)
    # K_S
    v=phi.copy(); KS=-1
    for j in range(4):
        if np.linalg.norm(v[~mask])/np.linalg.norm(v)>1e-10: break
        KS=j; v=Hs@v
    # Lambda_S (closed form)
    Uf=np.zeros((n,len(Em))); Uf[S,:]=Um; G=Hs@Uf; G[mask,:]=0.0
    Gram=G.T@G; D=Em[None,:]-Em[:,None]
    F=(4*eta*eta+2j*eta*D)/(D*D+4*eta*eta)
    Lam=float(np.sqrt(max(np.real(((c[:,None]*c[None,:])*Gram.T*F).sum()),0.0)))
    RHSp=2*nrm2*(1-wS); RHSn=(np.sqrt(nrm2)+np.sqrt(wS*nrm2))*(nphiQ+Lam/eta)
    r=dict(FR=FR,nS=int(len(S)),wS=wS,K_S=KS,relL1=LHS/normA,LHS=LHS,RHS_pap=RHSp,
           ratio=LHS/max(RHSp,1e-300),Lam=Lam,Lam_over_eta=Lam/eta,RHS_new=RHSn,
           ok_new=bool(LHS<=RHSn),LB=nrm2*(1-wS))
    out.append(r)
    log(f"FR={FR:.2f} |S|={r['nS']:5d} w_S={wS:.10f} K_S={KS:2d} relL1={r['relL1']:.5f} "
        f"L1={LHS:.5f} | RHSpap={RHSp:.3e} VIOLATED x{r['ratio']:.1f} | Lam/eta={Lam/eta:.2f} RHSnew={RHSn:.3f} ok={r['ok_new']}")
    del G,Gram,F,Uf
json.dump(out,open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"thm3_L8b.json"),'w'))
log("done")
