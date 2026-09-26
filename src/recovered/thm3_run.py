# -*- coding: utf-8 -*-
import numpy as np, json, time, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thm3_core import Spec, hubbard_problem, log

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "thm3_results.json")

def subspaces(sp_, rng):
    n = sp_.n; fam = {}
    w = sp_.phi**2
    order_w = np.argsort(w)[::-1]
    # (a) top-Born (static seed weight)
    for f in [0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.667,0.7,0.8,0.9,0.95,0.99,1.0]:
        k=max(1,int(round(f*n))); fam[("topPhi", f)] = order_w[:k]
    # (b) top time-averaged Born mixture (exactly sampled_akw.py's oracle ordering)
    K,dt=16,0.5
    coef=sp_.coef; wc=np.zeros(n)
    for k in range(K+1):
        vk=sp_.V@(np.exp(-1j*sp_.E*k*dt)*coef); wc+=np.abs(vk)**2
    om=np.argsort(wc)[::-1]
    for f in [0.1,0.2,0.3,0.5,0.7,0.85,0.9,0.95,0.99]:
        k=max(1,int(round(f*n))); fam[("topTE", f)] = om[:k]
    # (c) uniform random
    for f in [0.1,0.3,0.5,0.7,0.9]:
        for s in range(4):
            k=max(1,int(round(f*n))); fam[("rand%d"%s, f)] = rng.choice(n,size=k,replace=False)
    # (d) TE-QSCI finite-shot sampling (protocol of spectral_validation_max.py)
    shots=4000; seen={int(np.argmax(np.abs(sp_.phi)**2))}
    for k in range(0,13):
        vk=sp_.V@(np.exp(-1j*sp_.E*k*0.4)*coef); pr=np.abs(vk)**2; pr/=pr.sum()
        for d in np.unique(rng.choice(n,size=shots,p=pr)): seen.add(int(d))
        fam[("teqsci", k)] = np.array(sorted(seen))
    # (e) exact Krylov-support subspaces: union of supports of H^j phi  (paper's own construction)
    v=sp_.phi.copy(); sup=set()
    for j in range(0,8):
        sup |= set(np.where(np.abs(v)>1e-12)[0].tolist())
        fam[("krylovSup", j)] = np.array(sorted(sup))
        v = sp_.H@v
        if len(sup)==n: break
    return fam

def run(tag, Hs, phi, E0, eta, rng):
    sp_=Spec(Hs,phi,E0,eta)
    log(f"{tag}: dim={sp_.n} ||phi||^2={sp_.nrm2:.6f} normA={sp_.normA:.6f} eta={eta}")
    fam=subspaces(sp_,rng); rows=[]
    for (kind,par),S in fam.items():
        m=sp_.metrics(S); m['kind']=kind; m['par']=par; m['tag']=tag; m['eta']=eta
        m['norm_phi2']=sp_.nrm2; m['dim']=sp_.n
        rows.append(m)
    return sp_, rows

if __name__=='__main__':
    rng=np.random.default_rng(20260918)
    allrows=[]
    for (L,U,pbc,eta) in [(6,8.0,False,0.15),(6,8.0,False,0.05),(6,4.0,False,0.15),(6,8.0,True,0.18)]:
        Hs,phi,E0=hubbard_problem(L=L,U=U,pbc=pbc)
        tag=f"Hub L={L} U={U} {'PBC' if pbc else 'OBC'}"
        _,rows=run(tag,Hs,phi,E0,eta,rng); allrows+=rows
    json.dump(allrows,open(OUT,'w'))
    log("wrote",OUT,len(allrows),"rows")
