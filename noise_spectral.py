# -*- coding: utf-8 -*-
r"""NOISE ROBUSTNESS (real computation for fig:noise). L=6 Hubbard spectral function A(omega) from
bitstrings sampled from the time-evolved seed, passed through a per-qubit BIT-FLIP channel (rate eps),
then reconstructed two ways: (naive) post-select exact Hamming weight/Sz; (S-CoRe) self-consistent
configuration recovery of out-of-sector strings to the (N+1,Sz=+1) sector using mean occupations.
Reports spectral rel-L1 vs exact as a function of eps, seed-averaged. Writes noise_spectral.json."""
import json, time, numpy as np, scipy.sparse as sp
t0=time.time(); log=lambda *a: print(f"[{time.time()-t0:6.1f}s]",*a,flush=True)
L=6; U=4.0; thop=1.0; eta=0.15; M=2*L; dim=1<<M

def c_op(p):
    r=[];c=[];d=[]
    for s in range(dim):
        if (s>>p)&1:
            sg=(-1)**bin(s&((1<<p)-1)).count('1'); r.append(s&~(1<<p)); c.append(s); d.append(float(sg))
    return sp.csr_matrix((d,(r,c)),shape=(dim,dim))
C=[c_op(p) for p in range(M)]; Cd=[c.T.conj() for c in C]
Nocc=np.array([bin(s).count('1') for s in range(dim)])
Szocc=np.array([sum(((s>>(2*i))&1)-((s>>(2*i+1))&1) for i in range(L)) for s in range(dim)])
H=sp.csr_matrix((dim,dim))
for i in range(L):
    j=(i+1)%L
    for spin in (0,1):
        a=2*i+spin; b=2*j+spin; H=H-thop*(Cd[a]@C[b]+Cd[b]@C[a])
for i in range(L): H=H+U*(Cd[2*i]@C[2*i])@(Cd[2*i+1]@C[2*i+1])
H=H.tocsr()
gi=np.where((Nocc==L)&(Szocc==0))[0]
w,v=np.linalg.eigh(H[gi][:,gi].toarray()); E0=w[0]; psi0=np.zeros(dim,complex); psi0[gi]=v[:,0]
phi=Cd[0]@psi0
si=np.where((Nocc==L+1)&(Szocc==1))[0]; nS=len(si); pos={int(x):i for i,x in enumerate(si)}
H7=H[si][:,si].toarray(); En,Vn=np.linalg.eigh(H7); phi7=phi[si]; coef=Vn.conj().T@phi7
grid=np.linspace((En-E0).min()-1,(En-E0).max()+1,600)
def spec(pw,ww):
    A=np.zeros_like(grid)
    for a,b in zip(pw,ww): A+=b*(eta/np.pi)/((grid-a)**2+eta**2)
    return A
A_ex=spec(En-E0,np.abs(coef)**2); norm=np.trapezoid(A_ex,grid)
def relL1_from_configs(cfgset):
    S=np.array(sorted(pos[x] for x in cfgset if x in pos))
    if len(S)<1: return 1.0
    HS=H7[np.ix_(S,S)]; Em,Um=np.linalg.eigh(HS); a=Um.conj().T@phi7[S]
    return float(np.trapezoid(np.abs(spec(Em-E0,np.abs(a)**2)-A_ex),grid)/norm)
log(f"L={L} U={U}: (N+1) sector={nS}")

# ideal sampling distribution over si configs: union weight from time evolution of the seed
K=16; dt=0.5; wcfg=np.zeros(nS)
for k in range(K+1):
    vk=Vn@(np.exp(-1j*En*k*dt)*coef); wcfg+=np.abs(vk)**2
pr=wcfg/wcfg.sum()
def bits_of(cfg): return np.array([(cfg>>p)&1 for p in range(M)],dtype=np.int8)
sicfg=[int(x) for x in si]
# mean occupation per spin-orbital (for recovery), from the ideal distribution
meanocc=np.zeros(M)
for i,cfg in enumerate(sicfg): meanocc+=pr[i]*bits_of(cfg)

def recover(bitrow):
    """S-CoRe: push a noisy bitstring to the (N+1=7, Sz=+1) sector using mean occupations."""
    b=bitrow.copy()
    up=[2*i for i in range(L)]; dn=[2*i+1 for i in range(L)]
    for orbs,target in ((up,4),(dn,3)):   # na=4 up, nb=3 down (N+1=7, Sz=+1)
        occ=[o for o in orbs if b[o]==1]; emp=[o for o in orbs if b[o]==0]
        while len(occ)>target:            # remove electron from lowest-mean-occ set bit
            o=min(occ,key=lambda x:meanocc[x]); b[o]=0; occ.remove(o)
        while len(occ)<target:            # add electron to highest-mean-occ empty bit
            o=max(emp,key=lambda x:meanocc[x]); b[o]=1; emp.remove(o); occ.append(o)
    return int(sum(int(bb)<<p for p,bb in enumerate(b)))

epslist=[0.0,0.01,0.02,0.04,0.06,0.08,0.10,0.13,0.16,0.20]; SEEDS=range(6); shots=6000
out={'L':L,'U':U,'nSector':nS,'eps':epslist,'naive':[],'naive_std':[],'score':[],'score_std':[]}
for eps in epslist:
    rn=[]; rs=[]
    for sd in SEEDS:
        rng=np.random.default_rng(sd)
        draws=rng.choice(nS,size=shots,p=pr)
        bits=np.array([bits_of(sicfg[i]) for i in draws],dtype=np.int8)
        flips=rng.random(bits.shape)<eps; bits=bits^flips.astype(np.int8)
        naive=set(); score=set()
        for row in bits:
            cfg=int(sum(int(bb)<<p for p,bb in enumerate(row)))
            if Nocc[cfg]==L+1 and Szocc[cfg]==1: naive.add(cfg); score.add(cfg)
            else: score.add(recover(row))
        rn.append(relL1_from_configs(naive)); rs.append(relL1_from_configs(score))
    out['naive'].append(float(np.mean(rn))); out['naive_std'].append(float(np.std(rn)))
    out['score'].append(float(np.mean(rs))); out['score_std'].append(float(np.std(rs)))
    log(f"eps={eps:.2f}: naive rel-L1={np.mean(rn):.3f}  S-CoRe rel-L1={np.mean(rs):.3f}")
json.dump(out,open('/w/noise_spectral.json','w'),indent=1)
log("WROTE noise_spectral.json")
