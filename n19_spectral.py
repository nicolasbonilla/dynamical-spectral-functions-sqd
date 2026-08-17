# -*- coding: utf-8 -*-
r"""SUITE-WIDE SPECTRAL RECONSTRUCTION — the method (Contribution #1) demonstrated on all 19 molecules.
For each molecule at its dissociation (high-magic) geometry: build the validated CAS Hamiltonian
(JW full-Fock), get the exact ground state; seed the most-addable frontier orbital with c^dagger; compute
the EXACT orbital-projected addition spectral function A(w) (Lehmann in the (N+1) sector) AND the
bitstring-sampled TE-QSCI reconstruction; report the converged relative-L1 error and the sampled-subspace
fraction. Reuses the validated machinery of n19_suite.py / nitrogen_Aw.py. Writes n19_spectral.json.
"""
import json, time, numpy as np
import n19_suite as NS
from pyscf import gto, scf, mcscf, ao2mo, fci
t0=time.time(); log=lambda *a: print(f"[{time.time()-t0:6.1f}s]",*a,flush=True)
rng=np.random.default_rng(0)

def spectral_point(atom,basis,ncore,ncas,nelecas,spin,eta=0.10,K=26,dt=0.5,shots=8000,ng=700):
    na,nb=nelecas
    mol=gto.M(atom=atom,basis=basis,spin=na-nb,verbose=0)
    mf=(scf.RHF(mol) if spin==0 else scf.ROHF(mol)).run()
    cas=mcscf.CASCI(mf,ncas,nelecas); cas.ncore=ncore
    h1,ec=cas.get_h1cas(); h2=ao2mo.restore(1,cas.get_h2cas(),ncas)
    e_fci,_=fci.direct_spin1.FCI().kernel(h1,h2,ncas,nelecas,ecore=ec)
    C,Cd=NS.get_ops(ncas); M=2*ncas
    Nocc=np.array([bin(s).count('1') for s in range(1<<M)])
    Szocc=np.array([sum(((s>>(2*i))&1)-((s>>(2*i+1))&1) for i in range(ncas)) for s in range(1<<M)])
    gi=np.where((Nocc==na+nb)&(Szocc==na-nb))[0]
    H=NS.build_H(h1,h2,C,Cd,ncas)
    w,v=np.linalg.eigh(H[gi][:,gi].toarray()); E0=w[0]+ec
    dFCI=abs(E0-e_fci)
    psi0=np.zeros(1<<M,complex); psi0[gi]=v[:,0]
    # pick the most-addable spin-up orbital (smallest up-occupation) as the frontier seed
    occ=[float((np.abs(C[2*p]@psi0)**2).sum()) for p in range(ncas)]   # <n_{p,up}> = ||c_{p,up}|psi>||^2 ... use c: n=Cd C
    occ=[float(np.vdot(psi0,Cd[2*p]@(C[2*p]@psi0)).real) for p in range(ncas)]
    pseed=int(np.argmin(occ))
    phi=Cd[2*pseed]@psi0
    si=np.where((Nocc==na+nb+1)&(Szocc==(na+1)-nb))[0]; nS=len(si)
    phiS=phi[si]; wnorm=float(np.vdot(phiS,phiS).real)
    if wnorm<1e-10:   # fallback: seed a removable orbital instead (removal branch)
        pseed=int(np.argmax(occ)); phi=C[2*pseed]@psi0
        si=np.where((Nocc==na+nb-1)&(Szocc==(na-1)-nb))[0]; nS=len(si)
        phiS=phi[si]; wnorm=float(np.vdot(phiS,phiS).real); branch='removal'
    else:
        branch='addition'
    En,Vn=np.linalg.eigh(H[si][:,si].toarray()); coef=Vn.conj().T@phiS
    poles=En-E0; wts=np.abs(coef)**2
    grid=np.linspace(poles.min()-0.6,poles.max()+0.6,ng)
    def spec(pw,ww):
        A=np.zeros_like(grid)
        for a,b in zip(pw,ww): A+=b*(eta/np.pi)/((grid-a)**2+eta**2)
        return A
    A_ex=spec(poles,wts); norm=np.trapezoid(A_ex,grid)
    # TE-QSCI sampling: evolve phi in the (N+1) sector, accumulate discovered configs
    seen=set([int(si[int(np.argmax(np.abs(phiS)**2))])]); conv=[]; A_S=A_ex
    for k in range(K+1):
        vk=Vn@(np.exp(-1j*En*k*dt)*coef); pr=np.abs(vk)**2; pr/=pr.sum()
        for dd in np.unique(rng.choice(nS,size=shots,p=pr)): seen.add(int(si[dd]))
        Sset=np.array(sorted(seen)); Em,Um=np.linalg.eigh(H[Sset][:,Sset].toarray())
        aS=Um.conj().T@phi[Sset]; A_S=spec(Em-E0,np.abs(aS)**2)
        l1=float(np.trapezoid(np.abs(A_S-A_ex),grid)/norm)
        conv.append({'k':k,'S':int(len(Sset)),'rel_l1':l1})
    return dict(branch=branch,pseed=pseed,nSector=int(nS),sumrule=wnorm,dFCI=float(dFCI),
                relL1_final=conv[-1]['rel_l1'],frac_final=conv[-1]['S']/nS,
                grid=grid.tolist(),A_exact=A_ex.tolist(),A_sampled=A_S.tolist(),conv=conv)

if __name__=='__main__':
    out={'note':'orbital-projected addition A(w) at dissociation; exact (Lehmann) vs bitstring-sampled TE-QSCI',
         'mols':[]}
    for (lab,bond,basis,ncore,ncas,nelecas,spin,eqA,dsA) in NS.MOLS:
        try:
            r=spectral_point(dsA,basis,ncore,ncas,nelecas,spin)
            r['mol']=lab; r['bonding']=bond; out['mols'].append(r)
            log(f"{lab:5s} {r['branch']:8s} sec={r['nSector']:4d} sum={r['sumrule']:.3f} "
                f"dFCI={r['dFCI']:.1e} | rel-L1={r['relL1_final']:.2e} at frac={r['frac_final']:.2f}")
        except Exception as e:
            log(f"{lab:5s} FAILED: {type(e).__name__}: {e}")
            out['mols'].append(dict(mol=lab,bonding=bond,error=f"{type(e).__name__}: {e}"))
    json.dump(out,open('/w/n19_spectral.json','w'))
    ok=sum(1 for m in out['mols'] if m.get('relL1_final',9)<1e-2)
    log(f"WROTE n19_spectral.json  ({ok}/{len(NS.MOLS)} reach rel-L1<1e-2)")
