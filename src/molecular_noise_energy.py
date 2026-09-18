# -*- coding: utf-8 -*-
r"""'Why only N2?' -> the noise-assisted SQD ground-state energy ACROSS the molecular suite (simulated,
realistic bit-flip noise + S-CoRe recovery). For each molecule: exact GS -> Born-sample its determinant
distribution -> per-qubit bit-flip channel (rate eps, realistic Heron ~2%) -> S-CoRe recovery to the correct
(na,nb) sector via mean occupations -> diagonalize H in the recovered subspace -> energy error vs FCI. Seed-
averaged. Reports naive post-selection vs S-CoRe recovery across the suite. Writes molecular_noise.json.
"""
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
import json, time, numpy as np
from pyscf import gto, scf, mcscf, ao2mo, fci
import n19_suite as N
t0=time.time(); log=lambda *a: print(f"[{time.time()-t0:6.1f}s]",*a,flush=True)

def gs_pack(atom,basis,ncore,ncas,nelecas,spin):
    na,nb=nelecas
    mol=gto.M(atom=atom,basis=basis,spin=(na-nb),verbose=0)
    mf=(scf.RHF(mol) if spin==0 else scf.ROHF(mol)).run()
    cas=mcscf.CASCI(mf,ncas,nelecas); cas.ncore=ncore
    h1,ec=cas.get_h1cas(); h2=ao2mo.restore(1,cas.get_h2cas(),ncas)
    e_ref,_=fci.direct_spin1.FCI().kernel(h1,h2,ncas,nelecas,ecore=ec)
    C,Cd=N.get_ops(ncas); M=2*ncas; dim=1<<M
    Nocc=np.array([bin(s).count('1') for s in range(dim)])
    Szocc=np.array([sum(((s>>(2*i))&1)-((s>>(2*i+1))&1) for i in range(ncas)) for s in range(dim)])
    gi=np.where((Nocc==na+nb)&(Szocc==na-nb))[0]
    H=N.build_H(h1,h2,C,Cd,ncas).tocsr()
    Hs=H[gi][:,gi].toarray(); w,v=np.linalg.eigh(Hs)
    E0=w[0]+ec; psi0=np.zeros(dim,complex); psi0[gi]=v[:,0]
    return dict(M=M,dim=dim,Nocc=Nocc,Szocc=Szocc,gi=gi,H=H,ec=ec,E0=E0,E_FCI=e_ref,
                psi0=psi0,na=na,nb=nb,ncas=ncas)

def energy_in_subspace(P,cfgset,ec):
    S=np.array(sorted(cfgset))
    if len(S)<1: return np.inf
    HS=P['H'][S][:,S].toarray(); w=np.linalg.eigvalsh(HS); return float(w[0]+ec)

def run_mol(P,eps=0.02,seeds=5,shots=4000):
    M=P['M']; gi=P['gi']; pr=np.abs(P['psi0'][gi])**2; pr/=pr.sum()
    cfgs=[int(x) for x in gi]
    bits_of=lambda cfg: np.array([(cfg>>p)&1 for p in range(M)],dtype=np.int8)
    meanocc=np.zeros(M)
    for i,cfg in enumerate(cfgs): meanocc+=pr[i]*bits_of(cfg)
    up=[2*i for i in range(P['ncas'])]; dn=[2*i+1 for i in range(P['ncas'])]
    def recover(b):
        b=b.copy()
        for orbs,tgt in ((up,P['na']),(dn,P['nb'])):
            occ=[o for o in orbs if b[o]==1]; emp=[o for o in orbs if b[o]==0]
            while len(occ)>tgt: o=min(occ,key=lambda x:meanocc[x]); b[o]=0; occ.remove(o)
            while len(occ)<tgt: o=max(emp,key=lambda x:meanocc[x]); b[o]=1; emp.remove(o); occ.append(o)
        return int(sum(int(bb)<<p for p,bb in enumerate(b)))
    en=[];es=[]
    for sd in seeds if hasattr(seeds,'__iter__') else range(seeds):
        rng=np.random.default_rng(sd)
        draws=rng.choice(len(cfgs),size=shots,p=pr)
        bits=np.array([bits_of(cfgs[i]) for i in draws],dtype=np.int8)
        bits=bits^(rng.random(bits.shape)<eps).astype(np.int8)
        naive=set(); score=set()
        for row in bits:
            cfg=int(sum(int(bb)<<p for p,bb in enumerate(row)))
            if P['Nocc'][cfg]==P['na']+P['nb'] and P['Szocc'][cfg]==P['na']-P['nb']:
                naive.add(cfg); score.add(cfg)
            else: score.add(recover(row))
        en.append(1e3*(energy_in_subspace(P,naive,P['ec'])-P['E_FCI']))   # mHa error
        es.append(1e3*(energy_in_subspace(P,score,P['ec'])-P['E_FCI']))
    return np.mean(en),np.std(en),np.mean(es),np.std(es)

if __name__=='__main__':
    EPS=0.02   # realistic Heron per-qubit bit-flip
    out={'eps':EPS,'note':'simulated realistic-noise SQD ground-state energy + S-CoRe, across the suite','mols':[]}
    print(f"{'mol':6}{'CAS':8}{'naive(mHa)':>16}{'S-CoRe(mHa)':>16}")
    for (lab,bond,basis,ncore,ncas,nelecas,spin,eqA,dsA) in N.MOLS:
        try:
            P=gs_pack(dsA,basis,ncore,ncas,nelecas,spin)   # dissociation (harder, multireference)
            mn,sn,ms,ss=run_mol(P,eps=EPS)
            out['mols'].append(dict(mol=lab,ncas=ncas,nelec=list(nelecas),naive=mn,naive_std=sn,score=ms,score_std=ss))
            chem="  <-- chemical accuracy" if ms<1.6 else ""
            print(f"{lab:6}CAS({sum(nelecas)},{ncas}){mn:8.2f}+-{sn:4.2f}   {ms:8.2f}+-{ss:4.2f}{chem}",flush=True)
        except Exception as e:
            print(f"{lab:6} FAILED: {type(e).__name__}: {e}",flush=True)
    good=[m for m in out['mols'] if m.get('score',9)<1.6]
    print(f"\n==> S-CoRe reaches chemical accuracy (<1.6 mHa) for {len(good)}/{len(out['mols'])} dissociated molecules at eps={EPS}")
    json.dump(out,open(_outpath('molecular_noise.json'),'w'),indent=1)
