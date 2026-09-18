# -*- coding: utf-8 -*-
r"""Noise-robustness of SQD ground-state energy ACROSS chemistry, as a function of the per-qubit bit-flip
rate. For a representative set spanning easy->hard multireference (N2, H2O, NH3, H6, CO, C2), sweep eps and
compare naive post-selection vs S-CoRe recovery. Shows: (i) SQD energy is intrinsically noise-robust because
subspace diagonalization is variational; (ii) S-CoRe extends chemical accuracy to noise rates where naive
post-selection collapses (too few surviving in-sector configurations). Writes molecular_noise_sweep.json.
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
import n19_suite as N
from molecular_noise_energy import gs_pack, energy_in_subspace
t0=time.time(); log=lambda *a: print(f"[{time.time()-t0:6.1f}s]",*a,flush=True)

PICK=["CO","N2","C2","H2O","NH3","H6"]          # representative, spanning multireference strength
EPS=[0.0,0.02,0.05,0.10,0.15,0.20,0.30]
SEEDS=range(8); SHOTS=4000

def sweep_mol(P):
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
    rows=[]
    for eps in EPS:
        en=[];es=[];kept=[]
        for sd in SEEDS:
            rng=np.random.default_rng(1000*sd+int(eps*100))
            draws=rng.choice(len(cfgs),size=SHOTS,p=pr)
            bits=np.array([bits_of(cfgs[i]) for i in draws],dtype=np.int8)
            if eps>0: bits=bits^(rng.random(bits.shape)<eps).astype(np.int8)
            naive=set(); score=set()
            for row in bits:
                cfg=int(sum(int(bb)<<p for p,bb in enumerate(row)))
                if P['Nocc'][cfg]==P['na']+P['nb'] and P['Szocc'][cfg]==P['na']-P['nb']:
                    naive.add(cfg); score.add(cfg)
                else: score.add(recover(row))
            en.append(1e3*(energy_in_subspace(P,naive,P['ec'])-P['E_FCI']))
            es.append(1e3*(energy_in_subspace(P,score,P['ec'])-P['E_FCI']))
            kept.append(len(naive))
        rows.append(dict(eps=eps,naive=float(np.mean(en)),naive_std=float(np.std(en)),
                         score=float(np.mean(es)),score_std=float(np.std(es)),
                         kept_naive=float(np.mean(kept))))
    return rows

if __name__=='__main__':
    MB={m[0]:m for m in N.MOLS}
    out={'seeds':len(list(SEEDS)),'shots':SHOTS,'eps_grid':EPS,'mols':[]}
    for lab in PICK:
        (l,bond,basis,ncore,ncas,nelecas,spin,eqA,dsA)=MB[lab]
        P=gs_pack(dsA,basis,ncore,ncas,nelecas,spin)
        rows=sweep_mol(P)
        out['mols'].append(dict(mol=lab,ncas=ncas,rows=rows))
        log(f"{lab}: naive@0.10={rows[3]['naive']:.2f}  S-CoRe@0.10={rows[3]['score']:.2f}  "
            f"naive@0.20={rows[5]['naive']:.2f}  S-CoRe@0.20={rows[5]['score']:.2f} mHa")
    json.dump(out,open(_outpath('molecular_noise_sweep.json'),'w'),indent=1)
    log("WROTE molecular_noise_sweep.json")
