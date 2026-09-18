# -*- coding: utf-8 -*-
r"""THE NINETEEN-MOLECULE FERMIONIC-MAGIC SUITE — real, reproducible, self-verified.
Generalizes the validated N2 pipeline (nitrogen_spectral.py) to arbitrary spin and active space.
For each molecule at EQUILIBRIUM and at DISSOCIATION:
  * RHF/ROHF -> CASCI valence active space (the bond(s) being broken) -> h1,h2, active-space FCI.
  * Build the Hamiltonian in the Jordan-Wigner full-Fock basis; diagonalize in the (N,Sz) sector;
    VERIFY E0 == E_FCI (proves the JW basis is correct so the magic is meaningful).
  * Fermionic AntiFlatness FAF (k=1, Majorana covariance; =0 iff Gaussian/mean-field),
    bipartite entanglement S_ent, determinant cost n99. Degeneracy gap flagged.
Nothing is claimed that a verified run does not show. Writes n19_suite.json.
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
import scipy.sparse as sp
from pyscf import gto, scf, mcscf, ao2mo, fci
t0=time.time(); log=lambda *a: print(f"[{time.time()-t0:6.1f}s]",*a,flush=True)

# ---- JW operators / Hamiltonian / resource measures (arbitrary NCAS) ----
def jw_ops(M):
    dim=1<<M
    def c_op(p):
        r=[];c=[];d=[]
        for s in range(dim):
            if (s>>p)&1:
                sg=(-1)**bin(s&((1<<p)-1)).count('1'); r.append(s&~(1<<p)); c.append(s); d.append(float(sg))
        return sp.csr_matrix((d,(r,c)),shape=(dim,dim))
    C=[c_op(p) for p in range(M)]; return C,[c.T.conj() for c in C]

def build_H(h1,h2,C,Cd,ncas):
    M=2*ncas; dim=1<<M; H=sp.csr_matrix((dim,dim))
    for P in range(ncas):
        for Q in range(ncas):
            if abs(h1[P,Q])>1e-12:
                for s in (0,1): H=H+h1[P,Q]*(Cd[2*P+s]@C[2*Q+s])
    for P in range(ncas):
     for Q in range(ncas):
      for R in range(ncas):
       for S in range(ncas):
        v=h2[P,Q,R,S]
        if abs(v)>1e-12:
            for s in (0,1):
                for s2 in (0,1):
                    H=H+0.5*v*(Cd[2*P+s]@Cd[2*R+s2]@C[2*S+s2]@C[2*Q+s])
    return H.tocsr()

def faf(psi,C,Cd,M):
    G=[]
    for p in range(M): G.append(C[p]@psi+Cd[p]@psi); G.append(-1j*(C[p]@psi-Cd[p]@psi))
    G=np.array(G); Mc=np.zeros((2*M,2*M))
    for a in range(2*M):
        for b in range(a+1,2*M):
            v=-1j*np.vdot(G[a],G[b]); Mc[a,b]=v.real; Mc[b,a]=-v.real
    sv=np.linalg.svd(Mc,compute_uv=False); return float(M-0.5*np.sum(sv**2))
def ent(full,M):
    A=full.reshape(1<<(M//2),1<<(M//2)); s=np.linalg.svd(A,compute_uv=False)
    p=s[s>1e-12]**2; p/=p.sum(); return float(-(p*np.log(p)).sum())
def chi_eps(full,M,eps=0.05):
    """eps-truncated Schmidt rank (= minimal MPS bond dimension) across the central spin-orbital cut."""
    A=full.reshape(1<<(M//2),1<<(M//2)); s=np.linalg.svd(A,compute_uv=False)
    p=np.sort((s[s>1e-12]**2))[::-1]; p/=p.sum(); tail=np.cumsum(p[::-1])[::-1]
    for r in range(len(p)+1):
        if (tail[r] if r<len(tail) else 0.0)<=eps*eps: return r
    return len(p)
def n99(psi,thr=0.99):
    w=np.sort(np.abs(psi)**2)[::-1]; c=np.cumsum(w)/w.sum(); return int(np.searchsorted(c,thr)+1)

# ---- the suite: (label, bonding, basis, ncore, ncas, nelecas, spin=2S, eq-atoms, diss-atoms) ----
# geometries in Angstrom; equilibrium near experimental, dissociation ~2-2.5x the broken bond(s).
def diat(a,b,R): return f"{a} 0 0 0; {b} 0 0 {R}"
MOLS=[
 ("H2","covalent single","cc-pvdz",0,2,(1,1),0, diat("H","H",0.741), diat("H","H",2.0)),
 ("LiH","metal hydride","cc-pvdz",1,2,(1,1),0, diat("Li","H",1.596), diat("Li","H",3.2)),
 ("HF","polar single","cc-pvdz",4,2,(1,1),0, diat("H","F",0.917), diat("H","F",2.3)),
 ("F2","covalent single","cc-pvdz",8,2,(1,1),0, diat("F","F",1.412), diat("F","F",3.0)),
 ("LiF","ionic","cc-pvdz",5,2,(1,1),0, diat("Li","F",1.564), diat("Li","F",3.2)),
 ("BH","metal hydride","cc-pvdz",1,4,(2,2),0, diat("B","H",1.232), diat("B","H",2.6)),
 ("CO","polar triple","cc-pvdz",4,6,(3,3),0, diat("C","O",1.128), diat("C","O",2.4)),
 ("N2","covalent triple","cc-pvdz",4,6,(3,3),0, diat("N","N",1.098), diat("N","N",2.4)),
 ("C2","covalent double","cc-pvdz",2,6,(4,4),0, diat("C","C",1.243), diat("C","C",2.6)),
 ("BeH2","multicentre","cc-pvdz",1,4,(2,2),0, "Be 0 0 0; H 0 0 1.334; H 0 0 -1.334",
                                              "Be 0 0 0; H 0 0 2.9; H 0 0 -2.9"),
 ("H2O","multicentre","cc-pvdz",3,4,(2,2),0, "O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587",
                                             "O 0 0 0; H 0 1.817 1.409; H 0 -1.817 1.409"),
 ("NH3","multicentre","cc-pvdz",2,6,(3,3),0,
   "N 0 0 0; H 0 -0.9377 -0.3816; H 0.8121 0.4689 -0.3816; H -0.8121 0.4689 -0.3816",
   "N 0 0 0; H 0 -2.157 -0.878; H 1.868 1.079 -0.878; H -1.868 1.079 -0.878"),
 ("H4","H-chain","cc-pvdz",0,4,(2,2),0, "H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7",
                                        "H 0 0 0; H 0 0 2.0; H 0 0 4.0; H 0 0 6.0"),
 ("H6","H-chain","cc-pvdz",0,6,(3,3),0,
   "H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7; H 0 0 3.6; H 0 0 4.5",
   "H 0 0 0; H 0 0 2.0; H 0 0 4.0; H 0 0 6.0; H 0 0 8.0; H 0 0 10.0"),
 # --- open-shell / radicals (ROHF; degeneracy flagged) ---
 ("BeH","metal hydride (radical)","cc-pvdz",1,3,(2,1),1, diat("Be","H",1.343), diat("Be","H",2.7)),
 ("OH","polar (radical)","cc-pvdz",2,4,(3,2),1, diat("O","H",0.970), diat("O","H",2.4)),
 ("O2","covalent double (triplet)","cc-pvdz",2,6,(5,3),2, diat("O","O",1.208), diat("O","O",2.6)),
 ("CN","polar (radical)","cc-pvdz",3,6,(4,3),1, diat("C","N",1.172), diat("C","N",2.5)),
 ("NO","polar (radical)","cc-pvdz",4,6,(4,3),1, diat("N","O",1.151), diat("N","O",2.5)),
]

# cache JW ops by ncas
_cache={}
def get_ops(ncas):
    if ncas not in _cache: _cache[ncas]=jw_ops(2*ncas)
    return _cache[ncas]

def one_point(atom,basis,ncore,ncas,nelecas,spin):
    na,nb=nelecas
    mol=gto.M(atom=atom,basis=basis,spin=(na-nb),verbose=0)
    mf=(scf.RHF(mol) if spin==0 else scf.ROHF(mol)).run()
    cas=mcscf.CASCI(mf,ncas,nelecas); cas.ncore=ncore
    h1,ec=cas.get_h1cas(); h2=ao2mo.restore(1,cas.get_h2cas(),ncas)
    e_ref,_=fci.direct_spin1.FCI().kernel(h1,h2,ncas,nelecas,ecore=ec)
    C,Cd=get_ops(ncas); M=2*ncas
    Nocc=np.array([bin(s).count('1') for s in range(1<<M)])
    Szocc=np.array([sum(((s>>(2*i))&1)-((s>>(2*i+1))&1) for i in range(ncas)) for s in range(1<<M)])
    gi=np.where((Nocc==na+nb)&(Szocc==na-nb))[0]
    H=build_H(h1,h2,C,Cd,ncas)
    Hs=H[gi][:,gi].toarray()
    w,v=np.linalg.eigh(Hs); E0=w[0]+ec; gap=float(w[1]-w[0]) if len(w)>1 else float('inf')
    psi0=np.zeros(1<<M,complex); psi0[gi]=v[:,0]
    return dict(E0=float(E0),E_FCI=float(e_ref),dFCI=float(abs(E0-e_ref)),
                FAF=faf(psi0,C,Cd,M),S_ent=ent(psi0,M),n99=n99(psi0),chi=chi_eps(psi0,M),
                gap=gap,secdim=int(len(gi)),fafmax=float(M))

def run_all():
 out={"note":"real self-verified suite; FAF k=1 Majorana antiflatness, range [0,2*ncas]","mols":[]}
 for (lab,bond,basis,ncore,ncas,nelecas,spin,eqA,dsA) in MOLS:
    try:
        eq=one_point(eqA,basis,ncore,ncas,nelecas,spin)
        ds=one_point(dsA,basis,ncore,ncas,nelecas,spin)
        okE = eq['dFCI']<1e-5 and ds['dFCI']<1e-5
        rec=dict(mol=lab,bonding=bond,ncas=ncas,nelec=list(nelecas),spin=spin,
                 FAF_eq=eq['FAF'],FAF_diss=ds['FAF'],fafmax=eq['fafmax'],
                 S_eq=eq['S_ent'],S_diss=ds['S_ent'],n99_eq=eq['n99'],n99_diss=ds['n99'],
                 chi_eq=eq['chi'],chi_diss=ds['chi'],
                 dFCI_eq=eq['dFCI'],dFCI_diss=ds['dFCI'],gap_eq=eq['gap'],gap_diss=ds['gap'],
                 secdim=eq['secdim'],energyOK=bool(okE))
        out['mols'].append(rec)
        flag="" if okE else "  <<< ENERGY MISMATCH"
        dgn="" if min(eq['gap'],ds['gap'])>1e-4 else "  [near-degenerate GS]"
        log(f"{lab:5s} CAS({sum(nelecas)}e,{ncas}o) FAF {eq['FAF']:5.3f}->{ds['FAF']:5.3f} "
            f"(max {eq['fafmax']:.0f}) | dFCI {eq['dFCI']:.1e}/{ds['dFCI']:.1e}{flag}{dgn}")
    except Exception as e:
        log(f"{lab:5s} FAILED: {type(e).__name__}: {e}")
        out['mols'].append(dict(mol=lab,bonding=bond,error=f"{type(e).__name__}: {e}"))
 json.dump(out,open(_outpath('n19_suite.json'),'w'),indent=1)
 log(f"WROTE n19_suite.json  ({sum(1 for m in out['mols'] if m.get('energyOK'))}/{len(MOLS)} energy-verified)")
 return out

if __name__=='__main__':
    run_all()
