# -*- coding: utf-8 -*-
r"""N2 HERO DATA: structure, spectrum, and resources move as one along the triple-bond dissociation.
For N2 at a grid of bond lengths R (CAS(6,6), cc-pVDZ), compute (all from the validated JW/CAS pipeline,
energy-checked vs FCI at every R):
  * the orbital-projected addition spectral function A(w) for a FIXED frontier (LUMO) orbital -> the single
    quasiparticle peak at equilibrium splits into a strongly-correlated multiplet on stretch;
  * the fermionic magic FAF, the determinant cost n99, the bipartite entanglement S_ent.
Writes n2_hero.json. Reuses n19_suite.py (jw_ops, build_H, faf, ent, n99, get_ops).
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
import n19_suite as NS
from pyscf import gto, scf, mcscf, ao2mo, fci
t0=time.time(); log=lambda *a: print(f"[{time.time()-t0:6.1f}s]",*a,flush=True)
NCAS=6; NELEC=(3,3); na,nb=NELEC; eta=0.06; ng=1000
C,Cd=NS.get_ops(NCAS); M=2*NCAS
Nocc=np.array([bin(s).count('1') for s in range(1<<M)])
Szocc=np.array([sum(((s>>(2*i))&1)-((s>>(2*i+1))&1) for i in range(NCAS)) for s in range(1<<M)])
gi=np.where((Nocc==na+nb)&(Szocc==na-nb))[0]
si=np.where((Nocc==na+nb+1)&(Szocc==(na+1)-nb))[0]
PSEED=NCAS//2                      # fixed frontier (LUMO) orbital, spin up -> consistent A(w) across R

Rgrid=[1.00,1.10,1.25,1.45,1.70,2.00,2.35,2.75]
# common frequency grid (fixed) so the A(w) curves are directly comparable across R
gmin,gmax=-0.2,1.4; grid=np.linspace(gmin,gmax,ng)
def spec(poles,wts):
    A=np.zeros_like(grid)
    for a,b in zip(poles,wts): A+=b*(eta/np.pi)/((grid-a)**2+eta**2)
    return A

out={'ncas':NCAS,'pseed':PSEED,'eta':eta,'grid':grid.tolist(),'points':[]}
for R in Rgrid:
    mol=gto.M(atom=f'N 0 0 0; N 0 0 {R}',basis='cc-pvdz',verbose=0); mf=scf.RHF(mol).run()
    cas=mcscf.CASCI(mf,NCAS,NELEC); cas.ncore=4
    h1,ec=cas.get_h1cas(); h2=ao2mo.restore(1,cas.get_h2cas(),NCAS)
    e_fci,_=fci.direct_spin1.FCI().kernel(h1,h2,NCAS,NELEC,ecore=ec)
    H=NS.build_H(h1,h2,C,Cd,NCAS)
    w,v=np.linalg.eigh(H[gi][:,gi].toarray()); E0=w[0]+ec
    psi0=np.zeros(1<<M,complex); psi0[gi]=v[:,0]
    F=NS.faf(psi0,C,Cd,M); Sent=NS.ent(psi0,M); nd=NS.n99(psi0)
    # gauge-invariant multireference index from natural occupations (eigenvalues of the spin-summed 1-RDM)
    rho=np.zeros((NCAS,NCAS),complex)
    for P in range(NCAS):
        for Q in range(NCAS):
            rho[P,Q]=sum(np.vdot(psi0,Cd[2*P+s]@(C[2*Q+s]@psi0)) for s in (0,1))
    nocc=np.clip(np.linalg.eigvalsh((rho+rho.conj().T).real/2),0,2)
    Nu=float(np.sum(nocc*(2-nocc)))              # effective unpaired electrons (0 = single determinant)
    phi=Cd[2*PSEED]@psi0; phiS=phi[si]; wnorm=float(np.vdot(phiS,phiS).real)
    En,Vn=np.linalg.eigh(H[si][:,si].toarray()); coef=Vn.conj().T@phiS
    poles=(En-w[0]); wts=np.abs(coef)**2               # w[0] is E0 without ecore (ecore cancels in E_n-E_0)
    A=spec(poles,wts)
    npk=int((wts>0.02*wts.max()).sum())                # # of resolved poles (quasiparticle -> multiplet)
    out['points'].append({'R':R,'E0':float(E0),'dFCI':float(abs(E0-e_fci)),'FAF':F,'S_ent':Sent,
                          'n99':nd,'Nu':Nu,'nocc':nocc.tolist(),'sumrule':wnorm,'A':A.tolist(),'npeaks':npk})
    log(f"R={R:.2f}: dFCI={abs(E0-e_fci):.1e} FAF={F:.2f} Nu={Nu:.2f} S={Sent:.2f} peaks~{npk} sum={wnorm:.3f}")
json.dump(out,open(_outpath('n2_hero.json'),'w'))
log("WROTE n2_hero.json")
