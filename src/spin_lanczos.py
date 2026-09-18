# -*- coding: utf-8 -*-
r"""Dynamical SPIN structure factor S^{zz}(q,omega) of the 1D Hubbard model (PBC), exact, by the
Haydock continued-fraction (Lanczos) method -- the spin-spin response, a NEUTRAL, spin-conserving
excitation in the half-filled N-particle sector. Companion of the charge response S(q,omega).

  Sz_q = sum_j e^{i q j} (n_{j,up} - n_{j,dn})/2          (q != 0, so <Sz_q>=0)
  Szz(q,w) = -(1/pi) Im <psi0| Sz_{-q} [w-(H-E0)+i eta]^{-1} Sz_q |psi0>
The spin sector -> Heisenberg with J=4t^2/U at strong coupling: gapless 2-spinon continuum,
lower edge (des Cloizeaux-Pearson) (pi J/2)|sin q|, upper edge pi J |sin(q/2)|.
"""
# --- DEPOSIT PATHS (repaired 2026-09-18, second pass) -----------------------
# Eighteen scripts were repaired earlier today because they hard-coded their output
# under '/w/', the working directory of the Docker container the published runs were
# made in.  THIS FILE WAS NOT AMONG THEM, and it was broken in a quieter way: it wrote
# to the RELATIVE path 'data/...', which lands in whatever directory the reader happens
# to be standing in, and raises FileNotFoundError from anywhere except the repository
# root.  The Data Availability Statement claims that every deposited computation script
# writes into the repository's own data/ directory; that sentence was FALSE for this
# file until now.  Repaired exactly like the other eighteen: the default is computed
# from THIS FILE's location, and --out overrides it.
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


def _outpath(name, sub='data'):
    """Absolute path to write `name` to: --out if given, else <repo>/<sub>/<name>."""
    p = _os.path.abspath(_flag('out') or _os.path.join(_REPO, sub, name))
    _os.makedirs(_os.path.dirname(p), exist_ok=True)
    return p


def _repopath(name, sub='data'):
    """Absolute path <repo>/<sub>/<name>.  Never overridden: for a script's SECOND
    output, which --out (a single flag) cannot address unambiguously."""
    p = _os.path.abspath(_os.path.join(_REPO, sub, name))
    _os.makedirs(_os.path.dirname(p), exist_ok=True)
    return p
# ---------------------------------------------------------------------------
import json, time, numpy as np
# --- NUMPY_TRAPEZOID_BRIDGE ------------------------------------------------
# numpy 2.0 ADDED np.trapezoid and REMOVED np.trapz.  Files in this repository use
# both names, so without this bridge no single numpy version runs the whole deposit:
# numpy 1.x breaks the files that call trapezoid, numpy 2.x breaks the files that call
# trapz (this guardian included).  requirements.txt asks for numpy>=1.24; with the
# bridge that is true again.
if not hasattr(np, "trapezoid"):
    np.trapezoid = np.trapz          # numpy < 2.0
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid          # numpy >= 2.0
# ---------------------------------------------------------------------------
from scipy.sparse.linalg import eigsh
# --- DETERMINISTIC ARPACK START VECTOR (added 2026-09-18, second pass) ------
# eigsh() with no v0= lets ARPACK draw its own start vector from an UNSEEDED
# generator, so E0 converges to a slightly different point on every run (the last
# few digits move) and every quantity derived from it moves with it.  Measured, not
# hypothetical: three consecutive calls on the same matrix gave
# -2.0481308860914536 / ...504 / ...522, and in the leakage certificate a small
# subspace selection moved by 3.2%.  The eigenpair is the same to ARPACK's
# tolerance -- no physics changes -- but a deposit must be bit-reproducible.
# Deliberately NOT a numpy global seed: this touches only the ARPACK start vector.
# The seed 20260918 is the one used by scaling_lanczos.py and the certificate suite.
def _v0(n):
    """Fixed, dimension-dependent ARPACK start vector (never orthogonal to the GS)."""
    return np.random.default_rng(20260918).standard_normal(n)
# ---------------------------------------------------------------------------
from akw_lanczos import build_H_explicit, strings, haydock
t0=time.time(); log=lambda *a: print(f"[{time.time()-t0:6.1f}s]",*a,flush=True)
def occ_matrix(L,n):
    S,_=strings(L,n); return np.array([[(m>>i)&1 for i in range(L)] for m in S],dtype=float)
L=12;U=8.0;t=1.0;eta=0.10;nw=600;nl=200; nup=nd=L//2
H,Du,Dd=build_H_explicit(L,U,nup,nd,t)
E0,V0=eigsh(H,k=1,which='SA',v0=_v0(H.shape[0])); E0=float(E0[0]); Psi=V0[:,0].reshape(Du,Dd)
log(f"E0={E0:.6f}")
upocc=occ_matrix(L,nup); dnocc=occ_matrix(L,nd)
wg=np.linspace(-0.1,3.5,nw); zt=wg+1j*eta
Happly=lambda x: H@x - E0*x
out={}
for n in range(1,L):
    q=2*np.pi*n/L; ph=np.exp(1j*q*np.arange(L))
    nq_up=upocc@ph; nq_dn=dnocc@ph
    seed=(((nq_up[:,None]-nq_dn[None,:])/2.0)*Psi).ravel()   # Sz_q |psi0>
    Szz=-(1.0/np.pi)*np.imag(haydock(Happly,seed,nl,zt))
    out[round(q/np.pi,5)]=Szz.tolist()
    log(f"  q/pi={2*n/L:.3f} sum={np.trapz(Szz,wg):.3f}")
d=dict(L=L,U=U,eta=eta,E0=E0,J=4*t*t/U,wg=wg.tolist(),S={f"{qq:.5f}":v for qq,v in out.items()})
json.dump(d,open(_outpath('spinqw_L12.json'),'w'))
log("WROTE spinqw_L12.json (eta=0.10)")
# lower/upper 2-spinon edges from the data (3% threshold) for the figure
qs=np.array(sorted(out.keys())); M=np.array([out[q] for q in qs]); wgn=np.array(wg)
lines=["qpi lo hi peak"]
for q in qs:
    A=np.array(out[q]); thr=0.03*M.max(); sig=wgn[A>thr]
    lo=sig.min() if len(sig) else 0.0; hi=sig.max() if len(sig) else 0.0; pk=wgn[np.argmax(A)]
    lines.append(f"{q:.4f} {lo:.3f} {hi:.3f} {pk:.3f}")
open(_repopath('spinqw_edges.dat', sub='paper/figs'),'w').write("\n".join(lines)+"\n")
log("WROTE spinqw_edges.dat"); print("EDGES:\n"+"\n".join(lines))
