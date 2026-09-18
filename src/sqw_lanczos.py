# -*- coding: utf-8 -*-
r"""Dynamical structure factor S(q,omega) of the 1D Hubbard model (PBC), exact, by the
Haydock continued-fraction (Lanczos) method -- the density-density response, a NEUTRAL
(particle-conserving) excitation, computed in the half-filled N-particle sector.

  n_q = sum_j e^{i q j} (n_{j,up}+n_{j,dn})              (q != 0, so <n_q>=0)
  S(q,w) = -(1/pi) Im <psi0| n_{-q} [w-(H-E0)+i eta]^{-1} n_q |psi0>
         = sum_m |<m|n_q|psi0>|^2 delta(w-(E_m-E0))
Same sector engine / Haydock as akw_lanczos (exact spectral functions).
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
        if a[i] in ('--force-eta',):
            i += 1
            continue
        if a[i] in ('--out', '-o', '--in', '-i', '--eta'):
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
import sys
t0=time.time(); log=lambda *a: print(f"[{time.time()-t0:6.1f}s]",*a,flush=True)

def occ_matrix(L,n):
    S,_=strings(L,n)
    return np.array([[ (m>>i)&1 for i in range(L)] for m in S],dtype=float)   # (D,L)

def run(L=12,U=8.0,t=1.0,eta=0.20,nw=600,wwin=(-0.5,13.5),nl=200):
    nup=nd=L//2
    H,Du,Dd=build_H_explicit(L,U,nup,nd,t)
    E0,V0=eigsh(H,k=1,which='SA',v0=_v0(H.shape[0])); E0=float(E0[0]); Psi=V0[:,0].reshape(Du,Dd)
    log(f"L={L} gs dim={Du*Dd} E0={E0:.4f}")
    upocc=occ_matrix(L,nup); dnocc=occ_matrix(L,nd)
    wg=np.linspace(wwin[0],wwin[1],nw); zt=wg+1j*eta
    Happly=lambda x: H@x - E0*x                 # M = H - E0 ; poles at w = E_m-E0
    out={}
    for n in range(1,L):                        # q = 2 pi n / L, skip q=0
        q=2*np.pi*n/L
        ph=np.exp(1j*q*np.arange(L))
        nq_up=upocc@ph; nq_dn=dnocc@ph
        seed=((nq_up[:,None]+nq_dn[None,:])*Psi).ravel()
        G=haydock(Happly, seed, nl, zt)
        S=-(1.0/np.pi)*np.imag(G)
        out[round(q/np.pi,5)]=S.tolist()
        log(f"  q/pi={2*n/L:.3f}  sum S dw={np.trapezoid(S,wg):.3f}")
    return dict(L=L,U=U,eta=eta,E0=E0,wg=wg.tolist(),
                S={f"{qq:.5f}":v for qq,v in out.items()})

if __name__=="__main__":
    _a=_argv_positional(); L=int(_a[0]) if _a else 12
    # --- ETA GUARD (added 2026-09-18) --------------------------------------
    # KNOWN DISCREPANCY, deliberately NOT papered over: the default below is
    # eta=0.20, but data/sqw_L12.json was produced at eta=0.18 and the manuscript
    # quotes 0.18 in three places.  Changing the default would silently change the
    # DEPOSITED datum, so the default stays where the published runs left it and
    # the script says out loud when the run it is about to do cannot reproduce the
    # deposit.  Pass --eta 0.18 to reproduce the deposited file.
    # See docs/KNOWN_DISCREPANCIES.md, entry "sqw_lanczos.py eta default".
    import inspect as _inspect
    _eta_default = float(_inspect.signature(run).parameters['eta'].default)
    _eta = float(_flag('eta')) if _flag('eta') else _eta_default
    _dep = _inpath(f'sqw_L{L}.json')
    _dep_eta = None
    if _os.path.exists(_dep):
        try:    _dep_eta = float(json.load(open(_dep))['eta'])
        except Exception: pass
    # It is not enough to WARN.  A warning scrolls past, the run takes 726 s, and the
    # deposited spectrum is then overwritten at the wrong eta by a reader who did
    # nothing wrong.  So the guard REFUSES: if this run would overwrite a deposited
    # file that was made at a different eta, it stops before computing anything.
    # Escape hatches, both explicit: --out PATH (write somewhere else) or
    # --force-eta (overwrite on purpose).
    _dest = _outpath(f'sqw_L{L}.json')
    _would_clobber = (_os.path.normcase(_os.path.abspath(_dest))
                      == _os.path.normcase(_os.path.abspath(_dep)))
    if _dep_eta is not None and abs(_dep_eta - _eta) > 1e-12:
        _msg = (f"running at eta={_eta} but the DEPOSITED {_os.path.basename(_dep)} "
                f"has eta={_dep_eta}.")
        if _would_clobber and '--force-eta' not in _sys.argv:
            print("=" * 78)
            print(f"REFUSING TO RUN: {_msg}")
            print(f"  This run would spend its whole runtime and then OVERWRITE the deposited")
            print(f"  spectrum at an eta the paper does not quote.  Nothing has been computed.")
            print(f"  To reproduce the deposit:   python src/sqw_lanczos.py {L} --eta {_dep_eta}")
            print(f"  To explore another eta:     python src/sqw_lanczos.py {L} --eta {_eta} "
                  f"--out <path outside data/>")
            print(f"  To overwrite anyway:        add --force-eta")
            print("=" * 78, flush=True)
            _sys.exit(2)
        print("=" * 78)
        print(f"WARNING: {_msg}")
        print(f"         This run will NOT reproduce the deposited file, nor the eta the paper")
        print(f"         quotes. Re-run with:  python src/sqw_lanczos.py {L} --eta {_dep_eta}")
        print(f"         (writing to {_dest}; the deposited file is not touched)"
              if not _would_clobber else
              "         (--force-eta given: the deposited file WILL be overwritten)")
        print("=" * 78, flush=True)
    d=run(L=L,eta=_eta)
    json.dump(d,open(_dest,'w'))
    log(f"WROTE sqw_L{L}.json ({len(d['S'])} q-points, eta={_eta})")
