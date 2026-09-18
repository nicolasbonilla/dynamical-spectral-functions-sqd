# -*- coding: utf-8 -*-
"""Max-rigor checks: first-moment sum rule, peak dispersion, Haydock convergence."""
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
import json, numpy as np, scipy.sparse as sp
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
trapz=np.trapezoid
import akw_lanczos as ak

d=json.load(open(_inpath('akw_lanczos_L12.json'))); wg=np.array(d['wg'])
ks=np.array(sorted(float(x) for x in d['A'].keys()))
M=np.array([d['A'][f'{k:.5f}'] for k in ks])
L=12; U=8.0; eta=0.18

print("=== FIRST-MOMENT sum rule:  int w A(k,w) dw  vs exact  -2t cos(k)  (m=up) ===")
print("    (exact Hubbard f-sum: M1(k)=eps_k-mu+U<n_dn> = -2cos k at half-filling, mu=U/2)")
maxerr=0
for k in ks:
    A=M[int(np.argmin(np.abs(ks-k)))]
    m1=trapz(wg*A,wg); ref=-2*np.cos(k*np.pi); e=m1-ref; maxerr=max(maxerr,abs(e))
    print(f"  k/pi={k:.3f}  M1_data={m1:+.4f}   -2cos={ref:+.4f}   diff={e:+.4f}")
print(f"  max |diff| = {maxerr:.4f}  (small = f-sum rule satisfied; residual from +-9t truncation)")

print("\n=== dominant-peak dispersion from exact data ===")
for k in ks:
    A=M[int(np.argmin(np.abs(ks-k)))]
    neg=wg<0; pos=wg>0
    wr=wg[neg][np.argmax(A[neg])]; wa=wg[pos][np.argmax(A[pos])]
    print(f"  k/pi={k:.3f}  removal-peak={wr:+.2f}t  addition-peak={wa:+.2f}t")

print("\n=== HAYDOCK convergence at k_F=pi/2 (nl=150 used vs nl=300 reference) ===")
nup=nd=L//2
Hexpl,Du,Dd=ak.build_H_explicit(L,U,nup,nd)
E0,V0=eigsh(Hexpl,k=1,which='SA',maxiter=5000,tol=1e-10,v0=_v0(Hexpl.shape[0])); E0=float(E0[0]); psi0=V0[:,0]
Psi=psi0.reshape(Du,Dd)
Ha,_,_,_,_,Dua,_=ak.sector_H(L,U,nup+1,nd)
Hr,_,_,_,_,Dur,_=ak.sector_H(L,U,nup-1,nd)
cdU=ak.cdag_map(L,nup); cU=ak.cdag_map(L,nup-1)
k=np.pi/2; ph=np.exp(1j*k*np.arange(L))/np.sqrt(L)
add=np.zeros((Dua,Dd),dtype=complex)
for j in range(L): add+=ph[j]*(cdU[j]@Psi)
rem=np.zeros((Dur,Dd),dtype=complex)
for j in range(L): rem+=np.conj(ph[j])*(cU[j].T@Psi)
mu=U/2.0; zt=(wg+mu)+1j*eta
def spec(nl):
    Gp=ak.haydock(lambda x:ak._apply(Ha,x)-E0*x, add.ravel(), nl, zt)
    Gm=ak.haydock(lambda x:E0*x-ak._apply(Hr,x), rem.ravel(), nl, zt)
    return -(1/np.pi)*np.imag(Gp+Gm)
A150=spec(150); A300=spec(300)
rel=trapz(np.abs(A150-A300),wg)/trapz(np.abs(A300),wg)
print(f"  rel-L1[A_nl150 - A_nl300] at k_F = {rel:.2e}   (tiny => nl=150 is converged)")
print(f"  sum rule nl=150: {trapz(A150,wg):.4f}   nl=300: {trapz(A300,wg):.4f}")
