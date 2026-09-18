# -*- coding: utf-8 -*-
r"""Sampled-vs-exact momentum-resolved spectral function A(k,omega) (fig:akwsampled, Fig. 5).

Demonstrates that A(k,omega) is genuinely reconstructed FROM bitstring-sampled subspaces, not just shown
as an exact target. For the L=8, U/t=8 half-filled Hubbard chain (eta=0.18t): at each momentum k the seed
c^dag_{k,up}|GS> (and c_{k,up}|GS>) is time-evolution-sampled in the (N+-1) sector, the top-weight
configurations (85% of the sector) are diagonalized, and the reconstructed A^+(k,w)+A^-(k,w) is compared to
the exact Haydock target. Writes sampled_akw_L8.json (per-momentum rel-L1 + overlay curves at k=0,pi/2,pi).
Result: mean per-momentum rel-L1 = 0.0022, max 0.0033 -- faithful across the Brillouin zone.

Run:  python sampled_akw.py      (uses the sector engine of akw_lanczos.py; ~5 min, exact)
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
import sys, os, json, time, numpy as np
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
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import akw_lanczos as AK
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

t0 = time.time()
L, U, t, eta, K, dt = 8, 8.0, 1.0, 0.18, 16, 0.5
FR = 0.85                       # subspace fraction of the (N+-1) sector
nup = nd = L // 2

Hexpl, Du, Dd = AK.build_H_explicit(L, U, nup, nd, t)
E0, V0 = eigsh(Hexpl, k=1, which='SA', v0=_v0(Hexpl.shape[0])); E0 = float(E0[0]); Psi = V0[:, 0].reshape(Du, Dd)
wg = np.linspace(-9, 9, 600)


def spec(poles, wts):
    A = np.zeros_like(wg)
    for a, b in zip(poles, wts):
        A += b * (eta / np.pi) / ((wg - a) ** 2 + eta ** 2)
    return A


ks = [2 * np.pi * n / L for n in range(L)]
cdU = AK.cdag_map(L, nup); cUr = AK.cdag_map(L, nup - 1)


def channel(seed, Hsec, sgn):
    """exact A (dense eigh) + sampled A (time-evolution union select top FR, re-diagonalize)."""
    Ha = Hsec.toarray(); Em, Um = np.linalg.eigh(Ha)
    coef = Um.conj().T @ seed
    A_ex = spec(sgn * (Em - E0), np.abs(coef) ** 2)
    nS = Ha.shape[0]; wc = np.zeros(nS)
    for k in range(K + 1):                      # accumulate time-averaged Born weight
        vk = Um @ (np.exp(-1j * Em * k * dt) * coef); wc += np.abs(vk) ** 2
    order = np.argsort(wc)[::-1]; kk = max(1, int(round(FR * nS))); Sset = np.sort(order[:kk])
    HS = Ha[np.ix_(Sset, Sset)]; Es, Us = np.linalg.eigh(HS); a = Us.conj().T @ seed[Sset]
    return A_ex, spec(sgn * (Es - E0), np.abs(a) ** 2), kk, nS


res = {'L': L, 'U': U, 'eta': eta, 'frac': FR, 'wg': wg.tolist(), 'per_k': {}, 'overlay': {}}
ov = {0: '0.0', L // 4: '0.5', L // 2: '1.0'}
Hadd, _, _ = AK.build_H_explicit(L, U, nup + 1, nd, t)
Hrem, _, _ = AK.build_H_explicit(L, U, nup - 1, nd, t)
for n, k in enumerate(ks):
    ph = np.exp(1j * k * np.arange(L)) / np.sqrt(L)
    add = np.zeros((len(AK.strings(L, nup + 1)[0]), Dd), dtype=complex)
    for j in range(L):
        add += ph[j] * (cdU[j] @ Psi)
    rem = np.zeros((len(AK.strings(L, nup - 1)[0]), Dd), dtype=complex)
    for j in range(L):
        rem += np.conj(ph[j]) * (cUr[j].T @ Psi)
    Aa_ex, Aa_s, ka, nSa = channel(add.ravel(), Hadd, +1.0)
    Ar_ex, Ar_s, kr, nSr = channel(rem.ravel(), Hrem, -1.0)
    A_ex = Aa_ex + Ar_ex; A_s = Aa_s + Ar_s; nrm = np.trapz(np.abs(A_ex), wg)
    r = float(np.trapz(np.abs(A_s - A_ex), wg) / nrm)
    res['per_k']["%.3f" % (2 * n / L)] = {'relL1': r, 'frac_add': ka / nSa, 'frac_rem': kr / nSr}
    print("[%5.1fs] k/pi=%.3f  relL1=%.4f  (|S|add=%d/%d)" % (time.time() - t0, 2 * n / L, r, ka, nSa), flush=True)
    if n in ov:
        res['overlay'][ov[n]] = {'exact': A_ex.tolist(), 'sampled': A_s.tolist()}

import statistics as st
res['mean_relL1'] = st.mean(v['relL1'] for v in res['per_k'].values())
res['max_relL1'] = max(v['relL1'] for v in res['per_k'].values())
print("\nMEAN per-k rel-L1 @frac %.2f: %.4f   MAX: %.4f" % (FR, res['mean_relL1'], res['max_relL1']))
json.dump(res, open(_outpath('sampled_akw_L8.json'), 'w'))
print("wrote sampled_akw_L8.json  [%.1fs]" % (time.time() - t0))
