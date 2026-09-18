# -*- coding: utf-8 -*-
r"""leakage_certificate.py -- a posteriori L1 certificate for subspace-projected spectral functions.

WHAT THIS COMPUTES
------------------
For a Hermitian H on a sector of dimension D, a seed vector phi (unnormalised), a COORDINATE
subspace S (the bitstring set a sampled-subspace method retains) and a broadening eta>0, define

    A  (w) = (eta/pi) || (w+E0+i eta - H  )^-1 phi   ||^2        (exact spectral function)
    A_S(w) = (eta/pi) || (w+E0+i eta - H_S)^-1 phi_S ||^2 ,  H_S = P_S H P_S,  phi_S = P_S phi

(the Rayleigh-Ritz object that every sampled-subspace script in this repository actually builds).
This module evaluates the two-sided bound

    ||phi||^2 (1-w_S)  <=  ||A - A_S||_1  <=  (||phi||+||phi_S||) ( ||phi_Q|| + Lambda_S(eta)/eta )

    relative form:  rel-L1  <=  (1 + sqrt(w_S)) ( sqrt(1-w_S) + Lambda_hat_S(eta)/eta )

with w_S = ||phi_S||^2/||phi||^2, phi_Q = (1-P_S) phi, Lambda_hat_S = Lambda_S/||phi||, and the
eta-averaged leakage

    Lambda_S(eta)^2 = (eta/pi) \int dw || Q H (w+E0+i eta - H_S)^-1 phi_S ||^2
                    = sum_{m,m'} c_m conj(c_m') <g_m', g_m> [4 eta^2 + 2 i eta (Et_m' - Et_m)]
                                                            / [(Et_m' - Et_m)^2 + 4 eta^2]

over the Ritz pairs (Et_m, |m~>) of H_S, with c_m = <m~|phi> and g_m = Q H |m~> (the Ritz RESIDUAL,
g_m = H|m~> - Et_m|m~>, hence only the Q-block H_QS is needed).  Lambda_S is A POSTERIORI and
computable WITHOUT the answer: it needs only H, S and the Ritz pairs the existing code already forms.

WHICH PAPER CLAIM THIS SUPPORTS
-------------------------------
It REPLACES Theorem 1(iii) (paper/theorem_b2.tex, invoked in paper/method.tex), whose published
statement  ||A-A_S||_1 <= 2||phi||^2 (1-w_S) + C rho^{-K_S}  is false (analytic 2x2 counterexample;
186/292 numerical violations; 126/126 in the w_S>=0.99 regime the paper invokes).  The lower bound
also corrects the proof sketch, which quotes ||phi||^2(1-w_S) <= ||A-A_S||_1 with the inequality
reversed.  The quantities w_S and K_S are computed here for the first time in this repository.

COMPLEXITY
----------
Krylov path (default for |S| > dense_max): n_kry Lanczos steps with full reorthogonalisation on the
sparse H_S, one sparse apply of the H_QS block to n_kry vectors, one n_kry x n_kry Gram.  The
resolvent (z-H_S)^-1 phi_S lies exactly in the Krylov space of (H_S, phi_S), so the Krylov path
converges to the exact Lambda_S; certify() reports a convergence residual (kry_conv).
Dense path (|S| <= dense_max): exact eigendecomposition of H_S.

HOW TO RUN
----------
    python src/leakage_certificate.py selftest [out.json]   # self-verification, ~45 s
    (the default out.json is <repo>/data/leakage_certificate_selftest.json)

  As a library:
    from leakage_certificate import certify
    out = certify(H_csr, phi, S_indices, eta, E0=E0, ref_poles=..., ref_weights=...)

SELF-VERIFICATION (aborts on failure)
-------------------------------------
  (a) analytic 2x2 counterexample  H=[[0,g],[g,0]], phi=e0, S={0}:  w_S=1, K_S=0, Lambda_S=g exactly,
      and ||A-A_S||_1 -> 2||phi||^2 as g/eta -> infinity  (the configuration that kills Theorem 1(iii));
  (b) closed-form Lambda_S^2 against direct omega-quadrature on a random dense case;
  (c) the two-sided bound on 40 random dense cases with fixed seeds;
  (d) closed-form Lambda_S^2 against the TIME-DOMAIN identity
      Lambda_S(eta)^2 = 2 eta int_0^inf e^{-2 eta s} ||Q H e^{-i H_S s} phi_S||^2 ds
      (an independent representation that uses no resolvent and no eigen-decomposition).

TIMING:  selftest ~45 s.  A single certificate at |S|=3.3e3 ~ 2 s; at |S|=4.5e4 ~ 2 min.

Repository note: NO container paths.  The default output path is resolved from __file__
(<repo>/data/); a path given on the command line is used as given, relative to the CURRENT
DIRECTORY.  All RNG seeds are fixed and explicit.
"""
from __future__ import annotations

import json
import os
import platform
import sys
import time

import numpy as np
import scipy
import scipy.sparse as sp

# --- NUMPY_TRAPEZOID_BRIDGE ------------------------------------------------
# numpy 2.0 ADDED np.trapezoid and REMOVED np.trapz.  Files in this repository use
# both names, so without this bridge no single numpy version runs the whole deposit.
if not hasattr(np, "trapezoid"):
    np.trapezoid = np.trapz          # numpy < 2.0
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid          # numpy >= 2.0
# ---------------------------------------------------------------------------

__all__ = ["certify", "lambda_leak", "weight_w_S", "krylov_order_K_S", "LineQuad",
           "spectral_from_pairs", "rel_L1", "ritz_pairs", "selftest", "provenance"]

# numpy 1.x / 2.x compatibility for the trapezoidal rule
_TRAPZ = getattr(np, "trapezoid", None) or np.trapz


# ----------------------------------------------------------------------------------
# basic quantities
# ----------------------------------------------------------------------------------
def weight_w_S(phi, S):
    """Retained Born weight w_S = ||P_S phi||^2/||phi||^2, plus ||phi||, ||phi_S||, ||phi_Q||."""
    n2 = float(np.vdot(phi, phi).real)
    nS2 = float(np.vdot(phi[S], phi[S]).real)
    return nS2 / n2, np.sqrt(n2), np.sqrt(nS2), np.sqrt(max(n2 - nS2, 0.0))


def krylov_order_K_S(H, phi, S, kmax=64, tol=0.0):
    """Largest K with supp(H^k phi) contained in S for all k<=K; -1 if supp(phi) is not inside S.

    This is the K_S that Theorem 1(iii) invokes.  For COORDINATE subspaces P_S is diagonal, so
    K_S >= 0 is equivalent to w_S == 1: the two terms of the published bound are mutually
    exclusive.  Computed on the sparsity pattern only (boolean propagation), so it is cheap.
    """
    D = H.shape[0]
    inS = np.zeros(D, dtype=bool)
    inS[np.asarray(S)] = True
    supp = np.abs(np.asarray(phi).ravel()) > tol
    if np.any(supp & ~inS):
        return -1
    pat = H.copy()
    pat.data = np.ones_like(pat.data, dtype=np.float64)
    K = 0
    cur = supp
    for _ in range(kmax):
        nxt = ((pat @ cur.astype(np.float64)) > 0) | cur
        if np.any(nxt & ~inS):
            return K
        if np.array_equal(nxt, cur):
            return kmax          # invariant subspace reached; K_S is unbounded
        cur = nxt
        K += 1
    return kmax


# ----------------------------------------------------------------------------------
# Ritz pairs of H_S (dense or Krylov), and the leakage Lambda_S
# ----------------------------------------------------------------------------------
def _lanczos_full_reorth(matvec, v0, n):
    """Hermitian Lanczos with FULL reorthogonalisation.  Returns (V (m x N), alpha, beta, ||v0||).

    Works for COMPLEX v0 (the A(k,omega) seeds c^dag_{k,up}|GS> are complex for k != 0, pi): the
    basis is carried in complex arithmetic, while alpha and beta stay real because H is Hermitian.
    """
    N = v0.shape[0]
    n = int(min(n, N))
    cplx = np.iscomplexobj(v0)
    V = np.zeros((n, N), dtype=np.complex128 if cplx else np.float64)
    beta0 = float(np.linalg.norm(v0))
    V[0] = v0 / beta0
    alpha = np.zeros(n)
    beta = np.zeros(max(n - 1, 1))
    w = matvec(V[0])
    alpha[0] = float(np.vdot(V[0], w).real)
    w = w - alpha[0] * V[0]
    m = n
    for j in range(1, n):
        for _ in range(2):                      # reorthogonalise twice for numerical safety
            w -= V[:j].T @ (V[:j].conj() @ w)
        b = float(np.linalg.norm(w))
        beta[j - 1] = b
        if b < 1e-11:
            m = j
            break
        V[j] = w / b
        w = matvec(V[j])
        alpha[j] = float(np.vdot(V[j], w).real)
        w = w - alpha[j] * V[j] - b * V[j - 1]
    return V[:m], alpha[:m], beta[:max(m - 1, 0)], beta0


def ritz_pairs(H_S, phi_S, dense_max=900, n_kry=260):
    """Ritz pairs of H_S needed by the certificate.

    Returns (theta, c, Uc, mode): theta the Ritz values, c_m = <m~|phi_S>, Uc the |S| x n matrix of
    Ritz vectors (columns).  In Krylov mode the pairs span the Krylov space of (H_S, phi_S), which
    contains (z-H_S)^-1 phi_S exactly, so nothing relevant to the certificate is lost.
    """
    nS = H_S.shape[0]
    if nS <= dense_max:
        Hd = H_S.toarray() if sp.issparse(H_S) else np.asarray(H_S)
        theta, U = np.linalg.eigh(Hd)
        c = U.conj().T @ phi_S
        return theta, c, U, "dense"
    V, alpha, beta, beta0 = _lanczos_full_reorth((lambda x: H_S @ x), phi_S, n_kry)
    m = len(alpha)
    T = np.diag(alpha) + np.diag(beta[:m - 1], 1) + np.diag(beta[:m - 1], -1)
    theta, Th = np.linalg.eigh(T)
    # <m~|phi_S> = sum_j Th[j,m] <V[j]|phi_S> = beta0 Th[0,m]  (Th is real; T is real symmetric)
    c = beta0 * Th[0, :]
    Uc = V.T @ Th                               # |S| x m Ritz vectors
    return theta, c, Uc, "krylov(%d)" % m


def lambda_leak(H_QS, theta, c, Uc, eta):
    """Lambda_S(eta) in closed form from the Ritz pairs.  H_QS is the |Q| x |S| block of H."""
    R = H_QS @ Uc                               # |Q| x n  residuals g_m = Q H |m~>
    G = R.conj().T @ R                          # n x n Gram <g_m', g_m>
    d = theta[:, None] - theta[None, :]         # d[m',m] = theta_m' - theta_m
    F = (4.0 * eta ** 2 + 2j * eta * d) / (d ** 2 + 4.0 * eta ** 2)
    M = np.outer(np.conj(c), c)                 # M[m',m] = conj(c_m') c_m
    lam2 = complex(np.sum(M * G * F))
    return float(np.sqrt(max(lam2.real, 0.0))), lam2


# ----------------------------------------------------------------------------------
# spectral functions
# ----------------------------------------------------------------------------------
def spectral_from_pairs(poles, weights, grid, eta):
    """A(w) = (eta/pi) sum_m weights_m / ((w - poles_m)^2 + eta^2)."""
    poles = np.asarray(poles, dtype=float)
    weights = np.asarray(weights, dtype=float)
    A = np.zeros_like(grid)
    keep = weights > 1e-16 * max(weights.max() if weights.size else 1.0, 1e-300)
    for a, b in zip(poles[keep], weights[keep]):
        A += b * (eta / np.pi) / ((grid - a) ** 2 + eta ** 2)
    return A


def rel_L1(A, A_ref, grid, denom=None):
    """||A - A_ref||_1 / denom   (denom defaults to ||A_ref||_1 on the same grid)."""
    num = float(_TRAPZ(np.abs(A - A_ref), grid))
    den = float(denom if denom is not None else _TRAPZ(np.abs(A_ref), grid))
    return num / den, num


class LineQuad(object):
    """Quadrature over the WHOLE real line via w = w0 + c tan(u), u uniform on (-pi/2, pi/2).

    The L1 norms in Theorem B are over all of omega.  A Lorentzian tail decays only as 1/w^2, so a
    finite window loses O(eta/W) of the sum rule -- enough to make the (true) lower bound
    ||phi||^2 (1-w_S) <= ||A-A_S||_1 appear violated at the 1e-3 level, and a naive trapezoid on the
    tan-spaced omega points over-counts the tail badly.  Integrating in u with the Jacobian
    c sec^2(u) cures both: the transformed integrand tends to a finite constant at u -> +-pi/2.

    ppw = target number of grid points per eta at the outermost pole.
    """

    __slots__ = ("w", "u", "jac", "n")

    def __init__(self, poles, eta, ppw=24, nmax=4_000_001, nmin=120_001):
        poles = np.asarray(poles, dtype=float)
        w0 = 0.5 * (float(poles.min()) + float(poles.max()))
        c = max(0.5 * (float(poles.max()) - float(poles.min())), 4.0 * eta, 1e-12)
        # d w / d u = c sec^2 u <= 2c at the outermost pole; require (2 c) du < eta/ppw
        n = int(np.clip(np.ceil(2.0 * np.pi * c * ppw / eta), nmin, nmax)) | 1
        eps = 1e-7
        self.u = np.linspace(-0.5 * np.pi + eps, 0.5 * np.pi - eps, n)
        tu = np.tan(self.u)
        self.w = w0 + c * tu
        self.jac = c * (1.0 + tu * tu)
        self.n = n

    def integrate(self, y):
        return float(_TRAPZ(y * self.jac, self.u))


def full_line_grid(poles, eta, ppw=24):
    """Convenience wrapper returning a LineQuad."""
    return LineQuad(poles, eta, ppw=ppw)


# ----------------------------------------------------------------------------------
# the certificate
# ----------------------------------------------------------------------------------
def certify(H, phi, S, eta, E0=0.0, ref_poles=None, ref_weights=None,
            dense_max=900, n_kry=260, grid=None, ppw=24,
            want_K_S=True, kry_check=True, return_pairs=False):
    """Full certificate for the coordinate subspace S.  Returns a JSON-serialisable dict.

    H   : (D,D) scipy.sparse CSR, real symmetric (sector Hamiltonian).
    phi : (D,) seed vector, NOT normalised (||phi||^2 is the sum rule).
    S   : integer index array (the retained bitstrings/configurations).
    eta : Lorentzian broadening.
    E0  : reference energy; the poles of A sit at w = E_m - E0.
    ref_poles / ref_weights : exact Lehmann pairs of H w.r.t. phi, if available.  When given, the
          TRUE rel-L1 and the slack (upper bound / true error) are reported.
    """
    H = H.tocsr() if sp.issparse(H) else sp.csr_matrix(H)
    D = H.shape[0]
    S = np.asarray(np.sort(np.unique(S)), dtype=np.int64)
    nS = int(S.size)
    phi = np.asarray(phi).ravel()

    w_S, nphi, nphiS, nphiQ = weight_w_S(phi, S)

    inS = np.zeros(D, dtype=bool)
    inS[S] = True
    Qidx = np.nonzero(~inS)[0]
    HcolS = H[:, S].tocsr()
    H_S = HcolS[S, :].tocsr()
    H_QS = HcolS[Qidx, :].tocsr()

    theta, c, Uc, mode = ritz_pairs(H_S, phi[S], dense_max=dense_max, n_kry=n_kry)
    lam, lam2 = lambda_leak(H_QS, theta, c, Uc, eta)

    kry_conv = 0.0
    if kry_check and mode.startswith("krylov") and len(theta) > 40:
        theta_h, c_h, Uc_h, _ = ritz_pairs(H_S, phi[S], dense_max=0, n_kry=len(theta) // 2)
        lam_h, _ = lambda_leak(H_QS, theta_h, c_h, Uc_h, eta)
        kry_conv = abs(lam - lam_h) / max(lam, 1e-300)
        del Uc_h

    K_S = krylov_order_K_S(H, phi, S) if want_K_S else None

    lo = nphi ** 2 * (1.0 - w_S)
    up_new = (nphi + nphiS) * (nphiQ + lam / eta)
    up_triv = nphi ** 2 * (1.0 + w_S)
    up = min(up_new, up_triv)

    out = dict(D=int(D), nS=nS, frac=nS / D, eta=float(eta), mode=mode,
               w_S=float(w_S), one_minus_wS=float(1.0 - w_S),
               K_S=(int(K_S) if K_S is not None else None),
               norm_phi=float(nphi), norm_phiS=float(nphiS), norm_phiQ=float(nphiQ),
               Lambda=float(lam), Lambda_hat_over_eta=float(lam / (nphi * eta)),
               Lambda_over_eta=float(lam / eta), kry_conv=float(kry_conv),
               lower=float(lo), upper=float(up), upper_unclipped=float(up_new),
               upper_trivial=float(up_triv),
               rel_lower=float(lo / nphi ** 2), rel_upper=float(up / nphi ** 2))

    if ref_poles is not None:
        if grid is None:
            allp = np.concatenate([np.asarray(ref_poles, dtype=float) - E0, theta - E0])
            quad = LineQuad(allp, eta, ppw=ppw)
            wgrid, integ = quad.w, quad.integrate
        else:
            wgrid = np.asarray(grid, dtype=float)
            integ = (lambda y: float(_TRAPZ(y, wgrid)))
        A_ex = spectral_from_pairs(np.asarray(ref_poles) - E0, ref_weights, wgrid, eta)
        A_S = spectral_from_pairs(theta - E0, np.abs(c) ** 2, wgrid, eta)
        l1_abs = integ(np.abs(A_S - A_ex))
        out["L1_true"] = l1_abs
        out["relL1_true"] = l1_abs / nphi ** 2
        out["relL1_window"] = l1_abs / integ(np.abs(A_ex))
        out["slack"] = out["upper"] / max(l1_abs, 1e-300)
        # QUADRATURE TOLERANCE.  ||A-A_S||_1 is obtained by trapezoidal quadrature on the
        # tangent-transformed full line; its relative error is ~1e-8 at the coarsest settings used
        # here.  The lower bound is EXACTLY SATURATED whenever A - A_S has one sign (large eta,
        # small w_S), so equality cases must not be reported as violations: the test is therefore
        # relative at 1e-6, three orders of magnitude tighter than any genuine violation seen
        # (the falsified Theorem 1(iii) is violated by factors of 9 to 7.6e7, never by 1e-6).
        # On the Krylov path A_S is itself a compression of the Rayleigh-Ritz object (an n_kry-node
        # Haydock quadrature), so its own representation error -- about 1e-6 of the sum rule at the
        # depths used here -- floors how tight any comparison can be.  Ignoring this reports a
        # spurious violation at S = the whole sector, where the bound is exactly 0 while the
        # measured L1 is pure Krylov truncation.
        rtol = 1e-6
        atol = (1e-10 + (5e-6 if mode.startswith("krylov") else 0.0)) * nphi ** 2
        out["lower_ok"] = bool(out["lower"] <= l1_abs * (1 + rtol) + atol)
        out["upper_ok"] = bool(l1_abs <= out["upper"] * (1 + rtol) + atol)
        out["lower_saturated"] = bool(out["lower"] >= l1_abs * (1 - 1e-6) - atol)
        out["bound_nontrivial"] = bool(out["upper_unclipped"] < out["upper_trivial"])
        out["S_is_full_sector"] = bool(nS == D)
        out["window_capture"] = integ(A_ex) / nphi ** 2
    if return_pairs:
        # Ritz poles/weights of A_S; the caller can rebuild A_S on any grid with
        # spectral_from_pairs(theta - E0, ritz_weights, grid, eta) -- no extra diagonalisation.
        out["_ritz_theta"] = theta
        out["_ritz_weight"] = np.abs(c) ** 2
    return out


# ----------------------------------------------------------------------------------
# provenance
# ----------------------------------------------------------------------------------
def provenance(params, controls, t_start, script=None):
    # the name is taken from __file__, not typed: this module was called certificate.py while it
    # lived outside the repository, and a stale hard-coded name is how a JSON ends up pointing at
    # a script that does not exist in the deposit.
    return dict(script=script or os.path.basename(__file__),
                library=os.path.basename(__file__),
                generated_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                python=sys.version.split()[0], numpy=np.__version__, scipy=scipy.__version__,
                platform=platform.platform(), runtime_s=round(time.time() - t_start, 2),
                params=params, controls_reproduced=controls)


# ----------------------------------------------------------------------------------
# MANDATORY self-verification
# ----------------------------------------------------------------------------------
def _twobytwo(g, eta):
    """Analytic 2x2 counterexample: H=[[0,g],[g,0]], phi=e0, S={0}, E0=0."""
    H = sp.csr_matrix(np.array([[0.0, g], [g, 0.0]]))
    phi = np.array([1.0, 0.0])
    S = np.array([0])
    return certify(H, phi, S, eta, E0=0.0,
                   ref_poles=np.array([-abs(g), abs(g)]), ref_weights=np.array([0.5, 0.5]),
                   dense_max=10, ppw=48)


def selftest(verbose=True):
    """Aborts (AssertionError) if any control fails.  Returns a dict of reproduced controls."""
    ok = {}

    # ---- (a) the analytic 2x2 counterexample to Theorem 1(iii) ----------------------
    eta = 0.1
    rows = []
    for g in (0.1, 1.0, 10.0, 100.0, 1000.0):
        o = _twobytwo(g, eta)
        rows.append((g, o["w_S"], o["K_S"], o["Lambda"], o["L1_true"], o["upper"]))
        assert abs(o["w_S"] - 1.0) < 1e-14, "2x2: w_S must be exactly 1"
        assert o["K_S"] == 0, "2x2: K_S must be exactly 0"
        assert abs(o["Lambda"] - abs(g)) < 1e-10 * max(1.0, abs(g)), "2x2: Lambda_S must equal g"
        assert o["upper_ok"] and o["lower_ok"], "2x2: Theorem B violated"
    l1_big = rows[-1][4]
    assert abs(l1_big - 2.0) < 1e-3, \
        "2x2: ||A-A_S||_1 must tend to 2||phi||^2 = 2 (got %.6f)" % l1_big
    ok["2x2_counterexample"] = dict(L1_at_g1000_eta0p1=l1_big, target=2.0,
                                    Lambda_equals_g=True, w_S=1.0, K_S=0,
                                    tol="|L1-2| < 1e-3, |Lambda-g| < 1e-10 rel, |w_S-1| < 1e-14")
    if verbose:
        print("  [a] 2x2 counterexample (Thm 1(iii) killer): w_S=1, K_S=0 exactly")
        for g, w, k, lam, l1, up in rows:
            print("      g=%-8g w_S=%.12f K_S=%d Lambda=%-10.4g L1=%.6f  ThmB upper=%.4g"
                  % (g, w, k, lam, l1, up))
        print("      -> L1 -> 2.000000 = 2||phi||^2 forces C >= 2||phi||^2 in the published bound")

    # ---- (b) closed-form Lambda^2 vs direct omega-quadrature ------------------------
    rng = np.random.default_rng(20260918)
    Dt, nSt, et = 40, 13, 0.13
    Ar = rng.normal(size=(Dt, Dt))
    Ht = (Ar + Ar.T) / 2.0
    pht = rng.normal(size=Dt)
    St = np.sort(rng.choice(Dt, nSt, replace=False))
    Qt = np.setdiff1d(np.arange(Dt), St)
    th, Th = np.linalg.eigh(Ht[np.ix_(St, St)])
    ct = Th.T @ pht[St]
    lam_cf, _ = lambda_leak(sp.csr_matrix(Ht[np.ix_(Qt, St)]), th, ct, Th, et)
    quad = LineQuad(th, et, ppw=200)                # full-line quadrature, Jacobian form
    Rq = Ht[np.ix_(Qt, St)] @ Th
    Gq = Rq.conj().T @ Rq
    val = np.empty(quad.n)
    for i0 in range(0, quad.n, 200000):             # chunked: keeps the footprint small
        sl = slice(i0, min(i0 + 200000, quad.n))
        X = ct[:, None] / ((quad.w[sl] + 1j * et)[None, :] - th[:, None])
        val[sl] = np.einsum('mw,nw,nm->w', X, X.conj(), Gq).real
    lam_q = float(np.sqrt((et / np.pi) * quad.integrate(val)))
    rel = abs(lam_cf - lam_q) / lam_q
    assert rel < 1e-4, "closed-form Lambda disagrees with quadrature (%.3e)" % rel
    ok["closed_form_vs_quadrature"] = dict(Lambda_closed=lam_cf, Lambda_quad=lam_q,
                                           rel_diff=rel, tol="rel < 1e-4")
    if verbose:
        print("  [b] closed-form Lambda = %.9f  quadrature = %.9f  rel diff = %.2e"
              % (lam_cf, lam_q, rel))

    # ---- (c) two-sided bound on 40 random dense cases -------------------------------
    nviol_lo = nviol_up = 0
    ncase = 0
    minslack = np.inf
    for seed in range(40):
        r = np.random.default_rng(5000 + seed)
        D = int(r.integers(12, 46))
        nS = int(r.integers(2, D))
        Ar = r.normal(size=(D, D))
        Hm = sp.csr_matrix((Ar + Ar.T) / 2.0)
        ph = r.normal(size=D)
        Sm = np.sort(r.choice(D, nS, replace=False))
        ev, V = np.linalg.eigh(Hm.toarray())
        et = float(10 ** r.uniform(-1.3, -0.2))
        o = certify(Hm, ph, Sm, et, E0=0.0, ref_poles=ev, ref_weights=np.abs(V.T @ ph) ** 2,
                    dense_max=10 ** 6)
        ncase += 1
        nviol_lo += (not o["lower_ok"])
        nviol_up += (not o["upper_ok"])
        minslack = min(minslack, o["slack"])
    assert nviol_lo == 0 and nviol_up == 0, "Theorem B violated on random controls"
    ok["random_bound_check"] = dict(n=ncase, viol_lower=nviol_lo, viol_upper=nviol_up,
                                    min_slack=minslack)
    if verbose:
        print("  [c] random dense controls: %d cases, 0 violations, min slack = %.4f"
              % (ncase, minslack))

    # ---- (d) closed form vs the TIME-DOMAIN identity (an independent representation) ------
    #     Lambda_S(eta)^2 = 2 eta \int_0^inf e^{-2 eta s} || H_QS e^{-i H_S s} phi_S ||^2 ds
    #     (Laplace representation of 2 i eta / (d + 2 i eta); no resolvent, no eigen-decomposition
    #      of H_S is used on this path, so a convention or sign error would show up here.)
    from scipy.linalg import expm as _expm
    r = np.random.default_rng(160924)
    Dt, nSt, et = 30, 11, 0.25
    Ar = r.normal(size=(Dt, Dt))
    Ht = (Ar + Ar.T) / 2.0
    pht = r.normal(size=Dt) + 1j * r.normal(size=Dt)
    St = np.sort(r.choice(Dt, nSt, replace=False))
    Qt = np.setdiff1d(np.arange(Dt), St)
    HS = Ht[np.ix_(St, St)]
    HQS = Ht[np.ix_(Qt, St)]
    th, Th = np.linalg.eigh(HS)
    lam_cf2, _ = lambda_leak(sp.csr_matrix(HQS), th, Th.conj().T @ pht[St], Th, et)
    ds, smax = 0.002, 60.0
    Ustep = _expm(-1j * ds * HS)
    v = pht[St].copy()
    acc = 0.0
    s = 0.0
    prev = float(np.linalg.norm(HQS @ v) ** 2)
    while s < smax:
        v = Ustep @ v
        s += ds
        cur = float(np.linalg.norm(HQS @ v) ** 2)
        acc += 0.5 * ds * (np.exp(-2 * et * (s - ds)) * prev + np.exp(-2 * et * s) * cur)
        prev = cur
    lam_td = float(np.sqrt(2 * et * acc))
    rel2 = abs(lam_cf2 - lam_td) / lam_td
    assert rel2 < 1e-4, "closed-form Lambda disagrees with the time-domain identity (%.3e)" % rel2
    ok["closed_form_vs_time_domain"] = dict(Lambda_closed=lam_cf2, Lambda_time=lam_td,
                                            rel_diff=rel2, tol="rel < 1e-4")
    if verbose:
        print("  [d] closed-form Lambda = %.9f  time-domain = %.9f  rel diff = %.2e"
              % (lam_cf2, lam_td, rel2))
    return ok


if __name__ == "__main__":
    import os
    _REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if len(sys.argv) > 1 and sys.argv[1] == "selftest":
        t0 = time.time()
        print("leakage_certificate.py SELF-VERIFICATION")
        res = selftest()
        print("ALL CONTROLS PASS  [%.1f s]" % (time.time() - t0))
        outp = (sys.argv[2] if len(sys.argv) > 2
                else os.path.join(_REPO, "data", "leakage_certificate_selftest.json"))
        outp = os.path.abspath(outp)
        os.makedirs(os.path.dirname(outp), exist_ok=True)
        json.dump(dict(selftest=res, provenance=provenance({}, sorted(res), t0)),
                  open(outp, "w"), indent=1, default=float)
        print("wrote", outp)
    else:
        print(__doc__)
