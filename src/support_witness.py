# -*- coding: utf-8 -*-
r"""support_witness.py -- the TIE-FREE witness subspace  S* = supp(phi),  and what it measures.

WHAT THIS COMPUTES
------------------
For the half-filled 1D Hubbard ring (PBC, U/t = 8) at L = 6 and L = 8, with the seed of the
paper's addition channel

    phi = c^dag_{0,up} |Psi_0>          in the (N+1, Sz = +1) sector,

this script takes the ONE subspace for which the retained Born weight is exactly one and for
which no ranking, no argsort and no tie-breaking is involved:

    S* = { determinants whose up-string has orbital 0 occupied }  =  supp(phi).

S* is built COMBINATORIALLY from the string tables, never by testing phi against zero: the
latter would make the subspace depend on ARPACK round-off (at L = 8 the smallest "nonzero"
amplitude is 1.3e-17).  That the two definitions agree is checked, not assumed (W3).

By construction  |S*| / D = 1/2 + 1/L  EXACTLY  and  w_{S*} = 1  EXACTLY.  It then measures,
on a dense frequency grid and by exact diagonalisation of both H and H_{S*} = P H P:

    rel-L1(eta) = ||A - A_{S*}||_1 / ||phi||^2 ,
    Lambda_hat_{S*}(eta)/eta  and the Theorem-B bound  (1+sqrt(w_S))(sqrt(1-w_S) + Lam_hat/eta),
    against the trivial bound  1 + w_S = 2.

WHY IT EXISTS (this replaces a retracted measurement)
-----------------------------------------------------
The "top-phi 70/90/99%" subspaces used earlier to argue that shots cannot buy the second term of
the bound were padded with determinants of amplitude EXACTLY ZERO, selected by argsort ties
(55 of 255 at L = 6, FR = 0.85).  Ten different tie-breaks gave rel-L1 between 0.272 and 0.409,
so those numbers were never citable to three digits and have been withdrawn.  S* = supp(phi) is
the canonical repair: the same w_S = 1 situation, reached without any ranking at all, hence
reproducible bit for bit.

WHAT IT SHOWS, AND WHAT IT MUST NOT BE READ AS
-----------------------------------------------
The weight summand is IDENTICALLY ZERO here and the error is still O(1).  That is the content:
finite-shot sampling, which only ever buys weight, cannot buy this error.

It is NOT a power law.  There is a hard ceiling rel-L1 <= 1 + w_S = 2 (self-check W5), and the
measured local slopes d log(rel-L1) / d log(eta) vary by more than a factor of two across the
grid, so "the error scales as 1/eta" is refutable and is not claimed anywhere in the output.
The local slopes are reported precisely so that the claim cannot be made by accident.

METHOD / DETERMINISM
--------------------
Engine: src/akw_lanczos.py (build_H_explicit, cdag_map) -- the same PBC ring that produced
data/akw_lanczos_L12.json.  The ground state is an ARPACK eigsh started from a FIXED-seed
vector (SEED below): without it eigsh seeds itself at random and the last digits of E0, and with
them every rel-L1, move between runs.  Everything downstream is dense LAPACK, hence exact.
No container paths: the output path is resolved from __file__.

USAGE
-----
    python src/support_witness.py                 # L = 6, 8 -> data/support_witness.json
    python src/support_witness.py --L 6           # L = 6 only (seconds)
    python src/support_witness.py --check         # compute and verify, write nothing

SELF-VERIFICATION (the script exits non-zero and writes nothing if any of these fails)
    W1  |S*| == D * (1/2 + 1/L) as an EXACT integer identity
    W2  |w_{S*} - 1| <= 1e-14                            (support carries all the weight)
    W3  phi vanishes to machine precision OUTSIDE the combinatorially defined S*, i.e. the
        operator really does support phi where the counting argument says it does.  The
        smallest amplitude INSIDE S* is reported as a diagnostic and is NOT a criterion:
        at L = 8 it is 1.3e-17, which is ARPACK round-off, and a witness must not depend on it
    W4  the rel-L1 table reproduces the established values at eta = 0.18 t:
        L = 6 -> 0.431 and L = 8 -> 0.327, atol 0.005
    W5  rel-L1 <= 1 + w_S = 2 at every eta                (the ceiling that forbids a 1/eta law)
    W6  Lambda computed here agrees with src/leakage_certificate.py certify(), an INDEPENDENT
        implementation of the same quantity, to 1e-8 relative
    W7  the windowed sum rule: the deficit of int A dw against ||phi||^2 is at most the
        ANALYTIC truncated-Lorentzian tail (2/pi) arctan(eta/PAD), so the reported rel-L1
        is a property of the subspace and not of the frequency window

RUNTIME (measured from the repository root): L = 6 about 4 s; L = 8 about 184 s
(one dense eigh of 3920 and one of 2450, plus the certify() cross-check); total
188 s, under 1 GB.  No QPU.
"""
import argparse
import datetime
import json
import os
import platform
import sys
import time

import numpy as np
import scipy

# --- NUMPY_TRAPEZOID_BRIDGE ------------------------------------------------
# numpy 2.0 ADDED np.trapezoid and REMOVED np.trapz.  Files in this repository use
# both names, so without this bridge no single numpy version runs the whole deposit.
if not hasattr(np, "trapezoid"):
    np.trapezoid = np.trapz          # numpy < 2.0
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid          # numpy >= 2.0
# ---------------------------------------------------------------------------

from scipy.sparse.linalg import eigsh

SRC_DIR = os.path.dirname(os.path.abspath(__file__))
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)
REPO = os.path.dirname(SRC_DIR)
DEFAULT_OUT = os.path.join(REPO, "data", "support_witness.json")

import akw_lanczos as AK                        # noqa: E402
from leakage_certificate import certify         # noqa: E402  (independent Lambda, self-check W6)

U_HUB = 8.0
T_HOP = 1.0
SEED = 20260918            # fixed ARPACK start vector -> bit-reproducible E0
ETAS = (0.05, 0.10, 0.15, 0.18, 0.25, 0.50)
NGRID = 60001              # dense frequency grid; the window is padded by PAD on both sides
PAD = 12.0
AMP_TOL = 1e-12            # W3: smallest admissible nonzero amplitude
REF_RELL1_AT_018 = {6: 0.431, 8: 0.327}         # W4, from the established-truth record
REF_ATOL = 0.005

_T0 = time.time()


def log(*a):
    print("[%7.1fs]" % (time.time() - _T0), *a, flush=True)


def fail(msg):
    sys.stderr.write("FAIL: " + msg + "\n")
    sys.exit(1)


def _arpack_v0(n):
    v = np.random.default_rng(SEED).standard_normal(n)
    return v / np.linalg.norm(v)


def build(L):
    """Ground state of (L/2, L/2), seed phi = c^dag_{0,up}|Psi_0> and H in (L/2+1, L/2)."""
    nup = nd = L // 2
    H0, Du0, Dd0 = AK.build_H_explicit(L, U_HUB, nup, nd, T_HOP)
    d0 = Du0 * Dd0
    E0, V0 = eigsh(H0, k=1, which="SA", v0=_arpack_v0(d0), ncv=min(d0, 40),
                   tol=1e-11, maxiter=100000)
    E0 = float(E0[0])
    Psi = V0[:, 0].reshape(Du0, Dd0)
    cdU = AK.cdag_map(L, nup)
    phi = np.asarray((cdU[0] @ Psi).reshape(-1), dtype=float)
    H1, Du1, Dd1 = AK.build_H_explicit(L, U_HUB, nup + 1, nd, T_HOP)

    # S* is defined COMBINATORIALLY, not by testing phi against zero.  phi = c^dag_{0,up}|Psi_0>
    # is supported, by the structure of the operator, exactly on the determinants whose up-string
    # has orbital 0 occupied -- and the sector index of a basis state is i = iu * Dd + id.  Taking
    # S* = {i : phi_i != 0} instead would make the subspace depend on ARPACK round-off: at L = 8
    # the smallest "nonzero" amplitude that way is 1.3e-17, i.e. noise, not physics.
    Sup, _ = AK.strings(L, nup + 1)
    occ0 = np.array([i for i, m in enumerate(Sup) if (m >> 0) & 1], dtype=np.int64)
    Sstar = (occ0[:, None] * Dd1 + np.arange(Dd1)[None, :]).ravel()
    Sstar.sort()
    return E0, phi, H1.tocsr(), Du1 * Dd1, Sstar


def spectral(poles, weights, grid, eta):
    A = np.zeros_like(grid)
    for a, b in zip(poles, weights):
        A += b * (eta / np.pi) / ((grid - a) ** 2 + eta ** 2)
    return A


def run_L(L):
    E0, phi, H, D, Sstar = build(L)
    n2 = float(phi @ phi)
    a = np.abs(phi)
    nS = int(Sstar.size)
    log("L=%d  D=%d  |S*|=%d  |S*|/D=%.12f   (1/2+1/L = %.12f)"
        % (L, D, nS, nS / D, 0.5 + 1.0 / L))

    # ---- W1: the support size is an exact integer identity, not a rounded fraction
    expected = D * (L + 2) // (2 * L)
    if D * (L + 2) % (2 * L) or nS != expected:
        fail("(W1) L=%d: |S*| = %d, expected D*(1/2+1/L) = %d" % (L, nS, expected))
    # ---- W3: the combinatorial S* really does contain all the amplitude.  phi must vanish
    # OUTSIDE it to machine precision.  The smallest amplitude INSIDE is only a diagnostic: at
    # L = 8 it is 1.3e-17, i.e. ARPACK round-off on a component that is structurally zero, and
    # making a pass criterion out of it would make the witness depend on that noise.
    amin = float(a[Sstar].min())
    off = a[np.setdiff1d(np.arange(D), Sstar)] if nS < D else np.zeros(0)
    amax_off = float(off.max()) if off.size else 0.0
    if amax_off > AMP_TOL * max(1.0, float(a.max())):
        fail("(W3) L=%d: max|phi| OUTSIDE the combinatorial S* is %.3e -- phi is not supported "
             "where the operator says it is" % (L, amax_off))
    log("L=%d  |phi| outside S* <= %.3e (machine zero); smallest |phi| inside S* = %.3e "
        "(diagnostic only)" % (L, amax_off, amin))

    # ---- w_S is 1 by construction (W2)
    wS = float(phi[Sstar] @ phi[Sstar]) / n2
    if abs(wS - 1.0) > 1e-14:
        fail("(W2) L=%d: w_S* = %.16f is not 1" % (L, wS))
    log("L=%d  w_S* = %.15f  (exact by construction, no ties)" % (L, wS))

    # ---- exact Lehmann pairs of H and Rayleigh-Ritz pairs of H_S*
    Hd = H.toarray()
    Em, Um = np.linalg.eigh(Hd)
    c_ex = Um.T @ phi
    w_ex = c_ex ** 2
    del Um
    HS = Hd[np.ix_(Sstar, Sstar)]
    del Hd
    th, Us = np.linalg.eigh(HS)
    cs = Us.T @ phi[Sstar]
    w_s = cs ** 2
    log("L=%d  dense spectra ready (D=%d, |S*|=%d)" % (L, D, nS))

    # ---- leakage Gram: g_m = Q H |m~>
    M = np.zeros((D, nS))
    M[Sstar, :] = Us
    G = H @ M
    G[Sstar, :] = 0.0
    Gram = G.T @ G
    del M, G, Us, HS

    lo = min(Em.min(), th.min()) - E0 - PAD
    hi = max(Em.max(), th.max()) - E0 + PAD
    grid = np.linspace(lo, hi, NGRID)

    rows = []
    dE = th[None, :] - th[:, None]
    for eta in ETAS:
        A = spectral(Em - E0, w_ex, grid, eta)
        AS = spectral(th - E0, w_s, grid, eta)
        l1 = float(np.trapezoid(np.abs(A - AS), grid))
        win = float(np.trapezoid(A, grid))
        rel_full = l1 / n2
        rel_win = l1 / win
        Kmat = 2j * eta / (dE.T + 2j * eta)
        val = np.einsum("m,n,nm,nm->", cs, cs, Gram, Kmat)
        lam = float(np.sqrt(abs(val.real)))
        lam_hat_eta = lam / np.sqrt(n2) / eta
        B = (1.0 + np.sqrt(wS)) * (np.sqrt(max(0.0, 1.0 - wS)) + lam_hat_eta)
        triv = 1.0 + wS
        # ---- W7: the sum-rule deficit must be exactly the truncated Lorentzian tail.
        # Every pole sits at least PAD away from both window edges by construction, so the
        # mass outside the window is at most (2/pi) arctan(eta/PAD).  Checking against that
        # ANALYTIC bound -- rather than against a round number -- is what makes this a test:
        # a deficit larger than the tail would mean missing spectral weight, not a window.
        window_loss = abs(win - n2) / n2
        tail_bound = (2.0 / np.pi) * np.arctan(eta / PAD)
        if window_loss > tail_bound:
            fail("(W7) L=%d eta=%.3f: sum-rule deficit %.4f%% exceeds the analytic Lorentzian "
                 "tail bound %.4f%% -- weight is missing for a reason other than the window"
                 % (L, eta, 100 * window_loss, 100 * tail_bound))
        # ---- W5: the hard ceiling that forbids a power law
        if rel_full > triv + 1e-9:
            fail("(W5) L=%d eta=%.3f: rel-L1 = %.6f exceeds the ceiling 1+w_S = %.6f"
                 % (L, eta, rel_full, triv))
        rows.append(dict(eta=eta, relL1=rel_full, relL1_windowed=rel_win,
                         Lambda_hat_over_eta=lam_hat_eta, bound_theoremB=B,
                         bound_trivial=triv, certifies=bool(B < triv),
                         window_loss_rel=window_loss,
                         window_loss_analytic_bound=float(tail_bound)))
        log("L=%d  eta=%.3f  rel-L1=%.6f  Lam^/eta=%.4f  B=%.4f  trivial=%.4f  -> %s"
            % (L, eta, rel_full, lam_hat_eta, B, triv, "CERT" if B < triv else "vacuous"))

    # ---- W4: reproduce the established value at eta = 0.18
    r018 = [r for r in rows if abs(r["eta"] - 0.18) < 1e-12]
    if r018 and L in REF_RELL1_AT_018:
        got, ref_v = r018[0]["relL1"], REF_RELL1_AT_018[L]
        if abs(got - ref_v) > REF_ATOL:
            fail("(W4) L=%d: rel-L1(eta=0.18) = %.6f vs established %.3f +- %.3f"
                 % (L, got, ref_v, REF_ATOL))
        log("(W4) PASS  L=%d rel-L1(eta=0.18) = %.6f vs established %.3f" % (L, got, ref_v))

    # ---- W6: Lambda against the INDEPENDENT implementation in leakage_certificate.py
    eta_x = 0.18
    ref = certify(H, phi, Sstar, eta_x, E0=E0, dense_max=D + 1,
                  want_K_S=False, kry_check=False)
    mine = [r for r in rows if abs(r["eta"] - eta_x) < 1e-12][0]
    lam_mine = mine["Lambda_hat_over_eta"] * np.sqrt(n2) * eta_x
    rel = abs(lam_mine - ref["Lambda"]) / max(ref["Lambda"], 1e-300)
    if rel > 1e-8:
        fail("(W6) L=%d: Lambda here = %.12f vs leakage_certificate.certify = %.12f (rel %.2e)"
             % (L, lam_mine, ref["Lambda"], rel))
    log("(W6) PASS  L=%d Lambda = %.12f matches leakage_certificate.certify (rel %.2e)"
        % (L, lam_mine, rel))

    # ---- local log-slopes: reported so that "1/eta" cannot be claimed by accident
    e_arr = np.array([r["eta"] for r in rows])
    r_arr = np.array([r["relL1"] for r in rows])
    slopes = (np.diff(np.log(r_arr)) / np.diff(np.log(e_arr))).tolist()
    # the slopes are negative: the meaningful spread is the ratio of MAGNITUDES
    mags = [abs(x) for x in slopes]
    spread = float(max(mags) / min(mags)) if min(mags) > 0 else float("inf")
    log("L=%d  local d log(rel-L1)/d log(eta): %s   (spread x%.2f; a pure 1/eta law would be "
        "-1.000 everywhere)" % (L, " ".join("%.3f" % x for x in slopes), spread))

    return dict(L=L, D=int(D), nS=nS, frac=nS / D, frac_exact="1/2 + 1/L",
                E0=E0, norm_phi2=n2, w_S=wS,
                amplitude_gap=dict(min_abs_phi_on_S=amin, max_abs_phi_off_S=amax_off,
                                   note=("min_abs_phi_on_S is a DIAGNOSTIC: at L=8 it is 1.3e-17, "
                                         "ARPACK round-off on a structurally zero component. S* "
                                         "is combinatorial, so this number cannot move it.")),
                table=rows,
                local_log_slopes_relL1_vs_eta=slopes,
                local_slope_spread_factor=spread,
                slope_note=("reported to make the point NEGATIVELY: the slopes are not constant, "
                            "so the eta dependence is NOT a power law and must never be quoted "
                            "as one.  rel-L1 <= 1 + w_S = 2 is a hard ceiling."),
                cross_check_leakage_certificate=dict(
                    eta=eta_x, Lambda_here=lam_mine, Lambda_certify=ref["Lambda"],
                    rel_difference=rel))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--L", type=int, nargs="+", default=[6, 8])
    ap.add_argument("--out", default=DEFAULT_OUT,
                    help="output JSON (default: <repo>/data/support_witness.json); a relative "
                         "path is taken relative to the CURRENT DIRECTORY, never to the repository")
    ap.add_argument("--check", action="store_true", help="compute and verify, write nothing")
    args = ap.parse_args()

    log("support_witness: L=%s  U/t=%.1f  ARPACK seed=%d  etas=%s"
        % (args.L, U_HUB, SEED, list(ETAS)))
    log("numpy %s  scipy %s  python %s" % (np.__version__, scipy.__version__,
                                           platform.python_version()))
    results = [run_L(L) for L in args.L]

    print("")
    print("  L   |S*|/D      " + "".join("%9.2f" % e for e in ETAS))
    for r in results:
        print("%3d  %8.6f  " % (r["L"], r["frac"])
              + "".join("%9.4f" % row["relL1"] for row in r["table"]))
    print("  (entries are rel-L1 = ||A - A_S*||_1 / ||phi||^2 at w_S = 1 exactly)")

    if args.check:
        print("--check: nothing written.")
        return

    out = dict(
        provenance=dict(
            generated_by="src/support_witness.py",
            generated_utc=datetime.datetime.now(datetime.timezone.utc)
                            .strftime("%Y-%m-%dT%H:%M:%SZ"),
            runtime_seconds=round(time.time() - _T0, 1),
            python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
            platform=platform.platform(),
            engine="src/akw_lanczos.py build_H_explicit / cdag_map (1D Hubbard ring, PBC)",
            model="1D Hubbard, PBC, half filling, U/t=%.1f, t=%.1f" % (U_HUB, T_HOP),
            seed_operator="phi = c^dag_{0,up}|Psi_0> in the (N+1, Sz=+1) sector",
            subspace="S* = supp(phi): no ranking, no argsort, no tie-break; w_S = 1 exactly",
            arpack_seed=SEED, deterministic=True,
            parameters=dict(U_over_t=U_HUB, t=T_HOP, sizes=list(args.L), etas=list(ETAS),
                            n_grid=NGRID, window_pad=PAD, arpack_seed=SEED,
                            amplitude_tolerance=AMP_TOL),
            control_values_reproduced={
                str(r["L"]): dict(
                    quantity="rel-L1 at eta = 0.18 t",
                    established=REF_RELL1_AT_018.get(r["L"]),
                    recomputed=[t["relL1"] for t in r["table"] if abs(t["eta"] - 0.18) < 1e-12][0],
                    atol=REF_ATOL)
                for r in results if r["L"] in REF_RELL1_AT_018},
            n_grid=NGRID, window_pad=PAD, etas=list(ETAS),
            checks=dict(
                W1="|S*| == D*(1/2 + 1/L) as an exact integer identity",
                W2="|w_S* - 1| <= 1e-14",
                W3=("|phi| outside the combinatorial S* is <= %g relative to max|phi|; "
                    "the smallest amplitude inside S* is a diagnostic, not a criterion" % AMP_TOL),
                W4="rel-L1(eta=0.18) reproduces %s within %g" % (REF_RELL1_AT_018, REF_ATOL),
                W5="rel-L1 <= 1 + w_S = 2 at every eta",
                W6="Lambda agrees with src/leakage_certificate.py certify() to 1e-8 relative",
                W7=("sum-rule deficit at most the analytic Lorentzian tail "
                    "(2/pi) arctan(eta/PAD), PAD = %g" % PAD)),
            checks_status="ALL PASS",
            replaces=("the withdrawn 'top-phi 70/90/99%' measurement, whose subspaces were padded "
                      "with amplitude-zero determinants selected by argsort ties (55 of 255 at "
                      "L=6, FR=0.85; ten tie-breaks gave rel-L1 between 0.272 and 0.409)"),
            reading=("the weight summand of the Theorem-B bound is identically zero here and the "
                     "error is still O(1): finite-shot sampling, which only buys weight, cannot "
                     "buy this error.  Stated as a measured observation, NEVER as a power law."),
        ),
        results=results,
    )
    outp = os.path.abspath(args.out)
    os.makedirs(os.path.dirname(outp), exist_ok=True)
    with open(outp, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    log("WROTE %s" % outp)


if __name__ == "__main__":
    main()
