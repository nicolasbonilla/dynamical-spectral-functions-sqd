# -*- coding: utf-8 -*-
r"""Chemical potentials mu^+, mu^- and the charge (Mott) gap Delta of the half-filled
1D Hubbard RING (PBC, U/t = 8) by exact sector diagonalization, for L = 6, 8, 10, 12.

WHY THIS FILE EXISTS
--------------------
The manuscript prints Delta = 4.97 t and mu^± = ±2.484 t (measured from the particle-hole
symmetric Fermi level U/2), and the value 4.969 is carried into paper/figs/sqw_edges.dat.
Yet NO script in this repository computed Delta, mu^+ or mu^- : the numbers were correct but
UNPROVABLE from the deposit.  This script closes that hole.

DEFINITIONS (all energies in units of t, raw -- not shifted by the Fermi level)

    E0(N)        ground-state energy of the (n_up, n_dn) = (L/2, L/2) sector      [N   = L]
    E0(N+1)      ground-state energy of the (L/2 + 1, L/2) sector                 [N+1]
    E0(N-1)      ground-state energy of the (L/2 - 1, L/2) sector                 [N-1]

    mu^+ = E0(N+1) - E0(N)          (electron-addition chemical potential)
    mu^- = E0(N)   - E0(N-1)        (electron-removal  chemical potential)
    Delta = mu^+ - mu^-             (charge / Mott gap)

Measured from the half-filling Fermi level mu_F = U/2 (the convention used by
akw_lanczos.py:111 for the omega axis of A(k,omega)):

    mu^+_rel = mu^+ - U/2 = +Delta/2 ,   mu^-_rel = mu^- - U/2 = -Delta/2

which is where the manuscript's "+-2.484 t" comes from.  Particle-hole symmetry of the
half-filled bipartite (even-L) Hubbard ring forces mu^+ + mu^- = U exactly; this is checked.

PAPER CLAIMS THIS BACKS
-----------------------
  * "Delta = 4.97 t", "mu^± = ±2.484 t"  (Mott gap of the L=12 ring);
  * the 4.969 that the two FABRICATED rows of the old paper/figs/sqw_edges.dat carried --
    this script is what identifies that number as a ONE-PARTICLE-channel quantity, i.e. NOT a
    measurement of S(q,omega) (see src/make_sqw_edges.py);
  * the statement that the charge gap does NOT close with L: Delta = 5.358127 / 5.185842 /
    5.022235 / 4.968759 for L = 6, 8, 10, 12 -- monotonically decreasing and staying ABOVE the
    exact thermodynamic-limit value Delta_inf = 4.679517 t, which this script evaluates
    independently from the Lieb-Wu integral
        Delta_LW = U - 4t + 8t * int_0^inf dw J_1(w) / ( w (1 + exp(w U / 2t)) ).
    HONESTY NOTE: four sizes do not by themselves establish the limit, and this script performs
    NO extrapolation. What is demonstrated is monotone decrease bounded below by Delta_inf > 0,
    i.e. the charge channel does not close -- not the rate at which it converges.

ENGINE / DETERMINISM
--------------------
Hamiltonian from akw_lanczos.build_H_explicit (the same PBC ring engine that produced
data/akw_lanczos_L12.json, data/sqw_L12.json and data/spinqw_L12.json), ARPACK eigsh(k=1,'SA')
started from a FIXED seed (SEED = 12345) so every run is reproducible.  Every sector small
enough (dim <= DENSE_MAX) is cross-checked against an independent dense LAPACK diagonalization
of the same matrix.

NO container paths: everything is resolved relative to this file, so the script repopulates
data/ when run from anywhere.

USAGE
-----
    python src/charge_gap_ed.py                 # L = 6, 8, 10, 12  -> data/charge_gap.json
    python src/charge_gap_ed.py 6 8             # a subset (L = 12 must be present to verify)
    python src/charge_gap_ed.py --check         # compute and verify, write nothing

SELF-VERIFICATION (the script exits non-zero and writes nothing if any of these fails)
    (V1) |Delta(L=12) - 4.9688| <= 0.001                        <- the value the paper prints
    (V2) |mu^+ + mu^- - U| <= 1e-8  for every L                 <- particle-hole symmetry
    (V3) |E0(N, L=12) - E0 stored in data/sqw_L12.json| <= 1e-6 <- agreement with the deposit
    (V4) L <= 8: |E_ARPACK - E_dense| <= 1e-10                  <- independent dense check
    (V5) Delta(L) decreasing in L and > Delta_LW = 4.6795       <- saturation, not closure
"""
import json, os, sys, time, datetime
import numpy as np
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh as dense_eigh
from scipy.integrate import quad
from scipy.special import j1

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from akw_lanczos import build_H_explicit

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_JSON = os.path.join(REPO, 'data', 'charge_gap.json')
SQW_JSON = os.path.join(REPO, 'data', 'sqw_L12.json')

SEED = 12345
ARPACK_TOL = 1.0e-11        # ARPACK residual tolerance; Delta is then good to ~1e-10
ARPACK_NCV = 40             # Krylov basis size (fewer implicit restarts than the default 20)
DENSE_MAX = 5000            # sectors up to this dimension get an independent dense cross-check
U = 8.0
T_HOP = 1.0
DELTA_L12_REF = 4.9688      # (V1) reference value printed by the paper
DELTA_L12_TOL = 1.0e-3

t0 = time.time()


def log(*a):
    print("[%7.1fs]" % (time.time() - t0), *a, flush=True)


def fail(msg):
    sys.stderr.write("FAIL: " + msg + "\n")
    sys.exit(1)


def sector_ground_energy(L, nup, ndn):
    """Lowest eigenvalue of H in the (nup, ndn) sector of the L-site Hubbard ring."""
    H, Du, Dd = build_H_explicit(L, U, nup, ndn, T_HOP)
    dim = Du * Dd
    if dim == 1:
        return float(H.toarray()[0, 0]), dim, None
    v0 = np.random.default_rng(SEED).standard_normal(dim)
    ncv = min(dim, ARPACK_NCV)
    E = eigsh(H, k=1, which='SA', v0=v0, ncv=ncv, tol=ARPACK_TOL,
              maxiter=100000, return_eigenvectors=False)
    e_arpack = float(np.min(E))
    e_dense = None
    if dim <= DENSE_MAX:                                # (V4) independent dense cross-check
        e_dense = float(dense_eigh(H.toarray(), eigvals_only=True, subset_by_index=[0, 0])[0])
        if abs(e_arpack - e_dense) > 1e-10:
            fail("L=%d sector (%d,%d): ARPACK %.12f vs dense %.12f differ by %.2e"
                 % (L, nup, ndn, e_arpack, e_dense, abs(e_arpack - e_dense)))
    return e_arpack, dim, e_dense


def lieb_wu_gap(U_=U, t_=T_HOP):
    """Thermodynamic-limit Mott gap of the 1D Hubbard model (Lieb & Wu, PRL 20, 1445 (1968))."""
    def f(w):
        x = w * U_ / (2.0 * t_)
        return j1(w) / w * (np.exp(-x) / (1.0 + np.exp(-x))) if x > 0 else 0.0
    I, err = quad(f, 0.0, np.inf, limit=400)
    return U_ - 4.0 * t_ + 8.0 * t_ * I, err


def run(Ls):
    rows = []
    for L in Ls:
        h = L // 2
        eN, dN, dnN = sector_ground_energy(L, h, h)
        log("L=%2d  (N  ) sector (%d,%d) dim=%8d  E0=%+.12f" % (L, h, h, dN, eN))
        eP, dP, _ = sector_ground_energy(L, h + 1, h)
        log("L=%2d  (N+1) sector (%d,%d) dim=%8d  E0=%+.12f" % (L, h + 1, h, dP, eP))
        eM, dM, _ = sector_ground_energy(L, h - 1, h)
        log("L=%2d  (N-1) sector (%d,%d) dim=%8d  E0=%+.12f" % (L, h - 1, h, dM, eM))
        mu_p = eP - eN
        mu_m = eN - eM
        delta = mu_p - mu_m
        # (V2) particle-hole symmetry of the half-filled bipartite ring
        if abs(mu_p + mu_m - U) > 1e-8:
            fail("L=%d: mu+ + mu- = %.12f != U = %.1f (particle-hole symmetry broken)"
                 % (L, mu_p + mu_m, U))
        rows.append(dict(L=L, U=U, t=T_HOP,
                         dim_N=dN, dim_Np1=dP, dim_Nm1=dM,
                         E0_N=eN, E0_Np1=eP, E0_Nm1=eM,
                         E0_N_dense=dnN,
                         mu_plus=mu_p, mu_minus=mu_m,
                         mu_plus_rel=mu_p - U / 2.0, mu_minus_rel=mu_m - U / 2.0,
                         Delta=delta, L_times_Delta=L * delta))
        log("L=%2d  mu+=%.6f  mu-=%.6f  Delta=%.6f   (mu+_rel=%+.6f, mu-_rel=%+.6f)"
            % (L, mu_p, mu_m, delta, mu_p - U / 2.0, mu_m - U / 2.0))
    return rows


def verify(rows):
    by_L = {r['L']: r for r in rows}

    # (V1) the value the paper prints
    if 12 not in by_L:
        fail("L=12 was not computed, so the mandatory check (V1) cannot run")
    d12 = by_L[12]['Delta']
    if abs(d12 - DELTA_L12_REF) > DELTA_L12_TOL:
        fail("(V1) Delta(L=12) = %.6f is outside %.4f +- %.4f"
             % (d12, DELTA_L12_REF, DELTA_L12_TOL))
    log("(V1) PASS  Delta(L=12) = %.6f  vs reference %.4f +- %.4f  (residual %.2e)"
        % (d12, DELTA_L12_REF, DELTA_L12_TOL, abs(d12 - DELTA_L12_REF)))
    log("(V2) PASS  mu+ + mu- = U = %.1f for all %d sizes (particle-hole symmetry)" % (U, len(rows)))

    # (V3) agreement with the deposited S(q,omega) ground-state energy
    if os.path.exists(SQW_JSON):
        e_dep = float(json.load(open(SQW_JSON))['E0'])
        if abs(by_L[12]['E0_N'] - e_dep) > 1e-6:
            fail("(V3) E0(N,L=12) = %.12f disagrees with data/sqw_L12.json E0 = %.12f"
                 % (by_L[12]['E0_N'], e_dep))
        log("(V3) PASS  E0(N,L=12) = %.12f matches data/sqw_L12.json (%.12f)"
            % (by_L[12]['E0_N'], e_dep))
    else:
        log("(V3) SKIP  data/sqw_L12.json not present")

    log("(V4) PASS  dense cross-check agreed to 1e-10 on every sector small enough to diagonalize")

    # (V5) the charge gap must DECREASE with L and stay ABOVE the Lieb-Wu limit
    lw, lw_err = lieb_wu_gap()
    ds = [by_L[L]['Delta'] for L in sorted(by_L)]
    if len(ds) > 1 and any(b >= a for a, b in zip(ds, ds[1:])):
        fail("(V5) Delta(L) is not monotonically decreasing: %s" % ds)
    if min(ds) <= lw:
        fail("(V5) Delta(L) = %.6f fell below the Lieb-Wu limit %.6f" % (min(ds), lw))
    log("(V5) PASS  Delta(L) decreasing %s, all above Lieb-Wu Delta_inf = %.6f (quad err %.1e)"
        % (["%.4f" % x for x in ds], lw, lw_err))
    return lw, lw_err


def main():
    check_only = "--check" in sys.argv
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    Ls = [int(a) for a in args] if args else [6, 8, 10, 12]
    log("Hubbard ring (PBC), U/t=%.1f, ARPACK seed=%d, sizes %s" % (U, SEED, Ls))
    rows = run(Ls)
    lw, lw_err = verify(rows)

    print("")
    print(" L    dim(N+1)     mu^+        mu^-       Delta      mu^+_rel   mu^-_rel   L*Delta")
    for r in rows:
        print("%2d  %9d  %10.6f  %10.6f  %9.6f  %+9.6f  %+9.6f  %8.4f"
              % (r['L'], r['dim_Np1'], r['mu_plus'], r['mu_minus'], r['Delta'],
                 r['mu_plus_rel'], r['mu_minus_rel'], r['L_times_Delta']))
    print("Lieb-Wu thermodynamic limit: Delta_inf = %.6f" % lw)

    if check_only:
        print("--check: nothing written.")
        return

    out = dict(
        provenance=dict(
            generated_by="src/charge_gap_ed.py",
            generated_utc=datetime.datetime.now(datetime.timezone.utc)
                            .strftime("%Y-%m-%dT%H:%M:%SZ"),
            engine="akw_lanczos.build_H_explicit (1D Hubbard ring, PBC) + scipy.sparse.linalg.eigsh",
            model="1D Hubbard, PBC, half filling, U/t=%.1f, t=%.1f" % (U, T_HOP),
            definitions=dict(
                mu_plus="E0(n_up=L/2+1, n_dn=L/2) - E0(L/2, L/2)",
                mu_minus="E0(L/2, L/2) - E0(n_up=L/2-1, n_dn=L/2)",
                Delta="mu_plus - mu_minus",
                mu_rel="mu - U/2, the particle-hole symmetric Fermi level (akw_lanczos.py:111)"),
            arpack_seed=SEED, arpack_tol=ARPACK_TOL, arpack_ncv=ARPACK_NCV,
            deterministic=True,
            checks=dict(
                V1="|Delta(L=12) - %.4f| <= %.4f" % (DELTA_L12_REF, DELTA_L12_TOL),
                V2="|mu+ + mu- - U| <= 1e-8 (particle-hole symmetry)",
                V3="|E0(N,L=12) - data/sqw_L12.json E0| <= 1e-6",
                V4="ARPACK vs dense eigh <= 1e-10 on every sector of dim <= %d" % DENSE_MAX,
                V5="Delta(L) decreasing and above the Lieb-Wu limit"),
            checks_status="ALL PASS",
            backs=("Delta = 4.97 t and mu^± = ±2.484 t as printed in the manuscript; identifies "
                   "the 4.969 carried by the two fabricated rows of the former "
                   "paper/figs/sqw_edges.dat as a ONE-PARTICLE-channel quantity."),
        ),
        lieb_wu=dict(Delta_inf=lw, quad_abserr=lw_err,
                     formula="U - 4t + 8t * int_0^inf dw J_1(w)/(w(1+exp(wU/2t)))",
                     reference="Lieb & Wu, Phys. Rev. Lett. 20, 1445 (1968)"),
        results=rows,
    )
    json.dump(out, open(OUT_JSON, 'w'), indent=1)
    log("WROTE data/charge_gap.json")


if __name__ == '__main__':
    main()
