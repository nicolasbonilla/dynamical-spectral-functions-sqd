# -*- coding: utf-8 -*-
r"""HONEST FINITE-SHOT SAMPLING for the scaling figure (fig:scaling / Fig. 12) and Fig. 5.

WHAT THIS COMPUTES
------------------
The published scaling curve selects the subspace S with

    order = np.argsort(wc)[::-1]                # scaling_lanczos.py:66

i.e. by the EXACT Born distribution wc of the time-evolved state, which is the
infinite-shot limit of the declared protocol.  The manuscript calls those subspaces
"sampled", and the abstract says "we reconstruct, from sampled subspaces, ...".
This script replaces that selection by REAL finite-shot sampling: the subspace is the
union of the DISTINCT bitstrings actually observed in T shots split evenly over the
K+1 time slices -- exactly the template already used in the repository by
src/spectral_validation_max.py:50-62 and src/headtohead_ms.py:46-55.

Everything else is identical to src/scaling_lanczos.py:compute():
  1D Hubbard, PBC, U/t = 8, half filling, seed = c^dag_{0,up}|GS> in the (N+1, Sz+1)
  sector, eta = 0.15 t, K = 18 time slices, dt = 0.5, Haydock depth nl = 260,
  600-point frequency grid, acceptance threshold rel-L1 < 0.05.  The Hamiltonian,
  sector strings, c^dag map and the Haydock continued fraction are imported from
  src/akw_lanczos.py -- nothing is re-implemented.

WHICH CLAIM OF THE PAPER THIS SUPPORTS
--------------------------------------
  * the abstract's "from sampled subspaces" and the caption of Fig. 5
    ("not the exact target of Fig. 4, but the *sampled* output", results.tex:85);
  * the resource statement of Fig. 12 / results.tex:299 ("Fraction down, absolute cost
    up -- the honest two-sided picture"), which is here completed with its THIRD side:
    the shot budget T.
It produces, for L = 4, 6, 8:
  (1) sector fraction |S|/dim and rel-L1 versus shot budget T (log sweep, >= 5 values
      per L), mean +- s.d. over 8 independent RNG seeds;
  (2) the minimum T that reaches rel-L1 < 0.05, with its seed-to-seed band;
  (3) the SUBSPACE-MATCHED comparison: at exactly the same |S|, what rel-L1 does the
      infinite-shot ordering give and what does finite-shot sampling give;
  (4) the exact expected number of distinct configurations
      E[#distinct](T) = sum_i [1 - prod_k (1-p_{k,i})^{T/(K+1)}],
      evaluated in closed form (no Monte Carlo), so that larger L can be extrapolated
      without spending CPU, plus its inverse T(|S|).

SELF-CHECK (the run ABORTS if any of these fails)
-------------------------------------------------
  C1  infinite-shot ("oracle") fraction at eta = 0.15 must reproduce the published
      0.9167 / 0.82 / 0.5561 for L = 4 / 6 / 8   (tolerance below, ORACLE_TOL);
  C2  the Haydock reference A(omega) must agree with a dense Lehmann sum over the full
      sector to rel-L1 < 1e-6 (done for L <= 6, where dense diagonalisation is cheap);
  C3  the exact spectral weight integral must reproduce ||c^dag_{0,up}|GS>||^2 to
      better than 6% over the finite frequency window.  The residual IS the Lorentzian
      tail truncated by the window and is reported explicitly as window_loss_rel -- it
      is the same effect that makes the "0.500 to machine precision" sum-rule claim of
      the Fig. 4 caption untenable.

HOW TO RUN
----------
    python src/sampled_honest.py                   # L = 4, 6, 8 -> data/sampled_honest.json
    python src/sampled_honest.py --L 4 6 --out other.json
    python src/sampled_honest.py --nseeds 8 --smax-cap 300000

Run it from anywhere: the engine (src/akw_lanczos.py) and the output path are both
resolved from __file__.  No container paths.

RUNTIME (re-measured 2026-09-18 by running this file from the repository root)
-----------------------------------------------------------------------------
    L = 4 : 1 s     L = 6 : 9 s     L = 8 : 57 s      total 67 s, < 0.6 GB RAM.
    Peak memory is the shot streams: n_seeds x 19 x smax int32 (83 MB at L = 8).
    No QPU, no cluster, no GPU.

All RNG seeds are fixed and explicit (SEED_BASE below); the run is bit-for-bit
reproducible on a given numpy version.
"""
import argparse, json, os, platform, sys, time
import numpy as np
import scipy
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh, expm_multiply

# --- NUMPY_TRAPEZOID_BRIDGE ------------------------------------------------
# numpy 2.0 ADDED np.trapezoid and REMOVED np.trapz.  Files in this repository use
# both names, so without this bridge no single numpy version runs the whole deposit.
if not hasattr(np, "trapezoid"):
    np.trapezoid = np.trapz          # numpy < 2.0
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid          # numpy >= 2.0
# ---------------------------------------------------------------------------

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import akw_lanczos as AK

# Repository layout: this file lives in <repo>/src/, the deposit in <repo>/data/.
# NO container paths: the default output is resolved from __file__, so the script
# repopulates data/ when launched from anywhere.
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_OUT = os.path.join(REPO, "data", "sampled_honest.json")

# ----------------------------------------------------------------------------- parameters
U_HUB   = 8.0          # on-site repulsion, in units of t
ETA     = 0.15         # Lorentzian broadening, in units of t   (scaling_lanczos.py default)
K_SLICE = 18           # number of time slices AFTER t=0  -> K+1 = 19 distributions
DT      = 0.5          # time step
THR     = 0.05         # acceptance threshold on rel-L1
NL      = 260          # Haydock depth
NGRID   = 600          # frequency grid points
NSEEDS  = 8            # independent RNG seeds
SEED_BASE = 20260918   # fixed; change only if you WANT a different realisation

# published infinite-shot fractions this script must reproduce (control C1)
ORACLE_REF = {4: 0.9167, 6: 0.82, 8: 0.5561}
ORACLE_TOL = 0.02      # absolute, on the sector fraction

# subspace sizes at which the established-truth document reports the matched comparison
MATCHED_REF = {6: 247, 8: 2331}

t0 = time.time()
def log(*a):
    print(f"[{time.time()-t0:7.1f}s]", *a, flush=True)

_trapz = getattr(np, "trapezoid", None) or np.trapz


def _arpack_start(n):
    """Fixed, non-degenerate ARPACK start vector -> bit-reproducible eigenvalues."""
    v = np.random.default_rng(SEED_BASE).standard_normal(n)
    return v / np.linalg.norm(v)


# ----------------------------------------------------------------------------- physics
def build_problem(L):
    """Ground state, (N+1,Sz+1) sector Hamiltonian, seed vector, frequency grid, exact A."""
    nup = nd = L // 2
    H0, Du0, Dd0 = AK.build_H_explicit(L, U_HUB, nup, nd)
    # deterministic ARPACK start vector: without it eigsh seeds itself randomly and the
    # last ~3 digits of E0 (hence of every rel-L1) move from run to run.
    v0_0 = _arpack_start(H0.shape[0])
    w, v = eigsh(H0, k=1, which='SA', v0=v0_0)
    E0 = float(w[0]); psi_mat = v[:, 0].reshape(Du0, Dd0)

    cdU = AK.cdag_map(L, nup)                       # c^dag_j : nup -> nup+1
    H1, Du1, Dd1 = AK.build_H_explicit(L, U_HUB, nup + 1, nd)
    nS = H1.shape[0]
    seed = (cdU[0] @ psi_mat).reshape(-1).astype(complex)
    sumrule = float(np.vdot(seed, seed).real)       # = 1 - <n_{0,up}> = 0.5 at half filling

    v0_1 = _arpack_start(nS)
    emin = float(eigsh(H1, k=1, which='SA', v0=v0_1, return_eigenvectors=False)[0])
    emax = float(eigsh(H1, k=1, which='LA', v0=v0_1, return_eigenvectors=False)[0])
    grid = np.linspace(emin - E0 - 1.0, emax - E0 + 1.0, NGRID)

    def full_mv(x):
        return H1.dot(x)

    A_ex = spec(full_mv, seed, E0, grid, min(NL, nS))
    nrm = float(_trapz(np.abs(A_ex), grid))
    integral = float(_trapz(A_ex, grid))
    log(f"L={L}: E0={E0:.6f}  dim(N+1,Sz+1)={nS}  ||seed||^2={sumrule:.6f}  "
        f"int A_ex={integral:.6f}  window=[{grid[0]:.2f},{grid[-1]:.2f}]")
    return dict(L=L, nup=nup, nd=nd, E0=E0, H1=H1, nS=nS, seed=seed, grid=grid,
                A_ex=A_ex, nrm=nrm, sumrule=sumrule, integral=integral)


def spec(mv, v0, E0, grid, nl):
    """A(omega) = -Im <v0|(omega + i eta - (H - E0))^{-1}|v0> / pi  via Haydock."""
    z = grid + 1j * ETA
    G = AK.haydock(lambda x: mv(x) - E0 * x, v0, nl, z)
    return -G.imag / np.pi


def rel_l1(P, Sset):
    """rel-L1 of the Rayleigh-Ritz spectral function on the coordinate subspace Sset.

    Restricted matvec: nothing is ever materialised as a submatrix (same trick as
    scaling_lanczos.py:relL1_at)."""
    S = np.sort(np.asarray(Sset, dtype=np.int64))
    H1 = P['H1']; nS = P['nS']
    def rmv(xs):
        xf = np.zeros(nS, dtype=complex); xf[S] = xs
        return H1.dot(xf)[S]
    A_s = spec(rmv, P['seed'][S], P['E0'], P['grid'], min(NL, len(S)))
    return float(_trapz(np.abs(A_s - P['A_ex']), P['grid']) / P['nrm'])


def born_slices(P):
    """The K+1 Born distributions p_k = |e^{-i k dt H} seed|^2 / norm that the device
    would sample, and their accumulated weight wc (the infinite-shot ranking)."""
    nS = P['nS']
    vt = P['seed'] / np.sqrt(np.vdot(P['seed'], P['seed']).real)
    pr = np.abs(vt) ** 2
    slices = [pr / pr.sum()]
    wc = pr.copy()
    for k in range(1, K_SLICE + 1):
        vt = expm_multiply(-1j * DT * P['H1'], vt)
        pr = np.abs(vt) ** 2
        slices.append(pr / pr.sum())
        wc += pr
    return np.array(slices), wc


# ----------------------------------------------------------------------------- oracle
def oracle_fraction(P, order):
    """Reproduce scaling_lanczos.py's coarse scan + bisection on the infinite-shot order."""
    nS = P['nS']
    fracs = [0.02, 0.05, 0.08, 0.12, 0.16, 0.20, 0.25, 0.30, 0.40, 0.50,
             0.60, 0.70, 0.78, 0.85, 0.92, 1.0]
    ks = sorted(set(max(2, int(f * nS)) for f in fracs))
    frac, klo, kstar = 1.0, 2, nS
    for k in ks:
        r = rel_l1(P, order[:k])
        log(f"  oracle L={P['L']}: frac={k/nS:.4f} (k={k}) relL1={r:.4f}")
        if r < THR:
            khi = k
            while khi - klo > max(2, nS // 200):
                km = (klo + khi) // 2
                if rel_l1(P, order[:km]) < THR: khi = km
                else: klo = km
            frac = khi / nS; kstar = khi
            break
        klo = k
    return float(frac), int(kstar)


# ----------------------------------------------------------------------------- shots
def shot_streams(P, slices, rng, smax):
    """One ordered stream of smax draws per time slice (inverse-CDF sampling).

    Prefixes of these streams ARE the experiment at smaller budgets, so |S|(T) is a
    monotone nested family -- which makes the bisection on T well defined and costs
    no extra random numbers."""
    streams = []
    for pr in slices:
        cdf = np.cumsum(pr); cdf[-1] = 1.0
        u = rng.random(smax)
        streams.append(np.searchsorted(cdf, u).astype(np.int32))
    return streams


def observed_union(streams, s, dominant):
    """Distinct bitstrings seen in the first s shots of each slice, plus the dominant
    determinant (the house template seeds the set with it: spectral_validation_max.py:54)."""
    if s <= 0:
        return np.array([dominant], dtype=np.int64)
    seen = np.unique(np.concatenate([st[:s] for st in streams]))
    if dominant not in seen:
        seen = np.append(seen, dominant)
    return np.sort(seen.astype(np.int64))


def observed_topm(streams, s, dominant, m):
    """Size-matched finite-shot subspace: the m most FREQUENTLY observed bitstrings in
    the first s shots of each slice (the estimator used by headtohead_ms.py:53).
    Uses only observed data -- no access to the exact Born weights."""
    cat = np.concatenate([st[:s] for st in streams])
    vals, cnt = np.unique(cat, return_counts=True)
    cnt = cnt.astype(np.float64)
    if dominant in vals:
        cnt[np.searchsorted(vals, dominant)] += 1e9      # the seeded determinant ranks first
    else:
        vals = np.append(vals, dominant); cnt = np.append(cnt, 1e9)
        o = np.argsort(vals); vals, cnt = vals[o], cnt[o]
    keep = vals[np.argsort(-cnt, kind='stable')[:m]]
    return np.sort(keep.astype(np.int64))


def e_distinct(slices, s):
    """EXACT E[#distinct] after s shots per slice (no Monte Carlo):
    sum_i [1 - prod_k (1 - p_{k,i})^s]."""
    logmiss = np.zeros(slices.shape[1])
    for pr in slices:
        logmiss += s * np.log1p(-np.clip(pr, 0.0, 1.0 - 1e-16))
    return float(np.sum(1.0 - np.exp(logmiss)))


def invert_e_distinct(slices, target, smax=10 ** 12):
    """Smallest s per slice with E[#distinct] >= target (bisection on the closed form)."""
    if e_distinct(slices, smax) < target:
        return None
    lo, hi = 1, 1
    while e_distinct(slices, hi) < target:
        hi *= 2
        if hi > smax: return None
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if e_distinct(slices, mid) >= target: hi = mid
        else: lo = mid
    return int(hi)


# ----------------------------------------------------------------------------- driver
def run_L(L, nseeds, smax_cap, tsweep_points):
    P = build_problem(L)
    nS = P['nS']
    slices, wc = born_slices(P)
    order = np.argsort(wc)[::-1]
    dominant = int(order[0])

    # ---- control C1: the infinite-shot ("oracle") fraction ----
    frac_or, k_or = oracle_fraction(P, order)
    ref = ORACLE_REF[L]
    ok = abs(frac_or - ref) <= ORACLE_TOL
    log(f"CONTROL C1  L={L}: oracle frac={frac_or:.4f} (|S|={k_or}) vs published {ref} "
        f"-> {'OK' if ok else 'MISMATCH'}")
    if not ok:
        raise SystemExit(f"ABORT: control C1 failed at L={L}: {frac_or:.4f} vs {ref} "
                         f"(tol {ORACLE_TOL}). Do not trust anything downstream.")

    # ---- size the shot streams from the closed form, not by trial and error ----
    # headroom: finite-shot subspaces need MORE configurations than the infinite-shot
    # ranking to reach the same accuracy, and the coupon-collector tail is heavy.
    target_sz = min(0.97 * nS, k_or + 0.12 * nS)
    s_need = invert_e_distinct(slices, target_sz)
    if s_need is None:                      # target unreachable (configs of zero Born weight)
        s_need = 2 * (invert_e_distinct(slices, k_or) or 1024)
    smax = int(min(smax_cap, max(128, 6 * s_need)))
    log(f"L={L}: E[#distinct] inversion -> s={s_need} per slice for |S|={target_sz:.0f}; "
        f"using smax={smax} (T_max={smax*(K_SLICE+1)})")

    # ---- (4) exact expected-distinct curve ----
    s_grid = np.unique(np.round(np.logspace(0, np.log10(smax), 24)).astype(int))
    edist = [{'shots_per_slice': int(s), 'T': int(s * (K_SLICE + 1)),
              'E_distinct': e_distinct(slices, int(s)),
              'E_fraction': e_distinct(slices, int(s)) / nS} for s in s_grid]
    T_match_published = invert_e_distinct(slices, k_or)
    log(f"L={L}: exact E[#distinct] -> T needed to MATCH the published |S|={k_or} is "
        f"{None if T_match_published is None else T_match_published*(K_SLICE+1):.3g}")

    # ---- (1) T sweep, mean +- sd over seeds ----
    s_sweep = np.unique(np.round(np.logspace(0, np.log10(smax), tsweep_points)).astype(int))
    rngs = [np.random.default_rng(SEED_BASE + sd) for sd in range(nseeds)]
    streams_all = [shot_streams(P, slices, rngs[sd], smax) for sd in range(nseeds)]
    log(f"L={L}: shot streams drawn ({nseeds} seeds x {K_SLICE+1} slices x {smax} shots)")

    sweep = []
    for s in s_sweep:
        sizes, errs = [], []
        for sd in range(nseeds):
            Sset = observed_union(streams_all[sd], int(s), dominant)
            sizes.append(len(Sset)); errs.append(rel_l1(P, Sset))
        sweep.append({'shots_per_slice': int(s), 'T': int(s * (K_SLICE + 1)),
                      'S_mean': float(np.mean(sizes)), 'S_std': float(np.std(sizes)),
                      'frac_mean': float(np.mean(sizes)) / nS,
                      'frac_std': float(np.std(sizes)) / nS,
                      'relL1_mean': float(np.mean(errs)), 'relL1_std': float(np.std(errs)),
                      'relL1_min': float(np.min(errs)), 'relL1_max': float(np.max(errs)),
                      'n_pass': int(np.sum(np.array(errs) < THR))})
        log(f"L={L} sweep: T={s*(K_SLICE+1):>9d}  frac={np.mean(sizes)/nS:.4f}"
            f"+-{np.std(sizes)/nS:.4f}  relL1={np.mean(errs):.4f}+-{np.std(errs):.4f}"
            f"  pass {sweep[-1]['n_pass']}/{nseeds}")

    # ---- (2) minimum T reaching rel-L1 < THR, per seed ----
    per_seed_T, per_seed_frac, censored = [], [], 0
    for sd in range(nseeds):
        st = streams_all[sd]
        lo, hi = 0, None
        for s in s_sweep:                                   # coarse bracket from the sweep
            if rel_l1(P, observed_union(st, int(s), dominant)) < THR:
                hi = int(s); break
            lo = int(s)
        if hi is None:
            censored += 1
            per_seed_T.append(None); per_seed_frac.append(None)
            continue
        while hi - lo > max(1, hi // 100):                  # 1% resolution in T
            mid = (lo + hi) // 2
            if rel_l1(P, observed_union(st, mid, dominant)) < THR: hi = mid
            else: lo = mid
        Sfin = observed_union(st, hi, dominant)
        per_seed_T.append(int(hi * (K_SLICE + 1)))
        per_seed_frac.append(len(Sfin) / nS)
        log(f"L={L} minT seed {sd}: T*={hi*(K_SLICE+1)}  |S|={len(Sfin)} "
            f"frac={len(Sfin)/nS:.4f}")
    good_T = [x for x in per_seed_T if x is not None]
    good_f = [x for x in per_seed_frac if x is not None]

    # ---- (3) subspace-matched comparison, oracle vs finite shots ----
    matched = []
    targets = sorted(set([k_or] + ([MATCHED_REF[L]] if L in MATCHED_REF else [])))
    for m in targets:
        if m > nS: continue
        or_err = rel_l1(P, order[:m])
        # smallest budget whose union can supply m configurations (per seed), then the
        # frequency-ranked top-m of THAT budget: exactly matched, uses only observed data
        un_errs, tm_errs, tms = [], [], []
        for sd in range(nseeds):
            st = streams_all[sd]
            if len(observed_union(st, smax, dominant)) < m:
                continue
            lo, hi = 0, smax
            while hi - lo > 1:
                mid = (lo + hi) // 2
                if len(observed_union(st, mid, dominant)) >= m: hi = mid
                else: lo = mid
            tms.append(int(hi * (K_SLICE + 1)))
            un_errs.append(rel_l1(P, observed_union(st, hi, dominant)))
            tm_errs.append(rel_l1(P, observed_topm(st, hi, dominant, m)))
        if not tm_errs: continue
        # how many shots does it take to BUY BACK the infinite-shot ordering at fixed |S|?
        budget_curve = []
        s_min = int(np.mean(tms) / (K_SLICE + 1))
        for mult in (1, 2, 4, 8, 16):
            ss = int(min(smax, max(1, s_min * mult)))
            e = []
            for sd in range(nseeds):
                st = streams_all[sd]
                if len(observed_union(st, ss, dominant)) < m: continue
                e.append(rel_l1(P, observed_topm(st, ss, dominant, m)))
            if not e: continue
            budget_curve.append({'multiple_of_minimum': mult, 'shots_per_slice': ss,
                                 'T': ss * (K_SLICE + 1),
                                 'relL1_mean': float(np.mean(e)), 'relL1_std': float(np.std(e)),
                                 'penalty_factor': float(np.mean(e)) / or_err,
                                 'n_seeds_used': len(e)})
            log(f"L={L}   |S|={m} fixed, T={ss*(K_SLICE+1):>9d}: sampled "
                f"{np.mean(e):.4f}+-{np.std(e):.4f}  penalty x{np.mean(e)/or_err:.2f}")
            if ss >= smax: break
        matched.append({'S_matched': int(m), 'fraction': m / nS,
                        'budget_dependence': budget_curve,
                        'oracle_relL1': or_err,
                        'sampled_topm_relL1_mean': float(np.mean(tm_errs)),
                        'sampled_topm_relL1_std': float(np.std(tm_errs)),
                        'sampled_union_relL1_mean': float(np.mean(un_errs)),
                        'sampled_union_relL1_std': float(np.std(un_errs)),
                        'penalty_factor_topm': float(np.mean(tm_errs)) / or_err,
                        'penalty_factor_union': float(np.mean(un_errs)) / or_err,
                        'T_to_reach_mean': float(np.mean(tms)),
                        'n_seeds_used': len(tm_errs)})
        log(f"L={L} matched |S|={m} (frac {m/nS:.4f}): oracle {or_err:.4f} vs sampled "
            f"{np.mean(tm_errs):.4f}+-{np.std(tm_errs):.4f}  -> x{np.mean(tm_errs)/or_err:.2f}")

    del streams_all
    return dict(
        L=L, qubits=2 * L, dim_Np1=nS, E0=P['E0'], sumrule=P['sumrule'],
        integral_A_exact=P['integral'],
        oracle={'fraction': frac_or, 'S': k_or, 'published_reference': ref,
                'relL1_at_S': rel_l1(P, order[:k_or])},
        shot_sweep=sweep,
        min_T={'per_seed_T': per_seed_T, 'per_seed_frac': per_seed_frac,
               'T_mean': float(np.mean(good_T)) if good_T else None,
               'T_std': float(np.std(good_T)) if good_T else None,
               'T_min': int(np.min(good_T)) if good_T else None,
               'T_max': int(np.max(good_T)) if good_T else None,
               'frac_mean': float(np.mean(good_f)) if good_f else None,
               'frac_std': float(np.std(good_f)) if good_f else None,
               'n_censored': censored, 'T_budget_ceiling': int(smax * (K_SLICE + 1))},
        matched_subspace=matched,
        expected_distinct=edist,
        T_to_match_published_S=(None if T_match_published is None
                                else int(T_match_published * (K_SLICE + 1))),
    ), P


def controls_dense(P):
    """Control C2/C3: dense Lehmann cross-check of the Haydock reference."""
    L = P['L']; nS = P['nS']
    H = P['H1'].toarray()
    En, Vn = np.linalg.eigh(H)
    amp = Vn.conj().T @ P['seed']
    poles = En - P['E0']; wts = np.abs(amp) ** 2
    g = P['grid']
    A = np.zeros_like(g)
    for a, b in zip(poles, wts):
        A += b * (ETA / np.pi) / ((g - a) ** 2 + ETA ** 2)
    err = float(_trapz(np.abs(A - P['A_ex']), g) / P['nrm'])
    log(f"CONTROL C2  L={L}: Haydock(nl={min(NL,nS)}) vs dense Lehmann  rel-L1={err:.3e}")
    if err > 1e-6:
        raise SystemExit(f"ABORT: control C2 failed at L={L}: rel-L1={err:.3e} > 1e-6")
    return err


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--L', type=int, nargs='+', default=[4, 6, 8])
    ap.add_argument('--out', default=DEFAULT_OUT,
                    help='output JSON (default: <repo>/data/sampled_honest.json); '
                         'a relative path is taken relative to the CURRENT DIRECTORY, '
                         'never to the repository root')
    ap.add_argument('--nseeds', type=int, default=NSEEDS)
    ap.add_argument('--smax-cap', type=int, default=300000)
    ap.add_argument('--sweep-points', type=int, default=7)
    ap.add_argument('--dense-check-max-L', type=int, default=6)
    args = ap.parse_args()

    results, c2 = [], {}
    for L in args.L:
        res, P = run_L(L, args.nseeds, args.smax_cap, args.sweep_points)
        # C3: window loss on the sum rule
        loss = abs(P['integral'] - P['sumrule']) / P['sumrule']
        log(f"CONTROL C3  L={L}: int A_ex={P['integral']:.6f} vs ||seed||^2="
            f"{P['sumrule']:.6f}  window loss={100*loss:.2f}%")
        if loss > 0.06:
            raise SystemExit(f"ABORT: control C3 failed at L={L}: window loss {100*loss:.2f}%")
        res['window_loss_rel'] = float(loss)
        if L <= args.dense_check_max_L:
            c2[L] = controls_dense(P)
        results.append(res)
        del P

    out = {
        'provenance': {
            'script': os.path.basename(__file__),
            'purpose': 'finite-shot (honest) replacement for the argsort-on-exact-Born '
                       'subspace selection of scaling_lanczos.py:66 / sampled_akw.py:47',
            'generated_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
            'runtime_seconds': round(time.time() - t0, 1),
            'python': sys.version.split()[0], 'platform': platform.platform(),
            'numpy': np.__version__, 'scipy': scipy.__version__,
            'parameters': {'model': '1D Hubbard, PBC, half filling', 'U_over_t': U_HUB,
                           'eta_over_t': ETA, 'K_slices': K_SLICE, 'dt': DT,
                           'threshold_relL1': THR, 'lanczos_depth': NL, 'ngrid': NGRID,
                           'n_seeds': args.nseeds, 'seed_base': SEED_BASE,
                           'seeds': [SEED_BASE + i for i in range(args.nseeds)],
                           'operator': 'c^dag_{0,up} |GS>  in (N+1, Sz+1)',
                           'shot_allocation': 'T split evenly over the K+1 = %d time slices'
                                              % (K_SLICE + 1),
                           'sampler': 'inverse-CDF on |psi(t_k)|^2, ordered stream, '
                                      'prefixes = smaller budgets'},
            'controls': {
                'C1_oracle_fraction_reference': ORACLE_REF,
                'C1_tolerance_abs': ORACLE_TOL,
                'C2_haydock_vs_dense_lehmann_relL1': c2,
                'C2_tolerance': 1e-6,
                'C3_sum_rule_window_loss_tolerance': 0.06,
            },
        },
        'results': results,
    }
    # ---- exponential extrapolation of the shot budget (exact-formula anchors only) ----
    # L whose infinite-shot subspace already spans the whole Born support (rel-L1 == 0
    # there) is degenerate: its "required fraction" is set by the support, not by
    # accuracy, so it is excluded from every exponent fit and flagged.
    degen = [r['L'] for r in results if r['oracle']['relL1_at_S'] < 1e-10]
    use = [r for r in results if r['L'] not in degen]
    def expfit(xs, ys):
        if len(xs) < 2: return None
        b, a0 = np.polyfit(np.asarray(xs, float), np.log(np.asarray(ys, float)), 1)
        return {'slope_per_site': float(b), 'intercept': float(a0),
                'growth_per_two_sites': float(np.exp(2 * b)),
                'n_anchors': len(xs),
                'extrapolation': {int(L): float(np.exp(a0 + b * L)) for L in (10, 12, 14)}}
    Ls = [r['L'] for r in use]
    out['shot_budget_extrapolation'] = {
        'degenerate_L_excluded': degen,
        'anchors_L': Ls,
        'T_to_match_published_S': [r['T_to_match_published_S'] for r in use],
        'T_star_relL1_below_threshold': [r['min_T']['T_mean'] for r in use],
        'S_published': [r['oracle']['S'] for r in use],
        'S_star_honest': [r['min_T']['frac_mean'] * r['dim_Np1'] for r in use],
        'dim_Np1': [r['dim_Np1'] for r in use],
        'fit_T_match_vs_L': expfit(Ls, [r['T_to_match_published_S'] for r in use]),
        'fit_T_star_vs_L': expfit(Ls, [r['min_T']['T_mean'] for r in use]),
        'fit_S_star_vs_L': expfit(Ls, [r['min_T']['frac_mean'] * r['dim_Np1'] for r in use]),
        'fit_dim_vs_L': expfit(Ls, [r['dim_Np1'] for r in use]),
        'caveat': 'ONLY %d non-degenerate anchors: the exponents are indicative, not '
                  'determined. The L>=10 Born distributions were NOT computed here; '
                  'obtaining them needs the matrix-free ranking stage of '
                  'scaling_lanczos_mf.py.' % len(Ls),
    }

    outp = os.path.abspath(args.out)
    os.makedirs(os.path.dirname(outp), exist_ok=True)
    with open(outp, 'w') as f:
        json.dump(out, f, indent=1)
    log(f"WROTE {outp}")


if __name__ == '__main__':
    main()
