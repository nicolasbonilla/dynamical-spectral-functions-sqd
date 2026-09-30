# -*- coding: utf-8 -*-
r"""Static determinant support of the free-fermion (U = 0) ground state of the open chain.

Generator of the four values printed in Sec. S13 (sec:nogo:bilateral) and Sec. VI:
    |S|_{eps_w} = 35, 336, 3496, 37361   for L = 4, 6, 8, 10,
open chain, nearest-neighbour hopping t = 1, half filling (N_up = N_dn = L/2), eps_w = 2.5e-3.
Added 2026-09-28: until then these four values had no deposited generator.

DEFINITION (the one of Sec. S13 and of src/cost_vs_entanglement.py)
--------------------------------------------------------------------
|S|_{eps_w} is the number of site-basis Slater determinants of the ground state, taken in order
of decreasing weight, needed to carry all but eps_w of its weight: the smallest n with
    sum_{i<=n} w_(i) >= 1 - eps_w          (inclusive: a set whose discarded weight is EXACTLY
                                             eps_w counts as carrying all but eps_w).
cost_vs_entanglement.py writes the same cut as np.searchsorted(cumsum(w), 1 - eps^2) + 1 with
eps = 0.05, i.e. 1 - eps_w = 0.9975.

THE FREE-FERMION GROUND STATE
-----------------------------
At U = 0 the ground state is the product of one up-spin and one down-spin Slater determinant
built from the L/2 lowest orbitals phi_m(j) = sqrt(2/(L+1)) sin(pi m j/(L+1)) of the open chain
(non-degenerate for even L, so the half-filled ground state is unique).  The weight of the site
configuration (I_up, I_dn) is det(Phi[I_up])^2 det(Phi[I_dn])^2.  The weights are computed from
the determinants directly, not from a many-body diagonalization; for L = 4 and 6 the script also
cross-checks them against the exact diagonalization of src/gate1_ladder.py at U = 0.

THE TIE AT L = 4
----------------
At L = 4 the four smallest determinant weights are EXACTLY 1/400 = eps_w each (stored under
`smallest_weights`).  The inclusive cut therefore stops at 35 determinants, discarding one of
the four, with discarded weight exactly eps_w; a strict cut (discarded weight < eps_w) gives 36.  In
floating point the cumulative sum lands a few ulp either side of 0.9975, so a bare
np.searchsorted gives 35 or 36 depending on the summation order.  This script therefore
applies the inclusive cut with a relative tolerance of 1e-12 on the discarded weight and
records the strict count as well.  For L = 6, 8, 10 the discarded weight at the cut is strictly
below eps_w and both cuts agree.  (`n_weights_tied_at_boundary` counts the weights equal, by
symmetry, to the last one kept: it makes the kept SET non-unique but never the COUNT.)

Run:  python src/free_fermion_support.py      -> data/free_fermion_support.json  (< 1 min)
"""
import os, sys, json, itertools, hashlib, datetime, platform
import numpy as np

sys.dont_write_bytecode = True
_SRC = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_SRC)
sys.path.insert(0, _SRC)

EPS_W = 2.5e-3          # the lattice threshold of Sec. S13 (1 - eps^2 with eps = 0.05)
RTOL = 1e-12            # tie tolerance on the discarded weight
SIZES = (4, 6, 8, 10)


def orbitals(L):
    """Columns = open-chain single-particle orbitals, ascending energy -2 cos(pi m/(L+1))."""
    j = np.arange(1, L + 1)
    m = np.arange(1, L + 1)
    Phi = np.sqrt(2.0 / (L + 1)) * np.sin(np.pi * np.outer(j, m) / (L + 1))
    eps = -2.0 * np.cos(np.pi * m / (L + 1))
    order = np.argsort(eps)
    return Phi[:, order], eps[order]


def spin_weights(L, n):
    """det(Phi[I, :n])^2 for every n-subset I of the L sites (lexicographic order)."""
    Phi, _ = orbitals(L)
    occ = Phi[:, :n]
    subsets = list(itertools.combinations(range(L), n))
    w = np.array([np.linalg.det(occ[list(I), :]) ** 2 for I in subsets])
    return subsets, w


def support(L, eps_w=EPS_W):
    n = L // 2
    _, wu = spin_weights(L, n)
    w = np.sort(np.outer(wu, wu).ravel())[::-1]        # N_up = N_dn = L/2: same weights
    total = float(w.sum())
    # discarded weight after keeping the first k determinants, summed from the SMALL end
    # (tail sums of small numbers are accurate; 1 - cumsum is not)
    tail = np.concatenate([np.cumsum(w[::-1])[::-1], [0.0]])     # tail[k] = sum_{i>=k} w_i
    tail /= total
    thr = eps_w * (1.0 + RTOL)
    n_incl = int(np.argmax(tail <= thr))                          # smallest k with tail <= eps_w
    n_strict = int(np.argmax(tail < eps_w * (1.0 - RTOL)))        # smallest k with tail < eps_w
    wb = w[n_incl - 1] / total if n_incl > 0 else float("nan")
    tied = int(np.sum(np.abs(w / total - wb) <= RTOL * wb))
    # the bare float convention of cost_vs_entanglement.py, for the record
    c = np.cumsum(w / total)
    n_float = int(np.searchsorted(c, 1 - 0.05 * 0.05) + 1)
    return dict(L=L, N_up=n, N_dn=n, dim=int(w.size), total_weight=total,
                support_inclusive=n_incl, support_strict=n_strict,
                discarded_weight_inclusive=float(tail[n_incl]),
                boundary_weight=float(wb), n_weights_tied_at_boundary=tied,
                support_bare_float_searchsorted=n_float,
                smallest_weights=[float(x) for x in (w[-4:] / total)])


def ed_crosscheck(L):
    """Max |w_det - w_ED| against the U=0 exact diagonalization of gate1_ladder.py (L<=6)."""
    import gate1_ladder as G
    from scipy.sparse.linalg import eigsh
    n = L // 2
    H, Su, iu, Sd, idd, Du, Dd = G.build_H(L, 0.0, n, n, [(i, i + 1, 1.0) for i in range(L - 1)])
    e, v = eigsh(H, k=2, which='SA', v0=np.random.default_rng(20260918).standard_normal(H.shape[0]), tol=0)
    gap = float(sorted(e)[1] - sorted(e)[0])
    psi = v[:, int(np.argmin(e))]
    wED = np.sort(np.abs(psi) ** 2)[::-1]
    _, wu = spin_weights(L, n)
    wdet = np.sort(np.outer(wu, wu).ravel())[::-1]
    return dict(L=L, max_abs_weight_diff=float(np.max(np.abs(wED / wED.sum() - wdet / wdet.sum()))),
                many_body_gap=gap)


if __name__ == "__main__":
    rows = [support(L) for L in SIZES]
    for r in rows:
        print("L=%2d  dim=%6d  |S|_eps_w = %6d (inclusive)  %6d (strict)  tied at boundary: %d  "
              "discarded %.15g  bare-float searchsorted: %d" % (
                  r['L'], r['dim'], r['support_inclusive'], r['support_strict'],
                  r['n_weights_tied_at_boundary'], r['discarded_weight_inclusive'],
                  r['support_bare_float_searchsorted']))
    checks = []
    for L in (4, 6):
        try:
            checks.append(ed_crosscheck(L))
            print("ED cross-check L=%d: max |w_det - w_ED| = %.2e, many-body gap %.4f" % (
                L, checks[-1]['max_abs_weight_diff'], checks[-1]['many_body_gap']))
        except Exception as exc:                                  # pragma: no cover
            checks.append(dict(L=L, error=repr(exc)))
            print("ED cross-check L=%d not run: %r" % (L, exc))
    with open(os.path.abspath(__file__), "rb") as f:
        h = hashlib.sha256(f.read().replace(b"\r\n", b"\n")).hexdigest()
    out = dict(
        what="static determinant support |S|_eps_w of the U=0 (free-fermion) ground state, open chain, "
             "half filling, site basis; printed in Sec. S13 and Sec. VI",
        definition="smallest n such that the n largest determinant weights sum to >= 1 - eps_w (inclusive), "
                   "tie tolerance %g relative on the discarded weight" % RTOL,
        eps_w=EPS_W, hopping="t=1, open boundary conditions", rows=rows, ed_crosscheck=checks,
        published=dict(zip([str(L) for L in SIZES], [r['support_inclusive'] for r in rows])),
        provenance=dict(script="src/free_fermion_support.py", generator_sha256_lf=h,
                        generated_utc=datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                        python=platform.python_version(), numpy=np.__version__))
    p = os.path.join(_REPO, "data", "free_fermion_support.json")
    json.dump(out, open(p, "w"), indent=1)
    print("wrote", p)
