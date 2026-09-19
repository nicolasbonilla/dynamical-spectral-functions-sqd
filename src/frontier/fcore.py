# -*- coding: utf-8 -*-
r"""C3 core: the certificate frontier, measured.  MEMORY-LEAN.

Reuses (never rewrites):
  c0_lib.py     -> lanczos(reorth), born_order_mf, AK (src/akw_lanczos.py)
  c0_big.py     -> sector_mv (matrix-free Hubbard sector)
  c0_sweep.py   -> lambda_from_gram (closed form of Lambda(eta))
  (all three deposited beside this file in src/frontier/)

THE ONE NEW PIECE OF ALGEBRA (it is what makes C3 affordable at all):
  With full reorthogonalisation, V^T H V = T (tridiagonal, from alpha/beta) and
      H_S V = V T + beta_n v_{n+1} e_n^T ,        v_{n+1} _|_ V
  hence, with G1 := V^T H^2 V,
      Gram_K  = Sm^T [ G1 - T^2 ] Sm                          (= (Q_K H V)^T (Q_K H V))
      Gram_Ki = beta_n^2 * outer( Sm[n-1,:], Sm[n-1,:] )      RANK ONE, free from (alpha,beta)
      Gram_Ko = Gram_K - Gram_Ki
  so  HV IS NEVER STORED (only V), and every prefix depth n <= n_l is free once G1 is known,
  because G1_n / T_n / beta_n are just leading blocks.  G1 is built column by column with
  2 matvecs and ONE temporary vector.

  Consequences used below, all exact (not approximations):
    P_K phi = phi_S  =>  w_K = w_S exactly at every depth;  V^T phi = |phi_S| e_1.
    Lambda_K^2 = Lambda_K,in^2 + Lambda_K,out^2   (lambda^2 is linear in the Gram).
"""
import os, sys, time, gc, json, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PUERTAS = HERE          # c0_lib / c0_big / c0_sweep are deposited beside this file
sys.dont_write_bytecode = True
sys.path.insert(0, PUERTAS)
import c0_lib as C                     # noqa: E402
from c0_lib import AK                  # noqa: E402
from c0_big import sector_mv           # noqa: E402
from c0_sweep import lambda_from_gram  # noqa: E402
from scipy.sparse.linalg import eigsh, LinearOperator  # noqa: E402

T0 = time.time()
def log(*a):
    print("[%7.1fs]" % (time.time() - T0), *a, flush=True)

# ---- repo protocol, verbatim (release/src/scaling_lanczos_mf.py:90) --------
PROTO = dict(K=18, dt=0.5, sub=16, m=6)     # 19 snapshots, t = 0.0 .. 9.0
U_DEF = 8.0


# ===========================================================================
#  system
# ===========================================================================
def build_system(L, U=U_DEF):
    """Ground state of the half-filled sector, seed phi = c^dag_{0,up}|Psi_0>,
    matrix-free matvec of the (N+1) sector.  Identical conventions to the repo."""
    nup = nd = L // 2
    t = time.time()
    mv0, D0, Du0, Dd0 = sector_mv(L, U, nup, nd)
    E0, V0 = eigsh(LinearOperator((D0, D0), matvec=mv0, dtype=float),
                   k=1, which='SA', ncv=12, tol=1e-11)
    E0 = float(E0[0]); Psi = V0[:, 0].reshape(Du0, Dd0); del V0, mv0
    cdU = AK.cdag_map(L, nup)
    phi = np.asarray((cdU[0] @ Psi).reshape(-1), dtype=float)
    del Psi, cdU; gc.collect()
    mv1, D, _, _ = sector_mv(L, U, nup + 1, nd)
    return dict(L=L, U=U, E0=E0, phi=phi, mv=mv1, D=D, t_gs=time.time() - t)


def ranking(sysd, K=None, dt=None, sub=None, mkry=None):
    """Accumulated Born weight of the TIME-EVOLVED seed -> ranking.
    Defaults are the REPO protocol K=18, dt=0.5, sub=16, m=6 (19 snapshots)."""
    K = PROTO['K'] if K is None else K
    dt = PROTO['dt'] if dt is None else dt
    sub = PROTO['sub'] if sub is None else sub
    mkry = PROTO['m'] if mkry is None else mkry
    mv = sysd['mv']

    class _H:
        def __matmul__(self, x): return mv(x)
    t = time.time()
    order, wc = C.born_order_mf(_H(), sysd['phi'], K=K, dt=dt, sub=sub, m=mkry)
    return order, wc, time.time() - t


def subspace(order, wc, D, FR):
    """Top-FR*D configs.  Also returns the zero-amplitude / tie diagnostics that
    condition 5 of the gate verdict demands be DECLARED."""
    k = max(1, int(round(FR * D)))
    sel = order[:k]
    mask = np.zeros(D, bool); mask[sel] = True
    wsel = wc[sel]
    cut = float(wc[order[k - 1]])
    n_zero = int(np.count_nonzero(wsel == 0.0))
    n_at_cut_in = int(np.count_nonzero(wsel == cut))
    n_at_cut_tot = int(np.count_nonzero(wc == cut))
    diag = dict(k=int(k), FR_real=k / D, cut_weight=cut,
                n_zero_wc_in_S=n_zero,
                n_tied_at_cut_inside=n_at_cut_in,
                n_tied_at_cut_total=n_at_cut_tot,
                ambiguous=bool(n_at_cut_tot > n_at_cut_in),
                min_wc_in_S=float(wsel.min()),
                max_wc_out=float(wc[~mask].max()) if k < D else 0.0)
    return mask, diag


# ===========================================================================
#  the lean Lambda machinery
# ===========================================================================
def lanczos_S(mv, mask, phi, nl):
    """Full-reorthogonalisation Lanczos on H_S = P_S H P_S seeded with phi_S.
    CONDITION 1 of the gate verdict: reorth is NOT optional."""
    phiS = phi * mask

    def mvS(x):
        return mask * mv(mask * x)
    V, al, be, nrm = C.lanczos(mvS, phiS.astype(float), nl, reorth=True)
    return V, al, be, nrm, phiS


def gram_G1(mv, V, chunk_log=None):
    """G1 = V^T H^2 V, column by column: 2 matvecs + ONE temporary vector.
    This is the whole memory trick -- HV (a second D x n_l block) is never formed."""
    m = V.shape[1]
    G1 = np.empty((m, m))
    for j in range(m):
        G1[:, j] = V.T @ mv(mv(V[:, j]))
        if chunk_log and (j + 1) % chunk_log == 0:
            log(f"    G1 column {j+1}/{m}")
    return 0.5 * (G1 + G1.T)


def depth_eval(al, be, G1, n, nphiS2, etas, keep_pairs=False):
    """Everything at Krylov depth n (n <= n_l), free from leading blocks."""
    T = np.diag(al[:n]) + np.diag(be[:n - 1], 1) + np.diag(be[:n - 1], -1)
    th, Sm = np.linalg.eigh(T)
    c = np.sqrt(nphiS2) * Sm[0, :]                   # V^T phi = |phi_S| e_1  (exact)
    GK = Sm.T @ ((G1[:n, :n] - T @ T) @ Sm); GK = 0.5 * (GK + GK.T)
    last = Sm[n - 1, :]
    GKi = (be[n - 1] ** 2) * np.outer(last, last)    # rank one, exact
    GKo = GK - GKi
    out = dict(n=int(n), beta_n=float(be[n - 1]))
    if keep_pairs:
        out['theta'] = th; out['weights'] = c ** 2
    per = {}
    cc = c.astype(complex)
    for eta in etas:
        lk, _ = lambda_from_gram(th, cc, GK, eta)
        li, _ = lambda_from_gram(th, cc, GKi, eta)
        lo, _ = lambda_from_gram(th, cc, GKo, eta)
        per[float(eta)] = (float(lk), float(li), float(lo))
    out['per_eta'] = per
    return out


def certificate(wS, lamhat_over_eta):
    """Theorem B, relative form:  rel-L1 <= (1+sqrt(w))*(sqrt(1-w) + Lambda_hat/eta).
    Trivial bound 1+w.  The certificate is the min of the two."""
    sw = np.sqrt(wS)
    B = (1.0 + sw) * (np.sqrt(max(1.0 - wS, 0.0)) + lamhat_over_eta)
    triv = 1.0 + wS
    return float(B), float(triv), float(min(B, triv))


# ===========================================================================
#  rel-L1 between two Lorentzian sums (exact tails)
# ===========================================================================
def _lorentz(a, w, grid, eta, blk=2000):
    A = np.zeros_like(grid)
    for i0 in range(0, len(a), blk):
        aa = a[i0:i0 + blk]; ww = w[i0:i0 + blk]
        A += (eta / np.pi) * (ww[:, None] / ((grid[None, :] - aa[:, None]) ** 2 + eta ** 2)).sum(0)
    return A


def l1_lorentz(a1, w1, a2, w2, eta, pad=8.0, ppw=24.0):
    """int |A1 - A2| dw over the WHOLE real line: fine core grid + analytic tails."""
    lo = float(min(a1.min(), a2.min())) - pad
    hi = float(max(a1.max(), a2.max())) + pad
    npts = int((hi - lo) / (eta / ppw)) + 1
    npts = min(max(npts, 4000), 400000)
    grid = np.linspace(lo, hi, npts)
    A1 = _lorentz(a1, w1, grid, eta); A2 = _lorentz(a2, w2, grid, eta)
    core = float(np.trapz(np.abs(A1 - A2), grid))
    Rh1 = 0.5 - np.arctan((hi - a1) / eta) / np.pi
    Rh2 = 0.5 - np.arctan((hi - a2) / eta) / np.pi
    Rl1 = 0.5 + np.arctan((lo - a1) / eta) / np.pi
    Rl2 = 0.5 + np.arctan((lo - a2) / eta) / np.pi
    tail = abs(float(w1 @ Rh1 - w2 @ Rh2)) + abs(float(w1 @ Rl1 - w2 @ Rl2))
    return dict(l1=core + tail, core=core, tail=tail,
                int_A1=float(np.trapz(A1, grid)) + float(w1 @ Rh1) + float(w1 @ Rl1),
                win_absA1=float(np.trapz(np.abs(A1), grid)))


def ref_spectrum(mv, phi, nl_ref, D):
    """Exact A(w): full-reorthogonalisation Lanczos on the FULL sector, seed phi.
    Gauss quadrature exact in the first 2*nl_ref-1 moments."""
    V, al, be, nrm = C.lanczos(mv, phi.astype(float), min(nl_ref, D), reorth=True)
    n = V.shape[1]
    T = np.diag(al[:n]) + np.diag(be[:n - 1], 1) + np.diag(be[:n - 1], -1)
    th, Sm = np.linalg.eigh(T)
    w = (nrm ** 2) * Sm[0, :] ** 2
    del V; gc.collect()
    return th, w, n


def ladder(nl, base=(20, 30, 40, 60, 80, 100, 130, 160, 200, 250, 300, 375, 450,
                     550, 650, 800, 1000, 1250, 1500)):
    out = [n for n in base if n < nl]
    out.append(int(nl))
    return out
