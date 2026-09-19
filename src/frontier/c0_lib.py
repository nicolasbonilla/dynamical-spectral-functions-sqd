# -*- coding: utf-8 -*-
"""C0 gate: Theorem B for a general ORTHOGONAL projector.

Read-only w.r.t. the repo: we only IMPORT src/akw_lanczos.py.  Everything else is
written here.  Deposited 2026-09-19 from the C0/C3 working tree; the only change is
the SRC path below, which was an absolute Windows path.

Conventions match the repo exactly:
  H       : (N+1) sector Hubbard, PBC ring, U=8, t=1, half filling  (AK.build_H_explicit)
  phi     : c^dag_{0,up}|Psi_0>   (scaling_lanczos_mf.stage_exact, line 88)
  S       : top-|S| configs by accumulated Born weight of the TIME-EVOLVED seed
            (sampled_akw.channel lines 42-47:  K=16, dt=0.5, wc = sum_k |e^{-iHk dt}v|^2)
  A(w)    = (eta/pi) || (w+E0+i eta - H)^{-1} phi ||^2
  A_P(w)  = (eta/pi) || (w+E0+i eta - P H P)^{-1} P phi ||^2     (Rayleigh-Ritz on range(P))
"""
import os, sys, time, numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # release/src
sys.path.insert(0, SRC)
import akw_lanczos as AK

T0 = time.time()
def log(*a): print("[%7.1fs]" % (time.time()-T0), *a, flush=True)

# ---------------------------------------------------------------- system
def system(L, U=8.0, t=1.0):
    """ground state of the half-filled sector, seed phi in the (N+1) sector, H there."""
    nup = nd = L//2
    H0, Du, Dd = AK.build_H_explicit(L, U, nup, nd, t)
    E0, V0 = eigsh(H0, k=1, which='SA'); E0 = float(E0[0])
    Psi = V0[:, 0].reshape(Du, Dd)
    cdU = AK.cdag_map(L, nup)                     # c^dag_j : nup -> nup+1
    phi = np.asarray((cdU[0] @ Psi).reshape(-1), dtype=float)
    H1, Du1, Dd1 = AK.build_H_explicit(L, U, nup+1, nd, t)
    return dict(L=L, U=U, E0=E0, phi=phi, H=H1.tocsr(), D=H1.shape[0])

def born_order(H, phi, K=16, dt=0.5, dense=True):
    """accumulated Born weight of the time-evolved seed -> ranking (sampled_akw protocol)."""
    v0 = phi/np.linalg.norm(phi)
    if dense:
        Ha = H.toarray(); Em, Um = np.linalg.eigh(Ha)
        coef = Um.conj().T @ v0
        wc = np.zeros(H.shape[0])
        for k in range(K+1):
            vk = Um @ (np.exp(-1j*Em*k*dt)*coef); wc += np.abs(vk)**2
    else:
        raise NotImplementedError
    return np.argsort(wc)[::-1], wc

# ---------------------------------------------------------------- spectra
def A_exact(H, phi, E0, grid, eta):
    """(eta/pi)||R(z)phi||^2 by full eigendecomposition (small systems only)."""
    Em, Um = np.linalg.eigh(H.toarray() if sp.issparse(H) else H)
    c = Um.conj().T @ phi
    w = np.abs(c)**2
    A = np.zeros_like(grid)
    for a, b in zip(Em-E0, w):
        A += b*(eta/np.pi)/((grid-a)**2 + eta**2)
    return A

def A_from_pairs(theta, weights, E0, grid, eta):
    A = np.zeros_like(grid)
    for a, b in zip(theta-E0, weights):
        A += b*(eta/np.pi)/((grid-a)**2 + eta**2)
    return A

# ---------------------------------------------------------------- Lanczos
def lanczos(matvec, v0, nl, reorth=True):
    """Lanczos on a Hermitian operator. Returns V (n x m, columns orthonormal),
    alpha, beta.  reorth=True -> full reorthogonalization (exact P_K = V V^dag).
    reorth=False -> the plain 3-term recursion the repo's AK.haydock uses."""
    n = v0.shape[0]
    nrm = np.linalg.norm(v0)
    V = np.zeros((n, nl), dtype=v0.dtype)
    V[:, 0] = v0/nrm
    al = np.zeros(nl); be = np.zeros(nl)
    m = nl
    w = matvec(V[:, 0]); a = float(np.vdot(V[:, 0], w).real); al[0] = a
    w = w - a*V[:, 0]
    for j in range(1, nl):
        if reorth:                       # two-pass classical Gram-Schmidt
            for _ in range(2):
                Vj = V[:, :j]
                w -= Vj @ ((Vj.conj().T if np.iscomplexobj(Vj) else Vj.T) @ w)
        b = np.linalg.norm(w); be[j-1] = b
        if b < 1e-12:
            m = j; break
        V[:, j] = w/b
        w = matvec(V[:, j]); a = float(np.vdot(V[:, j], w).real); al[j] = a
        w = w - a*V[:, j] - b*V[:, j-1]
    if m == nl:
        if reorth:
            for _ in range(2):
                Vj = V[:, :nl]
                w -= Vj @ ((Vj.conj().T if np.iscomplexobj(Vj) else Vj.T) @ w)
        be[nl-1] = np.linalg.norm(w)      # beta_{n_l}: the overflow out of K
    return V[:, :m], al[:m], be[:m], nrm

def tridiag(al, be):
    m = len(al)
    return np.diag(al) + np.diag(be[:m-1], 1) + np.diag(be[:m-1], -1)

def haydock_cf(al, be, nrm2, z):
    """continued fraction <v0|(z-T)^{-1}|v0> * nrm2  (identical to AK.haydock)."""
    m = len(al)
    G = z - al[-1]
    for k in range(m-2, -1, -1):
        G = (z - al[k]) - be[k]**2/G
    return nrm2/G

# ---------------------------------------------------------------- Lambda
def F_kernel(Ep, Em, eta):
    """[4 eta^2 + 2 i eta (E' - E)] / [(E'-E)^2 + 4 eta^2], indexed [m', m]."""
    d = Ep[:, None] - Em[None, :]
    return (4*eta**2 + 2j*eta*d)/(d**2 + 4*eta**2)

def lambda_from_pairs(theta, c, G, eta, block=512):
    """Lambda^2 = sum_{m,m'} c_m conj(c_m') <g_m', g_m> F[m',m]; G columns are g_m."""
    nm = len(theta)
    tot = 0.0 + 0.0j
    cc = np.conj(c)
    for i0 in range(0, nm, block):
        i1 = min(i0+block, nm)
        Gram = G[:, i0:i1].conj().T @ G                      # (blk, nm) = <g_m', g_m>
        d = theta[i0:i1, None] - theta[None, :]
        F = (4*eta**2 + 2j*eta*d)/(d**2 + 4*eta**2)
        tot += np.sum(cc[i0:i1, None]*c[None, :]*Gram*F)
    lam2 = float(tot.real)
    return np.sqrt(max(lam2, 0.0)), lam2, float(tot.imag)

# ---------------------------------------------------------------- matrix-free helpers (L>=10)
_KEXP_WS = {}
def krylov_expm(mv, v, dt, m=6):
    """e^{-i dt H} v via a small Lanczos basis (copied verbatim in spirit from
    release/src/scaling_lanczos_mf.py:krylov_expm, so the ranking protocol is identical)."""
    from scipy.linalg import expm as dense_expm
    n = v.shape[0]
    key = (m, n)
    V = _KEXP_WS.get(key)
    if V is None:
        V = np.empty((m, n), dtype=complex); _KEXP_WS[key] = V   # reused across calls (no alloc churn)
    beta = np.linalg.norm(v)
    if beta < 1e-14: return v.copy()
    V[0] = v/beta; al = np.zeros(m); be = np.zeros(m)
    w = mv(V[0]); a = np.vdot(V[0], w).real; al[0] = a; w = w-a*V[0]; mm = m
    for j in range(1, m):
        b = np.linalg.norm(w); be[j-1] = b
        if b < 1e-12: mm = j; break
        V[j] = w/b; w = mv(V[j]); a = np.vdot(V[j], w).real; al[j] = a; w = w-a*V[j]-b*V[j-1]
    T = np.diag(al[:mm])+np.diag(be[:mm-1], 1)+np.diag(be[:mm-1], -1)
    E = dense_expm(-1j*dt*T)[:, 0]
    out = (beta*(V[:mm].T@E)); return out

def born_order_mf(H, phi, K=16, dt=0.5, sub=16, m=6):
    """matrix-free accumulated Born weight ranking (scaling_lanczos_mf.stage_rank protocol)."""
    mv = lambda x: H@x
    v = (phi/np.linalg.norm(phi)).astype(complex)
    wc = (np.abs(v)**2).astype(np.float32)
    ddt = dt/sub
    for s in range(1, K*sub+1):
        v = krylov_expm(mv, v, ddt, m)
        if s % sub == 0: wc += (np.abs(v)**2).astype(np.float32)
    return np.argsort(wc)[::-1], wc
