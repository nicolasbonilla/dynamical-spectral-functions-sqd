# -*- coding: utf-8 -*-
"""C6 / frontier instrumentation library.

Instruments the FOUR steps of the Theorem-B proof (VE section 1.5) and measures
where the slack is lost.  Everything EXACT (dense eigendecomposition of H and of
H_S = P_S H P_S) -- NO Krylov anywhere, so condition (1) full reorthogonalisation
and condition (3) n_l convergence are satisfied VACUOUSLY (there is no n_l), and
condition (4) is satisfied because what is certified is A_S itself, not A_K.

Ranking protocol: REPO protocol K=18, dt=0.5, sub=16 -> 19 snapshots (condition 2),
   release/src/scaling_lanczos_mf.py:90  stage_rank(L,U,K=18,dt=0.5,sub=16).

THE LADDER  (each rung is the value of the bound after one more step of the proof)
   R0 = int |A - A_S| dw                                  <- the truth
   R1 = int (eta/pi)(|u|+|u_S|) |u-u_S| dw                <- step (3), pointwise, BEFORE Cauchy-Schwarz
   R2 = N[F] N[G]                                         <- AFTER Cauchy-Schwarz in omega
   R3 = (|phi|+|phi_S|) N[G]                              <- AFTER Minkowski on the first factor
   R4 = (|phi|+|phi_S|)(|phi_Q| + N[R g])                 <- AFTER Minkowski on the second factor (Feshbach split)
   R5 = (|phi|+|phi_S|)(|phi_Q| + Lambda_S/eta)           <- AFTER |R(z)|<=1/eta  == THEOREM B, leakage branch
with N[f]^2 := (eta/pi) int |f(w)|^2 dw ,  F = |u|+|u_S| , G = |u-u_S| ,
     u  = R(z) phi,  u_S = R_S(z) phi_S,  z = w + E0 + i eta,
     g(w) = Q H R_S(z) phi_S,  a = R(z) phi_Q,  b = R(z) g(w),  G = a + b  (Feshbach, exact).

CLOSED FORMS USED (all verified numerically against the quadrature by the caller):
   N[u] = |phi| ,  N[u_S] = |phi_S|                        (int |R(z)psi|^2 dw = pi |psi|^2/eta)
   C_cross := (eta/pi) Re int <u,u_S> dw = sum_{m,n} c_m d_n O_mn * 4 eta^2/(D_mn^2+4 eta^2)
   N[G]^2  = |phi|^2 + |phi_S|^2 - 2 C_cross
   C_ab    := (eta/pi) Re int <a,b> dw   = sum_{m,n} q_m d_n gam_mn * D_mn/(D_mn^2+4 eta^2)
   N[b]^2  = N[G]^2 - |phi_Q|^2 - 2 C_ab
   Lambda_S^2 = Re sum_{n,n'} d_n d_n' Gram_n'n * 2 i eta/((Et_n' - Et_n) + 2 i eta)
   with c_m=<m|phi>, d_n=<nt|phi>, q_m=<m|phi_Q>, O_mn=<m|nt>, gam_mn=<m|Q H|nt>,
        D_mn = eps_m - epst_n,  eps = E - E0.
"""
import os, sys, gc, time, numpy as np, scipy.sparse as sp
sys.dont_write_bytecode = True
PUERTAS = os.path.dirname(os.path.abspath(__file__))  # c0_lib is deposited beside this file
sys.path.insert(0, PUERTAS)
import c0_lib as C          # reuse: system(), born_order(), born_order_mf(), lanczos(), A_from_pairs()
from c0_lib import AK

T0 = time.time()
def log(*a): print("[%7.1fs]" % (time.time()-T0), *a, flush=True)

# --------------------------------------------------------------- quadrature grid
def make_grid(poles, eta, per_eta=12.0, pad=None, ntail=1200, tailmax=1e8):
    """Composite trapezoid grid covering the WHOLE real line (the theorem's L1 is over R).
    core: uniform with spacing eta/per_eta;  tails: log-spaced to +-tailmax (integrands ~ C/w^2)."""
    if pad is None: pad = max(4.0, 40.0 * eta)
    lo = float(poles.min()) - pad; hi = float(poles.max()) + pad
    ncore = int(np.ceil((hi - lo) / (eta / per_eta))) + 1
    core = np.linspace(lo, hi, ncore)
    u = np.linspace(0.0, np.log(tailmax), ntail + 1)[1:]
    right = hi + np.expm1(u)
    left = lo - np.expm1(u)
    w = np.concatenate([left[::-1], core, right])
    dw = np.empty_like(w)
    dw[1:-1] = 0.5 * (w[2:] - w[:-2]); dw[0] = 0.5 * (w[1] - w[0]); dw[-1] = 0.5 * (w[-1] - w[-2])
    return w, dw, ncore

def lorentz_sum(poles, weights, w, eta, chunk=200):
    """sum_k weights_k * (eta/pi)/((w-poles_k)^2+eta^2)   -- memory-chunked over poles."""
    out = np.zeros_like(w)
    e2 = eta * eta; pref = eta / np.pi
    for i in range(0, len(poles), chunk):
        p = poles[i:i+chunk]; g = weights[i:i+chunk]
        out += ((g * pref)[:, None] / ((w[None, :] - p[:, None])**2 + e2)).sum(0)
    return out

def abs2_resolvent(poles, weights, w, eta, chunk=200):
    """|R(z)psi|^2 = sum_k weights_k/((w-poles_k)^2+eta^2)  = (pi/eta) * lorentz_sum"""
    return lorentz_sum(poles, weights, w, eta, chunk) * (np.pi / eta)

# --------------------------------------------------------------- the system + subspace
def build(L, U=8.0, t=1.0):
    s = C.system(L, U, t)
    return s

def repo_ranking(H, phi, K=18, dt=0.5, sub=16, m=6):
    """REPO protocol (scaling_lanczos_mf.stage_rank): K=18, dt=0.5, sub=16 -> 19 snapshots."""
    class _H:
        def __matmul__(self, x): return H @ x
    order, wc = C.born_order_mf(_H(), phi, K=K, dt=dt, sub=sub, m=m)
    return order, wc


# --------------------------------------------------------------- THE LADDER (exact, dense)
def ladder(Em, Um, phi, E0, H, Sset, etas, nb=128, per_eta=12.0, ntail=500, tag=""):
    """Exact 5-rung instrumentation of the Theorem-B proof for the coordinate subspace Sset."""
    D = len(phi); eps = Em - E0
    c = Um.T @ phi
    n2 = float(phi @ phi); nphi = np.sqrt(n2)
    mask = np.zeros(D, bool); mask[Sset] = True
    phiS = phi * mask; phiQ = phi - phiS
    wS = float(phiS @ phiS) / n2
    nphiS = float(np.linalg.norm(phiS)); nphiQ = float(np.linalg.norm(phiQ))
    q = Um.T @ phiQ
    # --- Rayleigh-Ritz on span(S)  (this is EXACTLY what the repo does: sampled_akw.py:42-44)
    HS = np.asarray(H[Sset][:, Sset].todense())
    Es, Us = np.linalg.eigh(HS); del HS
    m = len(Es); epst = Es - E0
    d = Us.T @ phi[Sset]
    # --- build G = Q H |nt>  column-blocked; accumulate O-sums and gam-sums per eta
    Gfull = np.empty((D, m))
    nE = len(etas)
    C_cross = np.zeros(nE); C_ab = np.zeros(nE)
    alpha = [np.zeros(m, complex) for _ in range(nE)]
    beta = [np.zeros(D, complex) for _ in range(nE)]
    Ublk = np.zeros((D, nb))
    for j0 in range(0, m, nb):
        j1 = min(j0 + nb, m); k = j1 - j0
        Ublk[:, :k] = 0.0; Ublk[Sset, :k] = Us[:, j0:j1]
        Oblk = Um.T @ Ublk[:, :k]                       # <m|nt>
        Gblk = H @ Ublk[:, :k]; Gblk[Sset, :] = 0.0     # Q H |nt>
        Gfull[:, j0:j1] = Gblk
        gamblk = Um.T @ Gblk                            # <m|QH|nt>
        dd = d[j0:j1]
        Mblk = (c[:, None] * Oblk) * dd[None, :]        # c_m d_n O_mn
        Nblk = (q[:, None] * gamblk) * dd[None, :]      # q_m d_n gam_mn
        Dlt = eps[:, None] - epst[None, j0:j1]          # Delta_mn
        D2 = Dlt * Dlt
        for ie, eta in enumerate(etas):
            den = D2 + 4.0 * eta * eta
            C_cross[ie] += float((Mblk * (4.0 * eta * eta / den)).sum())
            C_ab[ie] += float((Nblk * (Dlt / den)).sum())
            W = Mblk / ((-Dlt) - 2j * eta)              # (epst_n - eps_m) - 2 i eta
            alpha[ie][j0:j1] += W.sum(axis=0)
            beta[ie] += W.sum(axis=1)
        del Oblk, Gblk, gamblk, Mblk, Nblk, Dlt, D2
    del Ublk, Us; gc.collect()
    # --- Lambda_S^2 = Re sum_{n,n'} d_n d_n' <g_n',g_n> * 2 i eta/((Et_n'-Et_n)+2 i eta)
    lam2 = np.zeros(nE)
    for i0 in range(0, m, nb):
        i1 = min(i0 + nb, m)
        Gi = Gfull[:, i0:i1]
        for j0 in range(0, m, nb):
            j1 = min(j0 + nb, m)
            Gram = Gi.T @ Gfull[:, j0:j1]               # <g_{n'}, g_n>, n' in I, n in J
            dE = epst[i0:i1, None] - epst[None, j0:j1]  # Et_n' - Et_n
            P = (d[i0:i1, None] * d[None, j0:j1]) * Gram
            for ie, eta in enumerate(etas):
                lam2[ie] += float((P * (2j * eta / (dE + 2j * eta))).real.sum())
            del Gram, dE, P
    del Gfull; gc.collect()
    Lam = np.sqrt(np.maximum(lam2, 0.0))
    # --- omega-resolved rungs
    allp = np.concatenate([eps, epst])
    rows = []
    for ie, eta in enumerate(etas):
        w, dw, ncore = make_grid(allp, eta, per_eta=per_eta, ntail=ntail)
        A = lorentz_sum(eps, c * c, w, eta); AS = lorentz_sum(epst, d * d, w, eta)
        u2 = A * (np.pi / eta); uS2 = AS * (np.pi / eta)
        # Re<u,u_S>
        al = alpha[ie]; be = beta[ie]
        re = np.zeros_like(w)
        for i in range(0, m, 400):
            x = w[None, :] - epst[i:i+400, None]; den = x * x + eta * eta
            re += ((al.real[i:i+400, None] * x + al.imag[i:i+400, None] * eta) / den).sum(0)
        for i in range(0, D, 400):
            x = w[None, :] - eps[i:i+400, None]; den = x * x + eta * eta
            re -= ((be.real[i:i+400, None] * x - be.imag[i:i+400, None] * eta) / den).sum(0)
        G2 = np.maximum(u2 + uS2 - 2.0 * re, 0.0)
        Gw = np.sqrt(G2); Fw = np.sqrt(u2) + np.sqrt(uS2)
        pre = eta / np.pi
        R0 = float(np.dot(np.abs(A - AS), dw))
        R1 = float(pre * np.dot(Fw * Gw, dw))
        NF = np.sqrt(pre * float(np.dot(Fw * Fw, dw)))
        NG_q = np.sqrt(pre * float(np.dot(G2, dw)))
        NG = np.sqrt(max(n2 + nphiS**2 - 2.0 * C_cross[ie], 0.0))     # closed form
        Nb2 = NG**2 - nphiQ**2 - 2.0 * C_ab[ie]
        Nb = np.sqrt(max(Nb2, 0.0))
        R2 = NF * NG; R3 = (nphi + nphiS) * NG
        R4 = (nphi + nphiS) * (nphiQ + Nb)
        R5 = (nphi + nphiS) * (nphiQ + Lam[ie] / eta)
        triv = n2 * (1.0 + wS)
        # quadrature self-checks
        chk_u = pre * float(np.dot(u2, dw)) / n2            # must be 1
        chk_uS = pre * float(np.dot(uS2, dw)) / max(nphiS**2, 1e-300)
        chk_NG = NG_q / max(NG, 1e-300)
        rows.append(dict(eta=float(eta), R0=R0, R1=R1, R2=R2, R3=R3, R4=R4, R5=R5,
                         trivial=float(triv), cert=float(min(R5, triv)),
                         NF=float(NF), NG=float(NG), NG_quad=float(NG_q), Nb=float(Nb),
                         Lam=float(Lam[ie]), Lam_over_eta=float(Lam[ie] / eta),
                         Lamhat_over_eta=float(Lam[ie] / nphi / eta),
                         C_cross=float(C_cross[ie]), C_ab=float(C_ab[ie]),
                         cos_FG=float(R1 / R2) if R2 > 0 else float('nan'),
                         chk_normu=float(chk_u), chk_normuS=float(chk_uS), chk_NG=float(chk_NG),
                         ngrid=int(len(w))))
        del w, dw, A, AS, u2, uS2, re, G2, Gw, Fw
        gc.collect()
    return dict(tag=tag, D=int(D), S=int(len(Sset)), FR=float(len(Sset) / D), w_S=float(wS),
                norm_phi2=float(n2), norm_phiS=float(nphiS), norm_phiQ=float(nphiQ),
                m_ritz=int(m), rows=rows)


def curves_for(Em, Um, phi, E0, H, Sset, eta, nb=128, per_eta=16.0, ntail=1200):
    """Return the omega-resolved curves for ONE eta:  w, dw, A, A_S, F=|u|+|u_S|, G=|u-u_S|."""
    D = len(phi); eps = Em - E0; c = Um.T @ phi
    mask = np.zeros(D, bool); mask[Sset] = True
    HS = np.asarray(H[Sset][:, Sset].todense()); Es, Us = np.linalg.eigh(HS); del HS
    m = len(Es); epst = Es - E0; d = Us.T @ phi[Sset]
    alpha = np.zeros(m, complex); beta = np.zeros(D, complex)
    Ublk = np.zeros((D, nb))
    for j0 in range(0, m, nb):
        j1 = min(j0 + nb, m); k = j1 - j0
        Ublk[:, :k] = 0.0; Ublk[Sset, :k] = Us[:, j0:j1]
        Oblk = Um.T @ Ublk[:, :k]
        Mblk = (c[:, None] * Oblk) * d[None, j0:j1]
        Dlt = eps[:, None] - epst[None, j0:j1]
        W = Mblk / ((-Dlt) - 2j * eta)
        alpha[j0:j1] += W.sum(axis=0); beta += W.sum(axis=1)
        del Oblk, Mblk, Dlt, W
    del Ublk, Us; gc.collect()
    w, dw, _ = make_grid(np.concatenate([eps, epst]), eta, per_eta=per_eta, ntail=ntail)
    A = lorentz_sum(eps, c * c, w, eta); AS = lorentz_sum(epst, d * d, w, eta)
    u2 = A * (np.pi / eta); uS2 = AS * (np.pi / eta)
    re = np.zeros_like(w)
    for i in range(0, m, 400):
        x = w[None, :] - epst[i:i+400, None]; den = x * x + eta * eta
        re += ((alpha.real[i:i+400, None] * x + alpha.imag[i:i+400, None] * eta) / den).sum(0)
    for i in range(0, D, 400):
        x = w[None, :] - eps[i:i+400, None]; den = x * x + eta * eta
        re -= ((beta.real[i:i+400, None] * x - beta.imag[i:i+400, None] * eta) / den).sum(0)
    G = np.sqrt(np.maximum(u2 + uS2 - 2.0 * re, 0.0))
    Fw = np.sqrt(u2) + np.sqrt(uS2)
    return w, dw, A, AS, Fw, G, eps, epst
