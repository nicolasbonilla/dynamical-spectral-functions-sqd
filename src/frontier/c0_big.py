# -*- coding: utf-8 -*-
"""C0 step 4+5 at L=10 and L=12: MEASURED cost of Lambda_K (Krylov route), matrix-free.
No dense object of size |S| is ever formed.  Peak memory is ~3 x (D x n_l) float64.
"""
import sys, os, gc, json, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import c0_lib as C
from c0_lib import log, AK
from scipy.sparse.linalg import eigsh, LinearOperator

def sector_mv(L, U, nup, ndn, t=1.0):
    Tu, Su, iu = AK.hop(L, nup, t); Td, Sd, idd = AK.hop(L, ndn, t)
    Du, Dd = len(Su), len(Sd)
    upocc = np.array([[(m >> i) & 1 for i in range(L)] for m in Su], dtype=np.float64)
    dnocc = np.array([[(m >> i) & 1 for i in range(L)] for m in Sd], dtype=np.float64)
    diagU = U*(upocc@dnocc.T).ravel()
    def mv(x):
        X = x.reshape(Du, Dd); Y = Tu@X + (Td@X.T).T
        return Y.ravel()+diagU*x
    return mv, Du*Dd, Du, Dd

def run(L, FR, nl, U=8.0, etas=(0.05, 0.10, 0.15, 0.18, 0.25, 0.50), Kt=16, dt=0.5):
    R = dict(L=L, FR=FR, nl=nl); T = {}
    nup = nd = L//2
    t = time.time()
    mv0, D0, Du0, Dd0 = sector_mv(L, U, nup, nd)
    H0 = LinearOperator((D0, D0), matvec=mv0, dtype=float)
    E0, V0 = eigsh(H0, k=1, which='SA', ncv=12, tol=1e-10); E0 = float(E0[0])
    Psi = V0[:, 0].reshape(Du0, Dd0); del V0
    cdU = AK.cdag_map(L, nup)
    phi = np.asarray((cdU[0]@Psi).reshape(-1), dtype=float); del Psi, cdU
    mv1, D, Du1, Dd1 = sector_mv(L, U, nup+1, nd)
    T['gs'] = time.time()-t
    n2 = float(phi@phi); R.update(D=D, E0=E0, norm_phi2=n2)
    log(f"L={L}: D={D} E0={E0:.6f} |phi|^2={n2:.6f}  [gs {T['gs']:.1f}s]")

    t = time.time()
    class _H:
        def __matmul__(self, x): return mv1(x)
    order, _ = C.born_order_mf(_H(), phi, K=Kt, dt=dt, sub=16, m=6)
    T['rank'] = time.time()-t
    k = max(1, int(round(FR*D))); Sset = np.sort(order[:k]); del order
    mask = np.zeros(D, bool); mask[Sset] = True; del Sset
    phiS = phi*mask; wS = float(phiS@phiS)/n2
    R.update(S=int(k), w_S=wS, norm_phiQ=float(np.linalg.norm(phi-phiS)))
    log(f"L={L}: |S|={k} ({k/D:.3f})  w_S={wS:.12f}  [rank {T['rank']:.1f}s]")
    gc.collect()

    # ---- Lanczos on H_S with full reorthogonalization ----
    def mvS(x): return mask*mv1(mask*x)
    t = time.time()
    V, al, be, nrm = C.lanczos(mvS, phiS.astype(float), nl, reorth=True)
    T['lanczos'] = time.time()-t; m = V.shape[1]
    R['m'] = m; R['orth_defect'] = float(np.abs(V.T@V-np.eye(m)).max())
    R['res_PK_PS'] = float(np.abs(V[~mask, :]).max())
    PKphi = V@(V.T@phi)
    R['res_PKphi'] = float(np.abs(PKphi-phiS).max())
    R['w_K'] = float(PKphi@PKphi)/n2; del PKphi
    log(f"L={L}: lanczos m={m} orth={R['orth_defect']:.2e} PkPs={R['res_PK_PS']:.2e} "
        f"Pkphi={R['res_PKphi']:.2e} |wK-wS|={abs(R['w_K']-wS):.2e}  [{T['lanczos']:.1f}s]")

    # ---- leakage Gram on the Krylov subspace ----
    t = time.time()
    HV = np.empty_like(V)
    for j in range(m): HV[:, j] = mv1(V[:, j])
    T['HV'] = time.time()-t
    t = time.time()
    VHV = V.T@HV; VHV = 0.5*(VHV+VHV.T)
    th, Sm = np.linalg.eigh(VHV)
    cK = Sm.T@(V.T@phi)
    W = HV@Sm; del HV; gc.collect()
    W -= V@(V.T@W)                              # Q_K H |ritz>
    Gram_K = W.T@W; Gram_K = 0.5*(Gram_K+Gram_K.T)
    Wo = W[~mask, :]
    Gram_Ko = Wo.T@Wo; del Wo
    Gram_Ki = Gram_K-Gram_Ko
    T['gram'] = time.time()-t
    R['beta_last'] = float(be[m-1])
    # diagnostic (NOT part of the timed Lambda_K cost): V^T H V  ==  V^T H_S V ?
    A1 = V.T@np.column_stack([mvS(V[:, j]) for j in range(m)]); A1 = 0.5*(A1+A1.T)
    R['res_PKHPK'] = float(np.abs(VHV-A1).max()); del A1
    log(f"L={L}: HV {T['HV']:.1f}s   gram {T['gram']:.1f}s   res(PkHPk)={R['res_PKHPK']:.2e}")

    del W, V; gc.collect()

    # ---- eta sweep (free: only the tiny nl x nl kernel changes) ----
    import c0_sweep as SW
    rows = []
    nphi = np.sqrt(n2)
    for eta in etas:
        t = time.time(); lamK, _ = SW.lambda_from_gram(th, cK.astype(complex), Gram_K, eta)
        te = time.time()-t
        lamKo, _ = SW.lambda_from_gram(th, cK.astype(complex), Gram_Ko, eta)
        lamKi, _ = SW.lambda_from_gram(th, cK.astype(complex), Gram_Ki, eta)
        sw = np.sqrt(wS)
        rows.append(dict(eta=eta, Lambda_K=float(lamK), Lambda_K_outS=float(lamKo),
                         Lambda_K_inS=float(lamKi), LamHat_K_over_eta=float(lamK/nphi/eta),
                         boundB_K=float((1+sw)*(np.sqrt(max(1-wS, 0))+lamK/nphi/eta)),
                         bound_trivial=float(1+wS), t_eval=te))
        log("   eta=%.2f  Lam_K=%.5e  Lh_K/eta=%.4f  bK=%.4f  triv=%.4f  (inS %.3e / outS %.3e)  t=%.4fs"
            % (eta, lamK, lamK/nphi/eta, rows[-1]['boundB_K'], 1+wS, lamKi, lamKo, te))
    R['rows'] = rows; R['timings'] = T
    R['t_Lambda_total'] = T['lanczos']+T['HV']+T['gram']
    log(f"L={L} TOTAL Lambda_K cost (lanczos+HV+gram) = {R['t_Lambda_total']:.1f}s ; "
        f"+rank {T['rank']:.1f}s +gs {T['gs']:.1f}s")
    return R

if __name__ == '__main__':
    L = int(sys.argv[1]); FR = float(sys.argv[2]); nl = int(sys.argv[3])
    R = run(L, FR, nl)
    json.dump({k: v for k, v in R.items()},
              open(os.path.join(HERE, f"big_L{L}_FR{FR:.2f}_nl{nl}.json"), 'w'), indent=1)
    log("WROTE")
