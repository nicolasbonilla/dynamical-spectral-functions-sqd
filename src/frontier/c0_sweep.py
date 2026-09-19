# -*- coding: utf-8 -*-
"""C0 steps 3,4,5: Lambda_S vs Lambda_K over (FR, eta, n_l), measured cost, eta-freeness.

Everything eta-independent is done ONCE in prep();  eta only re-weights a cached Gram.
Dense Lambda_S uses  Gram_S = U^T [ (H^2)[S,S] - H_S^2 ] U   (no D x |S| matrix).
"""
import sys, os, json, time, numpy as np, scipy.sparse as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import c0_lib as C
from c0_lib import log

def lambda_from_gram(theta, c, Gram, eta, block=1024):
    nm = len(theta); tot = 0j; cc = np.conj(c)
    for i0 in range(0, nm, block):
        i1 = min(i0+block, nm)
        d = theta[i0:i1, None] - theta[None, :]
        F = (4*eta**2 + 2j*eta*d)/(d**2 + 4*eta**2)
        tot += np.sum(cc[i0:i1, None]*c[None, :]*Gram[i0:i1, :]*F)
    return np.sqrt(max(tot.real, 0.0)), float(tot.imag)

def prep(L, FR, nl, dense=True, K=16, dt=0.5):
    R = dict(L=L, FR=FR, nl=nl)
    s = C.system(L); H = s['H']; phi = s['phi']; E0 = s['E0']; D = s['D']
    n2 = float(phi@phi); R.update(D=D, E0=E0, norm_phi2=n2)
    t = time.time(); order, _ = C.born_order(H, phi, K=K, dt=dt); R['t_rank'] = time.time()-t
    k = max(1, int(round(FR*D))); Sset = np.sort(order[:k])
    mask = np.zeros(D, bool); mask[Sset] = True
    phiS = phi*mask; wS = float(phiS@phiS)/n2
    R.update(S=int(k), w_S=wS, norm_phiQ=float(np.linalg.norm(phi-phiS)))
    # exact poles/weights (needed for rel-L1 at any eta)
    Em, Um = np.linalg.eigh(H.toarray()); cE = Um.T@phi
    R['_ex'] = (Em-E0, cE**2)
    # ---------- dense Rayleigh-Ritz on span(S) ----------
    if dense:
        t = time.time()
        HS = np.asarray(H[Sset][:, Sset].todense())
        Es, Us = np.linalg.eigh(HS); cS = Us.T@phi[Sset]
        H2 = (H@H)[Sset][:, Sset].toarray()
        M = H2 - HS@HS
        Gram_S = Us.T@(M@Us); Gram_S = 0.5*(Gram_S+Gram_S.T)
        R['_S'] = (Es-E0, cS**2, Es, cS.astype(complex), Gram_S)
        R['t_prep_S'] = time.time()-t
        del HS, H2, M
    # ---------- Krylov on H_S, seed phi_S ----------
    def mvS(x): return mask*(H@(mask*x))
    t = time.time()
    V, al, be, nrm = C.lanczos(mvS, phiS.astype(float), nl, reorth=True)
    m = V.shape[1]; R['m'] = m; R['t_lanczos'] = time.time()-t
    t = time.time()
    HV = H@V                                   # D x m   (reuse: leakage needs full H)
    VHV = V.T@HV; VHV = 0.5*(VHV+VHV.T)
    th, Sm = np.linalg.eigh(VHV)
    cK = Sm.T@(V.T@phi)
    W = HV@Sm                                  # H |ritz_j>
    Gk = W - V@(V.T@W)                         # Q_K H |ritz_j>
    Gram_K = Gk.T@Gk; Gram_K = 0.5*(Gram_K+Gram_K.T)
    # split: out of span(S)  vs  in span(S) (= Krylov overflow)
    Gk_out = Gk.copy(); Gk_out[mask, :] = 0.0
    Gram_Kout = Gk_out.T@Gk_out
    Gram_Kin = Gram_K - Gram_Kout
    R['t_prep_K'] = time.time()-t
    R['_K'] = (th-E0, cK**2, th, cK.astype(complex), Gram_K, Gram_Kout, Gram_Kin)
    R['res_PK_PS'] = float(np.abs(V[~mask, :]).max())
    R['res_PKphi'] = float(np.abs(V@(V.T@phi)-phiS).max())
    R['w_K'] = float((V@(V.T@phi))@(V@(V.T@phi)))/n2
    R['orth_defect'] = float(np.abs(V.T@V-np.eye(m)).max())
    R['beta_last'] = float(be[m-1])
    VHSV = V.T@np.column_stack([mvS(V[:, j]) for j in range(m)])
    R['res_PKHPK'] = float(np.abs(VHV-0.5*(VHSV+VHSV.T)).max())
    del V, HV, W, Gk, Gk_out
    return R

def evaluate(R, etas):
    n2 = R['norm_phi2']; nphi = np.sqrt(n2); wS = R['w_S']
    pe, we = R['_ex']
    lo = min(pe.min(), 0)-30; hi = max(pe.max(), 0)+30
    grid = np.linspace(lo, hi, 30000)
    rows = []
    for eta in etas:
        A_ex = np.zeros_like(grid)
        for a, b in zip(pe, we): A_ex += b*(eta/np.pi)/((grid-a)**2+eta**2)
        nrmA = float(np.trapz(np.abs(A_ex), grid))
        row = dict(eta=eta, int_A=float(np.trapz(A_ex, grid)))
        if '_S' in R:
            ps, ws, Es, cS, GrS = R['_S']
            A_S = np.zeros_like(grid)
            for a, b in zip(ps, ws): A_S += b*(eta/np.pi)/((grid-a)**2+eta**2)
            row['relL1_S'] = float(np.trapz(np.abs(A_S-A_ex), grid)/nrmA)
            t = time.time(); lamS, im = lambda_from_gram(Es, cS, GrS, eta); row['t_eval_S'] = time.time()-t
            row['Lambda_S'] = float(lamS); row['LamHat_S_over_eta'] = float(lamS/nphi/eta)
        pk, wk, th, cK, GrK, GrKo, GrKi = R['_K']
        A_K = np.zeros_like(grid)
        for a, b in zip(pk, wk): A_K += b*(eta/np.pi)/((grid-a)**2+eta**2)
        row['relL1_K'] = float(np.trapz(np.abs(A_K-A_ex), grid)/nrmA)
        t = time.time()
        lamK, _ = lambda_from_gram(th, cK, GrK, eta)
        row['t_eval_K'] = time.time()-t
        lamKo, _ = lambda_from_gram(th, cK, GrKo, eta)
        lamKi, _ = lambda_from_gram(th, cK, GrKi, eta)
        row['Lambda_K'] = float(lamK); row['Lambda_K_outS'] = float(lamKo); row['Lambda_K_inS'] = float(lamKi)
        row['LamHat_K_over_eta'] = float(lamK/nphi/eta)
        if 'Lambda_S' in row: row['ratio_K_S'] = float(lamK/max(row['Lambda_S'], 1e-300))
        # Theorem B, relative form
        sw = np.sqrt(wS)
        row['boundB_S'] = (1+sw)*(np.sqrt(max(1-wS, 0))+row.get('LamHat_S_over_eta', np.nan))
        row['boundB_K'] = (1+sw)*(np.sqrt(max(1-wS, 0))+row['LamHat_K_over_eta'])
        row['bound_trivial'] = 1+wS
        rows.append(row)
    return rows

if __name__ == '__main__':
    L = int(sys.argv[1]); FR = float(sys.argv[2]); nl = int(sys.argv[3])
    dense = (len(sys.argv) < 5 or sys.argv[4] != 'nodense')
    etas = [0.05, 0.10, 0.15, 0.18, 0.25, 0.50]
    R = prep(L, FR, nl, dense=dense)
    log(f"L={L} FR={FR} nl={nl}: |S|={R['S']}/{R['D']} w_S={R['w_S']:.12f} "
        f"orth={R['orth_defect']:.1e} res(PkPs)={R['res_PK_PS']:.1e} "
        f"res(PkHPk)={R['res_PKHPK']:.1e} res(Pkphi)={R['res_PKphi']:.1e} "
        f"|wK-wS|={abs(R['w_K']-R['w_S']):.1e}")
    log(f"   t_rank={R['t_rank']:.1f}s t_prep_S={R.get('t_prep_S',float('nan')):.1f}s "
        f"t_lanczos={R['t_lanczos']:.1f}s t_prep_K={R['t_prep_K']:.1f}s")
    rows = evaluate(R, etas)
    for r in rows:
        log("   eta=%.2f  relL1_S=%.4e relL1_K=%.4e | Lam_S=%.4e Lam_K=%.4e  K/S=%.4f | "
            "Lhat_S/eta=%.4f Lhat_K/eta=%.4f | boundK=%.4f trivial=%.4f | t_evalK=%.3fs" % (
                r['eta'], r.get('relL1_S', float('nan')), r['relL1_K'], r.get('Lambda_S', float('nan')),
                r['Lambda_K'], r.get('ratio_K_S', float('nan')), r.get('LamHat_S_over_eta', float('nan')),
                r['LamHat_K_over_eta'], r['boundB_K'], r['bound_trivial'], r['t_eval_K']))
    meta = {k: v for k, v in R.items() if not k.startswith('_')}
    json.dump(dict(meta=meta, rows=rows),
              open(os.path.join(HERE, f"sweep_L{L}_FR{FR:.2f}_nl{nl}.json"), 'w'), indent=1)
