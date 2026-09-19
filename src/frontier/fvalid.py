# -*- coding: utf-8 -*-
r"""C3 validation gate.  Runs BEFORE any sweep point is believed.

(A) ALGEBRAIC CONTROLS at L=6 (dense reachable):
      A1  Gram_Ko(lean) == Gram_Ko(direct (HV)_out^T (HV)_out)     <- the memory trick
      A2  Gram_Ki(lean, rank one) == Gram_K - Gram_Ko(direct)
      A3  w_K == w_S exactly, P_K phi == phi_S, K subset span(S)
      A4  Lambda_K^2 == Lambda_in^2 + Lambda_out^2
      A5  at n_l = |S| (full):  Lambda_K == Lambda_S (dense Rayleigh-Ritz)
      A6  rel-L1 from Lorentz sums == rel-L1 from a dense eigendecomposition

(B) MANDATORY CROSS-VALIDATION at L=8, FR=0.85, eta=0.18 against the dense form:
      target of the gate verdict (protocol K=16):  Lambda_S = 0.0733666, Lambda_K = 0.0734851
      then REDONE with the repo protocol K=18 and the shift reported.
"""
import os, sys, time, gc, json, numpy as np
import scipy.sparse as sp
HERE = os.path.dirname(os.path.abspath(__file__))
# --- deposit paths (repaired 2026-09-19; see docs/C3_FRONTIER.md) ----------
SRC = os.path.dirname(HERE)                       # release/src
ROOT = os.path.dirname(SRC)                       # repository root
OUTDIR = os.environ.get('C3_OUT') or os.path.join(ROOT, 'data', 'c3_frontier')
TMP = os.environ.get('C3_TMP') or os.path.join(ROOT, 'build')
# --------------------------------------------------------------------------
sys.dont_write_bytecode = True
sys.path.insert(0, HERE)
import fcore as F
from fcore import log

ETAS = [0.10, 0.15, 0.18, 0.25]
OUT = {}


def dense_lambda_S(L, mask, phi, etas, U=8.0):
    """Dense Rayleigh-Ritz on span(S) -> Lambda_S, exactly as c0_sweep.prep does."""
    import c0_lib as C
    s = C.system(L, U)
    H = s['H']
    Sset = np.where(mask)[0]
    HS = np.asarray(H[Sset][:, Sset].todense())
    Es, Us = np.linalg.eigh(HS)
    cS = Us.T @ phi[Sset]
    H2 = (H @ H)[Sset][:, Sset].toarray()
    M = H2 - HS @ HS
    GramS = Us.T @ (M @ Us); GramS = 0.5 * (GramS + GramS.T)
    del HS, H2, M; gc.collect()
    from c0_sweep import lambda_from_gram
    res = {}
    for eta in etas:
        lam, _ = lambda_from_gram(Es, cS.astype(complex), GramS, eta)
        res[eta] = float(lam)
    return res, Es, cS ** 2, s


# ===========================================================================
def partA():
    L = 6; FR = 0.85
    log("=== (A) algebraic controls, L=6, FR=0.85, protocol K=18 ===")
    sysd = F.build_system(L)
    phi, mv, D = sysd['phi'], sysd['mv'], sysd['D']
    n2 = float(phi @ phi)
    order, wc, trank = F.ranking(sysd)
    mask, sdiag = F.subspace(order, wc, D, FR)
    k = sdiag['k']
    phiS = phi * mask; nphiS2 = float(phiS @ phiS); wS = nphiS2 / n2
    log(f"D={D} |S|={k} w_S={wS:.15f}  zero-wc in S={sdiag['n_zero_wc_in_S']} "
        f"ties_at_cut(in/total)={sdiag['n_tied_at_cut_inside']}/{sdiag['n_tied_at_cut_total']}")

    nl = k                                     # FULL depth -> P_K = P_S exactly
    V, al, be, nrm, _ = F.lanczos_S(mv, mask, phi, nl)
    m = V.shape[1]
    orth = float(np.abs(V.T @ V - np.eye(m)).max())
    resPKPS = float(np.abs(V[~mask, :]).max())
    PKphi = V @ (V.T @ phi)
    resPKphi = float(np.abs(PKphi - phiS).max())
    wK = float(PKphi @ PKphi) / n2
    log(f"A3: m={m} orth_defect={orth:.2e} res(K subset S)={resPKPS:.2e} "
        f"res(P_K phi - phi_S)={resPKphi:.2e}  |w_K-w_S|={abs(wK-wS):.2e}")

    # direct (expensive) reference: build HV explicitly, split by rows
    HV = np.column_stack([mv(V[:, j]) for j in range(m)])
    VHV = V.T @ HV; VHV = 0.5 * (VHV + VHV.T)
    G1_direct = HV.T @ HV; G1_direct = 0.5 * (G1_direct + G1_direct.T)
    HVo = HV[~mask, :]
    G2_direct = HVo.T @ HVo; G2_direct = 0.5 * (G2_direct + G2_direct.T)
    del HVo
    # lean
    G1_lean = F.gram_G1(mv, V)
    T = np.diag(al[:m]) + np.diag(be[:m - 1], 1) + np.diag(be[:m - 1], -1)
    e_G1 = float(np.abs(G1_lean - G1_direct).max() / max(np.abs(G1_direct).max(), 1e-300))
    e_T = float(np.abs(VHV - T).max())
    G2_lean = G1_lean - T @ T - (be[m - 1] ** 2) * np.outer(np.eye(m)[m - 1], np.eye(m)[m - 1])
    # compare in the COORDINATE basis (Sm cancels out of a similarity check)
    e_G2 = float(np.abs(G2_lean - (G2_direct - 0.0)).max() / max(np.abs(G2_direct).max(), 1e-300))
    log(f"A1: max|G1_lean-G1_direct|/max|G1| = {e_G1:.3e}")
    log(f"A1b: max|V^T H V - T(alpha,beta)| = {e_T:.3e}")
    log(f"A2: max|Gram_out(lean) - Gram_out(direct)|/max = {e_G2:.3e}")
    del HV, G1_direct, G2_direct, G2_lean, VHV; gc.collect()

    # A4 + A5
    dv = F.depth_eval(al, be, G1_lean, m, nphiS2, ETAS, keep_pairs=True)
    lamS, Es, wsD, s = dense_lambda_S(L, mask, phi, ETAS)
    a4 = []; a5 = []
    for eta in ETAS:
        lk, li, lo = dv['per_eta'][float(eta)]
        a4.append(abs(np.hypot(li, lo) - lk) / max(lk, 1e-300))
        a5.append(abs(lk - lamS[eta]) / max(lamS[eta], 1e-300))
        log(f"A4/A5 eta={eta:.2f}: Lam_K={lk:.10e}  sqrt(in^2+out^2)={np.hypot(li,lo):.10e}"
            f"  Lam_S(dense)={lamS[eta]:.10e}  rel.diff={a5[-1]:.2e}")

    # A6: rel-L1 from Lorentz sums vs dense eigendecomposition
    H = s['H']; E0 = s['E0']
    Em, Um = np.linalg.eigh(H.toarray()); cE = Um.T @ phi
    ref_th, ref_w, nref = F.ref_spectrum(mv, phi, min(400, D), D)
    a6 = []
    for eta in ETAS:
        d_dense = F.l1_lorentz(Em, cE ** 2, Es, wsD, eta)
        d_lanc = F.l1_lorentz(ref_th, ref_w, dv['theta'], dv['weights'], eta)
        a6.append(abs(d_dense['l1'] - d_lanc['l1']) / max(d_dense['l1'], 1e-300))
        log(f"A6 eta={eta:.2f}: relL1(dense)={d_dense['l1']/n2:.8e}  "
            f"relL1(lanczos-lean)={d_lanc['l1']/n2:.8e}  rel.diff={a6[-1]:.2e}  "
            f"tailfrac={d_dense['tail']/max(d_dense['l1'],1e-300):.1e}")
    OUT['A'] = dict(L=L, FR=FR, D=D, S=k, w_S=wS, orth_defect=orth,
                    res_K_in_S=resPKPS, res_PKphi=resPKphi, dw=abs(wK - wS),
                    e_G1=e_G1, e_T=e_T, e_G2=e_G2,
                    max_A4=float(max(a4)), max_A5=float(max(a5)), max_A6=float(max(a6)),
                    subspace_diag=sdiag)
    ok = (e_G1 < 1e-10 and e_G2 < 1e-10 and max(a4) < 1e-10 and max(a5) < 1e-8
          and max(a6) < 1e-4 and orth < 1e-10 and abs(wK - wS) < 1e-12)
    log(f"(A) VERDICT: {'PASS' if ok else 'FAIL'}")
    OUT['A']['pass'] = bool(ok)
    del V, G1_lean; gc.collect()
    return ok


# ===========================================================================
def partB(Kproto, nl=150):
    L = 8; FR = 0.85; eta_t = 0.18
    log(f"=== (B) L=8 cross-validation, FR=0.85, protocol K={Kproto}, n_l={nl} ===")
    sysd = F.build_system(L)
    phi, mv, D = sysd['phi'], sysd['mv'], sysd['D']
    n2 = float(phi @ phi)
    order, wc, trank = F.ranking(sysd, K=Kproto)
    mask, sdiag = F.subspace(order, wc, D, FR)
    phiS = phi * mask; nphiS2 = float(phiS @ phiS); wS = nphiS2 / n2
    log(f"|S|={sdiag['k']}/{D} w_S={wS:.12f}  zero-wc in S={sdiag['n_zero_wc_in_S']} "
        f"ties_at_cut(in/tot)={sdiag['n_tied_at_cut_inside']}/{sdiag['n_tied_at_cut_total']}")
    t = time.time()
    V, al, be, nrm, _ = F.lanczos_S(mv, mask, phi, nl)
    m = V.shape[1]; tl = time.time() - t
    orth = float(np.abs(V.T @ V - np.eye(m)).max())
    resPKPS = float(np.abs(V[~mask, :]).max())
    t = time.time(); G1 = F.gram_G1(mv, V); tg = time.time() - t
    del V; gc.collect()
    dv = F.depth_eval(al, be, G1, m, nphiS2, ETAS)
    t = time.time(); lamS, Es, wsD, s = dense_lambda_S(L, mask, phi, ETAS); tdense = time.time() - t
    rows = []
    for eta in ETAS:
        lk, li, lo = dv['per_eta'][float(eta)]
        rows.append(dict(eta=eta, Lambda_K=lk, Lambda_K_in=li, Lambda_K_out=lo,
                         Lambda_S_dense=lamS[eta], ratio_K_S=lk / lamS[eta]))
        log(f"  eta={eta:.2f}  Lambda_K={lk:.7e}  Lambda_S(dense)={lamS[eta]:.7e}  "
            f"K/S={lk/lamS[eta]:.6f}   (in={li:.2e} out={lo:.2e})")
    log(f"  cost: lanczos {tl:.1f}s  G1 {tg:.1f}s  dense Lambda_S {tdense:.1f}s  orth={orth:.1e}")
    res = dict(K=Kproto, nl=nl, S=sdiag['k'], D=D, w_S=wS, orth_defect=orth,
               res_K_in_S=resPKPS, rows=rows, subspace_diag=sdiag,
               t_lanczos=tl, t_G1=tg, t_dense=tdense)
    r18 = [r for r in rows if abs(r['eta'] - eta_t) < 1e-12][0]
    res['at_eta_0.18'] = dict(Lambda_S=r18['Lambda_S_dense'], Lambda_K=r18['Lambda_K'])
    return res


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if which in ('all', 'A'):
        partA()
    if which in ('all', 'B'):
        OUT['B_K16'] = partB(16, nl=150)
        log("TARGET (gate verdict, K=16): Lambda_S = 0.0733666 , Lambda_K = 0.0734851")
        g = OUT['B_K16']['at_eta_0.18']
        log("GOT    (K=16)              : Lambda_S = %.7f , Lambda_K = %.7f  -> dev %.2e / %.2e"
            % (g['Lambda_S'], g['Lambda_K'],
               abs(g['Lambda_S'] - 0.0733666) / 0.0733666,
               abs(g['Lambda_K'] - 0.0734851) / 0.0734851))
        OUT['B_K18'] = partB(18, nl=150)
        h = OUT['B_K18']['at_eta_0.18']
        log("GOT    (K=18, REPO proto)  : Lambda_S = %.7f , Lambda_K = %.7f  "
            "-> shift vs K=16: %.3f%% / %.3f%%"
            % (h['Lambda_S'], h['Lambda_K'],
               100 * (h['Lambda_S'] - g['Lambda_S']) / g['Lambda_S'],
               100 * (h['Lambda_K'] - g['Lambda_K']) / g['Lambda_K']))
        log("K=18 dense-vs-Krylov agreement: K/S = %.6f" % (h['Lambda_K'] / h['Lambda_S']))
    json.dump(OUT, open(os.path.join(OUTDIR, 'validation.json'), 'w'), indent=1, default=float)
    log("WROTE out/validation.json")
