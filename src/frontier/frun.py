# -*- coding: utf-8 -*-
r"""C3 sweep driver.  ONE L per process (so the OS reclaims memory between L).

  python frun.py <L> <nl> <FR,FR,...> [nref] [tag]

Per (L,FR) it writes out/pt_L{L}_FR{FR}_nl{nl}.json containing, for every eta and
every Krylov depth on the ladder:
    w_S, Lambda_K, Lambda_K_in, Lambda_K_out, leakage bound, trivial bound,
    effective certificate min{leak,trivial}, the measured rel-L1, the slack,
    and the convergence flag derived from Lambda_K_in.

Honesty rules enforced here:
  * full reorthogonalisation always (gate condition 1)
  * repo protocol K=18, dt=0.5, sub=16 (gate condition 2)
  * the certificate is reported as certifying A_K and ONLY A_K (gate condition 4)
  * tie / zero-Born-weight diagnostics of the ranking are recorded (gate condition 5)
  * rel-L1 is ||A - A_K||_1 / ||phi||^2 : the THEOREM's normalisation, over the whole
    real line (analytic tails), not the repo's window-normalised variant -- both are
    stored so the two objects are never silently confused.
"""
import os, sys, time, gc, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
# --- deposit paths (repaired 2026-09-19; see docs/C3_FRONTIER.md) ----------
SRC = os.path.dirname(HERE)                       # release/src
ROOT = os.path.dirname(SRC)                       # repository root
OUT = os.environ.get('C3_OUT') or os.path.join(ROOT, 'data', 'c3_frontier')
TMP = os.environ.get('C3_TMP') or os.path.join(ROOT, 'build')
# --------------------------------------------------------------------------
sys.dont_write_bytecode = True
sys.path.insert(0, HERE)
import fcore as F
from fcore import log
import c0_lib as C

ETAS = [0.10, 0.15, 0.18, 0.25]


def bare_ref(mv, phi, nl_ref):
    """3-term Lanczos WITHOUT reorthogonalisation: memory = 3 vectors.
    Ghosts duplicate eigenvalues but the Gauss quadrature for A(w) is unaffected;
    the difference against a reorthogonalised reference is measured, never assumed."""
    n = phi.shape[0]
    nrm = float(np.linalg.norm(phi))
    v = phi / nrm
    al = np.zeros(nl_ref); be = np.zeros(nl_ref)
    vm = np.zeros_like(v)
    w = mv(v); a = float(v @ w); al[0] = a; w = w - a * v
    m = nl_ref
    for j in range(1, nl_ref):
        b = float(np.linalg.norm(w)); be[j - 1] = b
        if b < 1e-12:
            m = j; break
        vm = v; v = w / b
        w = mv(v); a = float(v @ w); al[j] = a
        w = w - a * v - b * vm
    T = np.diag(al[:m]) + np.diag(be[:m - 1], 1) + np.diag(be[:m - 1], -1)
    th, Sm = np.linalg.eigh(T)
    return th, (nrm ** 2) * Sm[0, :] ** 2, m


def run(L, nl, FRs, nref=400, nref_reorth=0, tag=''):
    t_all = time.time()
    sysd = F.build_system(L)
    phi, mv, D, E0 = sysd['phi'], sysd['mv'], sysd['D'], sysd['E0']
    n2 = float(phi @ phi)
    log(f"L={L} D={D} E0={E0:.8f} |phi|^2={n2:.8f}  [gs {sysd['t_gs']:.1f}s]")
    order, wc, trank = F.ranking(sysd)
    log(f"ranking (repo protocol K=18,dt=0.5,sub=16, 19 snapshots) [{trank:.1f}s]; "
        f"exact-zero Born weights in sector: {int(np.count_nonzero(wc==0.0))}/{D}")

    # ---------------- exact reference A(w) ----------------
    t = time.time()
    th_r, w_r, m_r = bare_ref(mv, phi, min(nref, D))
    t_ref = time.time() - t
    th_r2, w_r2, _ = bare_ref(mv, phi, min(max(nref // 2, 20), D))
    ref_drift = {}
    for eta in ETAS:
        d = F.l1_lorentz(th_r, w_r, th_r2, w_r2, eta)
        ref_drift[eta] = d['l1'] / n2
    log(f"reference: bare Lanczos m={m_r} [{t_ref:.1f}s]; drift nref/2 -> nref (rel-L1/|phi|^2): "
        + " ".join("eta=%.2f:%.2e" % (e, ref_drift[e]) for e in ETAS))
    ref_reorth_gap = {}
    if nref_reorth:
        th_o, w_o, m_o = F.ref_spectrum(mv, phi, min(nref_reorth, D), D)
        for eta in ETAS:
            d = F.l1_lorentz(th_r, w_r, th_o, w_o, eta)
            ref_reorth_gap[eta] = d['l1'] / n2
        log(f"reference: bare vs REORTH(m={m_o}) gap: "
            + " ".join("eta=%.2f:%.2e" % (e, ref_reorth_gap[e]) for e in ETAS))
        del th_o, w_o
    gc.collect()

    meta = dict(L=L, D=D, E0=E0, norm_phi2=n2, nl_req=nl, nref=m_r,
                protocol=dict(F.PROTO), t_gs=sysd['t_gs'], t_rank=trank, t_ref=t_ref,
                ref_drift={str(k): v for k, v in ref_drift.items()},
                ref_reorth_gap={str(k): v for k, v in ref_reorth_gap.items()},
                n_zero_wc_sector=int(np.count_nonzero(wc == 0.0)))

    for FR in FRs:
        fn = os.path.join(OUT, f"pt_L{L}_FR{FR:.2f}_nl{nl}{tag}.json")
        if os.path.exists(fn):
            log(f"SKIP existing {os.path.basename(fn)}"); continue
        t0 = time.time()
        mask, sdiag = F.subspace(order, wc, D, FR)
        k = sdiag['k']
        phiS = phi * mask; nphiS2 = float(phiS @ phiS); wS = nphiS2 / n2
        nl_eff = int(min(nl, k))
        t = time.time()
        V, al, be, nrm, _ = F.lanczos_S(mv, mask, phi, nl_eff)
        m = V.shape[1]; t_lan = time.time() - t
        # cheap exact controls
        res_K_in_S = float(np.abs(V[~mask, :]).max()) if k < D else 0.0
        PKphi = V @ (V.T @ phi)
        res_PKphi = float(np.abs(PKphi - phiS).max()); wK = float(PKphi @ PKphi) / n2
        del PKphi
        idx = np.unique(np.linspace(0, m - 1, min(m, 48)).astype(int))
        Osub = V.T @ V[:, idx]; Osub[idx, np.arange(len(idx))] -= 1.0
        orth = float(np.abs(Osub).max()); del Osub
        t = time.time()
        G1 = F.gram_G1(mv, V, chunk_log=200 if D > 2e5 else None)
        t_g1 = time.time() - t
        del V; gc.collect()

        lad = F.ladder(m)
        depths = []
        for n in lad:
            depths.append(F.depth_eval(al, be, G1, n, nphiS2, ETAS,
                                       keep_pairs=(n == lad[-1] or n == lad[max(0, len(lad) - 3)])))
        # rel-L1 at the two deepest kept depths
        rel = {}
        for dv in depths:
            if 'theta' not in dv:
                continue
            rr = {}
            for eta in ETAS:
                d = F.l1_lorentz(th_r, w_r, dv['theta'], dv['weights'], eta)
                rr[str(eta)] = dict(relL1=d['l1'] / n2, relL1_win=d['l1'] / d['win_absA1'],
                                    tail_frac=d['tail'] / max(d['l1'], 1e-300),
                                    int_AK=float(dv['weights'].sum()))
            rel[str(dv['n'])] = rr
            dv.pop('theta'); dv.pop('weights')

        rows = []
        for dv in depths:
            for eta in ETAS:
                lk, li, lo = dv['per_eta'][float(eta)]
                sw = np.sqrt(wS); sq = np.sqrt(max(1 - wS, 0.0))
                lh = lk / np.sqrt(n2) / eta
                lho = lo / np.sqrt(n2) / eta
                B, triv, cert = F.certificate(wS, lh)
                Bo, _, certo = F.certificate(wS, lho)
                r = dict(n=dv['n'], eta=eta, beta_n=dv['beta_n'],
                         Lambda_K=lk, Lambda_K_in=li, Lambda_K_out=lo,
                         LamHat_over_eta=lh, LamHat_out_over_eta=lho,
                         bound_leak=B, bound_trivial=triv, cert=cert,
                         bound_leak_outonly=Bo, cert_outonly=certo,
                         in_over_out=li / max(lo, 1e-300),
                         in_frac_of_Lam2=(li ** 2) / max(lk ** 2, 1e-300))
                rr = rel.get(str(dv['n']), {}).get(str(eta))
                if rr:
                    r.update(relL1=rr['relL1'], relL1_win=rr['relL1_win'],
                             tail_frac=rr['tail_frac'], int_AK=rr['int_AK'],
                             slack=cert / max(rr['relL1'], 1e-300),
                             slack_leakbranch=B / max(rr['relL1'], 1e-300))
                rows.append(r)

        # Lambda_out plateau test (is the PHYSICAL part of the certificate converged?)
        plateau = {}
        for eta in ETAS:
            seq = [(d['n'], d['per_eta'][float(eta)][2]) for d in depths]
            if len(seq) >= 2:
                a_, b_ = seq[-2][1], seq[-1][1]
                plateau[str(eta)] = dict(n_prev=seq[-2][0], n_top=seq[-1][0],
                                         Lout_prev=a_, Lout_top=b_,
                                         rel_change=abs(b_ - a_) / max(b_, 1e-300),
                                         monotone_up=bool(all(seq[i][1] <= seq[i + 1][1] * (1 + 1e-9)
                                                              for i in range(len(seq) - 1))))
        pt = dict(meta=meta, FR=FR, S=k, FR_real=k / D, w_S=wS, w_K=wK, dw=abs(wK - wS),
                  nl=m, orth_defect_sampled=orth, res_K_in_S=res_K_in_S,
                  res_PKphi=res_PKphi, subspace_diag=sdiag,
                  t_lanczos=t_lan, t_G1=t_g1, t_point=time.time() - t0,
                  peak_vec_GiB=m * D * 8 / 2 ** 30,
                  ladder=lad, rows=rows, plateau=plateau)
        json.dump(pt, open(fn, 'w'), indent=1, default=float)
        top = [r for r in rows if r['n'] == m]
        log(f"L={L} FR={FR:.2f} |S|={k}({k/D:.3f}) w_S={wS:.9f} nl={m} "
            f"[lanczos {t_lan:.1f}s + G1 {t_g1:.1f}s = {time.time()-t0:.1f}s, "
            f"peak {m*D*8/2**30:.2f} GiB] ties@cut {sdiag['n_tied_at_cut_inside']}/"
            f"{sdiag['n_tied_at_cut_total']} zeroBorn {sdiag['n_zero_wc_in_S']}")
        for r in top:
            log("    eta=%.2f  Lam/eta=%.4f  B=%.4f  triv=%.4f  cert=%.4f  "
                "relL1=%.4e  slack=%.1f  in/out=%.3f  B_out=%.4f"
                % (r['eta'], r['LamHat_over_eta'], r['bound_leak'], r['bound_trivial'],
                   r['cert'], r.get('relL1', float('nan')), r.get('slack', float('nan')),
                   r['in_over_out'], r['bound_leak_outonly']))
        del G1, al, be; gc.collect()
    log(f"L={L} ALL DONE in {time.time()-t_all:.1f}s")


if __name__ == '__main__':
    L = int(sys.argv[1]); nl = int(sys.argv[2])
    FRs = [float(x) for x in sys.argv[3].split(',')]
    nref = int(sys.argv[4]) if len(sys.argv) > 4 else 400
    nrefo = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    tag = sys.argv[6] if len(sys.argv) > 6 else ''
    run(L, nl, FRs, nref, nrefo, tag)
