# -*- coding: utf-8 -*-
r"""C3 at L=12 (and any large L): the same measurement as frun.py, but with the
Krylov basis stored in the |S|-dimensional COORDINATES of S instead of the full
sector, and optionally backed by a memmap in $C3_TMP (default: build/).

  memory(V) = |S| * n_l * 8   instead of   D * n_l * 8      (factor FR)

Nothing about the physics or the protocol changes:
  * the Lanczos recursion below is the two-pass classical Gram-Schmidt of
    c0_lib.lanczos, algorithmically IDENTICAL -- only the storage layout differs
    (full reorthogonalisation is NOT relaxed; that is gate condition 1);
  * K = 18, dt = 0.5, sub = 16 ranking (gate condition 2);
  * K subset span(S) is exact by construction here (res = 0), so that control is
    replaced by the orthogonality defect and by P_K phi == phi_S.

  python frun12.py <L> <nl> <FR,...> [nref] [mmap:0|1] [tag]
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
from frun import bare_ref

ETAS = [0.10, 0.15, 0.18, 0.25]


def lanczos_coords(mvS_coord, v0, nl, V):
    """Two-pass classical Gram-Schmidt Lanczos, writing into the preallocated
    (possibly memmapped) V.  Same algorithm as c0_lib.lanczos(reorth=True)."""
    nrm = float(np.linalg.norm(v0))
    V[:, 0] = v0 / nrm
    al = np.zeros(nl); be = np.zeros(nl); m = nl
    w = mvS_coord(np.asarray(V[:, 0]))
    a = float(V[:, 0] @ w); al[0] = a
    w = w - a * np.asarray(V[:, 0])
    for j in range(1, nl):
        for _ in range(2):
            Vj = V[:, :j]
            w -= Vj @ (Vj.T @ w)
        b = float(np.linalg.norm(w)); be[j - 1] = b
        if b < 1e-12:
            m = j; break
        V[:, j] = w / b
        w = mvS_coord(np.asarray(V[:, j]))
        a = float(V[:, j] @ w); al[j] = a
        w = w - a * np.asarray(V[:, j]) - b * np.asarray(V[:, j - 1])
    if m == nl:
        for _ in range(2):
            Vj = V[:, :nl]
            w -= Vj @ (Vj.T @ w)
        be[nl - 1] = float(np.linalg.norm(w))
    return al[:m], be[:m], nrm, m


def run(L, nl, FRs, nref=800, use_mmap=1, tag=''):
    t_all = time.time()
    sysd = F.build_system(L)
    phi, mv, D, E0 = sysd['phi'], sysd['mv'], sysd['D'], sysd['E0']
    n2 = float(phi @ phi)
    log(f"L={L} D={D} E0={E0:.8f} |phi|^2={n2:.8f}  [gs {sysd['t_gs']:.1f}s]")
    order, wc, trank = F.ranking(sysd)
    log(f"ranking (repo protocol K=18,dt=0.5,sub=16) [{trank:.1f}s]; "
        f"exact-zero Born weights: {int(np.count_nonzero(wc==0.0))}/{D}")
    t = time.time()
    th_r, w_r, m_r = bare_ref(mv, phi, min(nref, D))
    th_r2, w_r2, _ = bare_ref(mv, phi, min(max(nref // 2, 20), D))
    ref_drift = {}
    for eta in ETAS:
        ref_drift[eta] = F.l1_lorentz(th_r, w_r, th_r2, w_r2, eta)['l1'] / n2
    del th_r2, w_r2; gc.collect()
    log(f"reference: bare Lanczos m={m_r} [{time.time()-t:.1f}s]; drift nref/2->nref: "
        + " ".join("eta=%.2f:%.2e" % (e, ref_drift[e]) for e in ETAS))

    meta = dict(L=L, D=D, E0=E0, norm_phi2=n2, nl_req=nl, nref=m_r,
                protocol=dict(F.PROTO), t_gs=sysd['t_gs'], t_rank=trank,
                ref_drift={str(k): v for k, v in ref_drift.items()},
                ref_reorth_gap={}, n_zero_wc_sector=int(np.count_nonzero(wc == 0.0)),
                storage='S-coordinates' + ('+memmap' if use_mmap else ''))

    for FR in FRs:
        fn = os.path.join(OUT, f"pt_L{L}_FR{FR:.2f}_nl{nl}{tag}.json")
        if os.path.exists(fn):
            log(f"SKIP existing {os.path.basename(fn)}"); continue
        t0 = time.time()
        mask, sdiag = F.subspace(order, wc, D, FR)
        k = sdiag['k']
        idx = np.where(mask)[0]
        phiS_c = phi[idx].copy()
        nphiS2 = float(phiS_c @ phiS_c); wS = nphiS2 / n2
        nl_eff = int(min(nl, k))
        gib = k * nl_eff * 8 / 2 ** 30
        log(f"L={L} FR={FR:.2f} |S|={k} ({k/D:.3f}) w_S={wS:.10f}  V = {gib:.2f} GiB "
            f"({'memmap' if use_mmap else 'RAM'})")

        buf = np.zeros(D)

        def mvS_coord(xs):
            buf[:] = 0.0; buf[idx] = xs
            return mv(buf)[idx]

        mmpath = os.path.join(TMP, f"_V_L{L}_FR{FR:.2f}.dat")
        if not os.path.isdir(TMP):
            os.makedirs(TMP)
        # COLUMN-MAJOR is essential: the Lanczos basis is accessed one column at a
        # time, and a C-ordered (|S| x n_l) block would stride by n_l*8 bytes per
        # element -- ruinous for a memmap (one page fault per element).
        if use_mmap:
            V = np.memmap(mmpath, dtype=np.float64, mode='w+', shape=(k, nl_eff), order='F')
        else:
            V = np.zeros((k, nl_eff), order='F')
        t = time.time()
        al, be, nrm, m = lanczos_coords(mvS_coord, phiS_c, nl_eff, V)
        t_lan = time.time() - t
        V = V[:, :m]
        # controls
        jsel = np.unique(np.linspace(0, m - 1, min(m, 40)).astype(int))
        Osub = np.asarray(V.T @ np.asarray(V[:, jsel]))
        Osub[jsel, np.arange(len(jsel))] -= 1.0
        orth = float(np.abs(Osub).max()); del Osub
        Vtphi = np.asarray(V.T @ phiS_c)
        res_PKphi = float(np.abs(np.asarray(V @ Vtphi) - phiS_c).max())
        wK = float(Vtphi @ Vtphi) / n2
        # G1 = V^T H^2 V, column by column, one temp vector
        t = time.time()
        G1 = np.empty((m, m))
        for j in range(m):
            buf[:] = 0.0; buf[idx] = np.asarray(V[:, j])
            y = mv(mv(buf))[idx]
            G1[:, j] = np.asarray(V.T @ y)
            if (j + 1) % 50 == 0:
                log(f"    G1 column {j+1}/{m}  [{time.time()-t:.0f}s]")
        G1 = 0.5 * (G1 + G1.T); t_g1 = time.time() - t
        del V; gc.collect()
        if use_mmap:
            try:
                os.remove(mmpath)
            except OSError:
                pass

        lad = F.ladder(m)
        depths = [F.depth_eval(al, be, G1, n, nphiS2, ETAS,
                               keep_pairs=(n == lad[-1] or n == lad[max(0, len(lad) - 3)]))
                  for n in lad]
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
                lh = lk / np.sqrt(n2) / eta; lho = lo / np.sqrt(n2) / eta
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
        plateau = {}
        for eta in ETAS:
            seq = [(d['n'], d['per_eta'][float(eta)][2]) for d in depths]
            a_, b_ = seq[-2][1], seq[-1][1]
            plateau[str(eta)] = dict(n_prev=seq[-2][0], n_top=seq[-1][0], Lout_prev=a_,
                                     Lout_top=b_, rel_change=abs(b_ - a_) / max(b_, 1e-300),
                                     monotone_up=False)
        pt = dict(meta=meta, FR=FR, S=k, FR_real=k / D, w_S=wS, w_K=wK, dw=abs(wK - wS),
                  nl=m, orth_defect_sampled=orth, res_K_in_S=0.0, res_PKphi=res_PKphi,
                  subspace_diag=sdiag, t_lanczos=t_lan, t_G1=t_g1,
                  t_point=time.time() - t0, peak_vec_GiB=gib, ladder=lad,
                  rows=rows, plateau=plateau)
        json.dump(pt, open(fn, 'w'), indent=1, default=float)
        log(f"  done [lanczos {t_lan:.0f}s + G1 {t_g1:.0f}s = {time.time()-t0:.0f}s] "
            f"orth={orth:.1e} |wK-wS|={abs(wK-wS):.1e} ties@cut "
            f"{sdiag['n_tied_at_cut_inside']}/{sdiag['n_tied_at_cut_total']} "
            f"zeroBorn {sdiag['n_zero_wc_in_S']}")
        for r in [x for x in rows if x['n'] == m]:
            log("    eta=%.2f  Lam/eta=%.4f  B=%.4f  triv=%.4f  cert=%.4f  relL1=%.4e "
                " slack=%.1f  in/out=%.3e  B_out=%.4f"
                % (r['eta'], r['LamHat_over_eta'], r['bound_leak'], r['bound_trivial'],
                   r['cert'], r.get('relL1', float('nan')), r.get('slack', float('nan')),
                   r['in_over_out'], r['bound_leak_outonly']))
        del G1, al, be, buf; gc.collect()
    log(f"L={L} ALL DONE in {time.time()-t_all:.1f}s")


if __name__ == '__main__':
    L = int(sys.argv[1]); nl = int(sys.argv[2])
    FRs = [float(x) for x in sys.argv[3].split(',')]
    nref = int(sys.argv[4]) if len(sys.argv) > 4 else 800
    mm = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    tag = sys.argv[6] if len(sys.argv) > 6 else ''
    run(L, nl, FRs, nref, mm, tag)
