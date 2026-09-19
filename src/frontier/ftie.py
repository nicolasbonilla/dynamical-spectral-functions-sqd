# -*- coding: utf-8 -*-
r"""CONDITION 5 of the gate verdict, MEASURED, not assumed.

"Nada de determinantes de amplitud cero.  Por encima de FR = 1/2 + 1/L el ranking
 rellena con determinantes que el dispositivo nunca muestrearia.  Declara como los tratas."

What is actually true for the REPO protocol (K=18, dt=0.5, sub=16) is measured here:

 1. how many configs of the (N+1) sector have EXACTLY zero accumulated Born weight
    (the quantity the ranking sorts) -- if none, no determinant in S is unsampleable;
 2. how many configs are outside supp(phi) (zero GROUND-STATE amplitude) -- these are
    the ones the C0/VE witness was accused of padding with, and they are sampleable
    here because the ranking sorts the TIME-EVOLVED Born distribution;
 3. how big the argsort tie at the cut is, and
 4. the actual sensitivity: re-rank with a float64 accumulator (the repo uses float32)
    and with a reversed tie-break, and remeasure Lambda / w_S / the certificate.
"""
import os, sys, time, json, numpy as np
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


def born_order_f64(sysd, K=18, dt=0.5, sub=16, mkry=6):
    """Identical protocol, float64 accumulator instead of the repo's float32."""
    mv = sysd['mv']; phi = sysd['phi']
    v = (phi / np.linalg.norm(phi)).astype(complex)
    wc = (np.abs(v) ** 2).astype(np.float64)
    ddt = dt / sub
    for s in range(1, K * sub + 1):
        v = C.krylov_expm(mv, v, ddt, mkry)
        if s % sub == 0:
            wc += np.abs(v) ** 2
    return np.argsort(wc)[::-1], wc


def one(L, FRs, nl):
    sysd = F.build_system(L)
    phi, mv, D = sysd['phi'], sysd['mv'], sysd['D']
    n2 = float(phi @ phi)
    order32, wc32, _ = F.ranking(sysd)
    order64, wc64 = born_order_f64(sysd)
    supp = np.abs(phi) > 0
    log(f"L={L} D={D}: |supp(phi)|={int(supp.sum())} ({supp.sum()/D:.6f}; 1/2+1/L={0.5+1.0/L:.6f})")
    log(f"   configs with EXACTLY zero accumulated Born weight: float32 "
        f"{int((wc32==0).sum())}/{D}   float64 {int((wc64==0).sum())}/{D}")
    log(f"   smallest nonzero wc: f32 {wc32[wc32>0].min():.3e}   f64 {wc64[wc64>0].min():.3e}")
    res = dict(L=L, D=D, supp=int(supp.sum()), supp_frac=float(supp.sum()) / D,
               zero_wc_f32=int((wc32 == 0).sum()), zero_wc_f64=int((wc64 == 0).sum()),
               points=[])
    for FR in FRs:
        k = max(1, int(round(FR * D)))
        S32 = np.zeros(D, bool); S32[order32[:k]] = True
        S64 = np.zeros(D, bool); S64[order64[:k]] = True
        diff = int(np.count_nonzero(S32 != S64)) // 2
        out_of_supp = int(np.count_nonzero(S32 & ~supp))
        cut = wc32[order32[k - 1]]
        tied = int(np.count_nonzero(wc32 == cut))
        tied_in = int(np.count_nonzero(wc32[order32[:k]] == cut))
        r = dict(FR=FR, S=k, n_swapped_f32_vs_f64=diff, out_of_supp_in_S=out_of_supp,
                 tied_at_cut_total=tied, tied_at_cut_inside=tied_in)
        # remeasure the certificate on BOTH subspaces
        for name, mask in (('f32', S32), ('f64', S64)):
            phiS = phi * mask; nphiS2 = float(phiS @ phiS); wS = nphiS2 / n2
            V, al, be, nrm, _ = F.lanczos_S(mv, mask, phi, min(nl, k))
            G1 = F.gram_G1(mv, V)
            m = V.shape[1]; del V
            dv = F.depth_eval(al, be, G1, m, nphiS2, ETAS)
            for eta in ETAS:
                lk, li, lo = dv['per_eta'][float(eta)]
                B, triv, cert = F.certificate(wS, lk / np.sqrt(n2) / eta)
                r[f'{name}_w_S'] = wS
                r[f'{name}_eta{eta}'] = dict(Lambda_K=lk, B=B, cert=cert, trivial=triv)
            del G1
        d18 = (r['f32_eta0.18']['B'], r['f64_eta0.18']['B'])
        log(f"   FR={FR:.2f} |S|={k}: f32-vs-f64 subspaces differ in {diff} determinants; "
            f"{out_of_supp} of S lie outside supp(phi); tie at cut {tied_in}/{tied}. "
            f"w_S {r['f32_w_S']:.10f} vs {r['f64_w_S']:.10f} ; "
            f"B(eta=.18) {d18[0]:.6f} vs {d18[1]:.6f} (rel {abs(d18[0]-d18[1])/d18[0]:.2e})")
        res['points'].append(r)
    return res


if __name__ == '__main__':
    out = []
    out.append(one(6, [0.50, 0.70, 0.85, 0.90], 290))
    out.append(one(8, [0.50, 0.70, 0.85, 0.90], 500))
    json.dump(out, open(os.path.join(OUT, 'ties.json'), 'w'), indent=1, default=float)
    log("WROTE out/ties.json")
