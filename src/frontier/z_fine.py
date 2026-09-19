# -*- coding: utf-8 -*-
"""ATK2 / attack 4: is FR*_5 a MEASUREMENT or an INTERPOLATION ARTEFACT?

C6 reports FR*_5 by log-linear interpolation on the FR grid.  At L=6 the crossing bracket is
[0.98, 0.99], across which R5 falls from 3.2e-1 to 3.8e-14 -- THIRTEEN orders of magnitude,
because at FR=0.99 the subspace is 297/300 and Lambda collapses.  Interpolating log-linearly
across that is meaningless.  Here I evaluate EVERY integer |S| in the bracket exactly.

Own implementation of Lambda_S:  I re-derive the eta-average myself.
   g(w) = sum_n d_n g_n /(w - e_n + i eta),   g_n = Q H |n~>,  e_n = Et_n - E0
   (eta/pi) int dw 1/[(w-a+i eta)(w-b-i eta)] = 2 i eta / ((b-a) + 2 i eta)   [residue at w=b+i eta]
   => Lambda^2 = Re sum_{n,n'} d_n d_n' <g_n',g_n> * 2 i eta/((e_n' - e_n) + 2 i eta)
and I check it against a DIRECT omega quadrature of ||g(w)||^2 at L=6.
"""
import os, sys, time, json, numpy as np
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
# --- deposit paths (repaired 2026-09-19; see docs/C3_FRONTIER.md) ----------
SRC = os.path.dirname(HERE)                       # release/src
ROOT = os.path.dirname(SRC)                       # repository root
OUT = os.environ.get('C3_OUT') or os.path.join(ROOT, 'data', 'c3_frontier')
TMP = os.environ.get('C3_TMP') or os.path.join(ROOT, 'build')
# --------------------------------------------------------------------------
ADV = os.path.join(OUT, 'adversarial')
if not os.path.isdir(ADV):
    os.makedirs(ADV)
PUERTAS = HERE
sys.path.insert(0, PUERTAS)
import c0_lib as C
T0 = time.time()
def log(*a): print("[%6.1fs]" % (time.time()-T0), *a, flush=True)

def ranking(H, phi, K=18, dt=0.5, sub=16, m=6):
    mv = lambda x: H @ x
    v = (phi/np.linalg.norm(phi)).astype(complex)
    wc = np.abs(v)**2
    for s in range(1, K*sub+1):
        v = C.krylov_expm(mv, v, dt/sub, m)
        if s % sub == 0: wc += np.abs(v)**2
    return np.argsort(wc)[::-1], wc

def grid(poles, eta, per=14.0, pad=None, ntail=800, tailmax=1e7):
    if pad is None: pad = max(4.0, 40.0*eta)
    lo = poles.min()-pad; hi = poles.max()+pad
    core = np.linspace(lo, hi, int(np.ceil((hi-lo)/(eta/per)))+1)
    u = np.linspace(0, np.log(tailmax), ntail+1)[1:]
    w = np.concatenate([(lo-np.expm1(u))[::-1], core, hi+np.expm1(u)])
    dw = np.empty_like(w); dw[1:-1] = .5*(w[2:]-w[:-2]); dw[0] = .5*(w[1]-w[0]); dw[-1] = .5*(w[-1]-w[-2])
    return w, dw

def lor(poles, wts, w, eta, chunk=256):
    out = np.zeros_like(w); pre = eta/np.pi
    for i in range(0, len(poles), chunk):
        p = poles[i:i+chunk]; g = wts[i:i+chunk]
        out += ((g*pre)[:, None]/((w[None, :]-p[:, None])**2+eta**2)).sum(0)
    return out

def evaluate(L, Em, Um, H, phi, E0, Sset, etas, direct_check=False):
    D = len(phi); n2 = float(phi@phi); nphi = np.sqrt(n2)
    eps = Em-E0; c = Um.T@phi
    mask = np.zeros(D, bool); mask[Sset] = True
    phiS = phi*mask; phiQ = phi-phiS
    wS = float(phiS@phiS)/n2; nS = np.linalg.norm(phiS); nQ = np.linalg.norm(phiQ)
    HS = np.asarray(H[Sset][:, Sset].todense())
    Es, Us = np.linalg.eigh(HS)
    epst = Es-E0; d = Us.T@phi[Sset]
    # G = Q H |n~>  (D x m)
    Ufull = np.zeros((D, len(Es))); Ufull[Sset, :] = Us
    G = H@Ufull; G[Sset, :] = 0.0
    Gram = G.T@G
    del Ufull
    out = []
    for eta in etas:
        dE = epst[:, None]-epst[None, :]                     # e_n' - e_n  (rows n', cols n)
        K = 2j*eta/(dE+2j*eta)
        lam2 = float(((d[:, None]*d[None, :])*Gram*K).real.sum())
        Lam = np.sqrt(max(lam2, 0.0))
        w, dw = grid(np.concatenate([eps, epst]), eta)
        A = lor(eps, c*c, w, eta); AS = lor(epst, d*d, w, eta)
        R0 = float(np.dot(np.abs(A-AS), dw))
        R5 = (nphi+nS)*(nQ+Lam/eta)
        rec = dict(k=int(len(Sset)), FR=len(Sset)/D, eta=float(eta), w_S=wS,
                   R0=R0, relR0=R0/n2, R5=R5, relR5=R5/n2, Lam=float(Lam),
                   Lamhat_over_eta=float(Lam/nphi/eta), phiQ_rel=float(nQ/nphi),
                   trivial_rel=float(1.0+wS))
        if direct_check:
            # direct omega quadrature of ||g(w)||^2 -- no closed form
            wq, dwq = grid(epst, eta, per=20.0, ntail=1200)
            acc = 0.0
            for i in range(0, len(wq), 512):
                ww = wq[i:i+512]
                coef = d[None, :]/(ww[:, None]-epst[None, :]+1j*eta)     # (nw, m)
                gg = coef@G.T                                            # (nw, D)
                acc += float(((np.abs(gg)**2).sum(1)*dwq[i:i+512]).sum())
            lam2q = (eta/np.pi)*acc
            rec['Lam_quad'] = float(np.sqrt(max(lam2q, 0.0)))
            rec['Lam_quad_over_closed'] = rec['Lam_quad']/max(Lam, 1e-300)
        out.append(rec)
    del G, Gram
    return out

if __name__ == '__main__':
    L = int(sys.argv[1])
    ks = [int(x) for x in sys.argv[2].split(',')]
    etas = [float(x) for x in sys.argv[3].split(',')]
    chk = len(sys.argv) > 4 and sys.argv[4] == 'check'
    s = C.system(L); H = s['H']; phi = s['phi']; E0 = s['E0']; D = s['D']
    order, wc = ranking(H, phi)
    Em, Um = np.linalg.eigh(H.toarray())
    log("L=%d D=%d  ready" % (L, D))
    res = []
    for k in ks:
        Sset = np.sort(order[:k])
        t = time.time()
        rows = evaluate(L, Em, Um, H, phi, E0, Sset, etas, direct_check=chk)
        for r in rows:
            r['L'] = L
            log("  |S|=%4d FR=%.5f eta=%.2f | w_S=%.12f relR0=%.6e relR5=%.6e Lamhat/eta=%.6e%s"
                % (r['k'], r['FR'], r['eta'], r['w_S'], r['relR0'], r['relR5'], r['Lamhat_over_eta'],
                   ("  Lam_quad/closed=%.9f" % r['Lam_quad_over_closed']) if 'Lam_quad_over_closed' in r else ""))
        res += rows
        log("     (%.1fs)" % (time.time()-t))
    json.dump(res, open(os.path.join(ADV, 'z_fine_L%d.json' % L), 'w'), indent=1)
