# -*- coding: utf-8 -*-
"""ATK2 / attack 5: the E0 objection, measured against the RIGHT baseline.

C6 reports that with E0 estimated from a sampled subspace of the N-electron sector the
certified object blows up by x30-x40, and concludes that Sec.5 would be certifying an object
no shot-based pipeline produces.  Three things to check:

 (1) A COMMON shift cancels exactly:  if A and A_P are both referred to the same E0hat, the
     L1 difference is invariant (both curves translate by the same delta).  So the blow-up is
     an artefact of comparing A_P(E0hat) against A(E0_exact), i.e. against a reference the
     pipeline never has.  Measured, not asserted.
 (2) The blow-up is measured against the TRUE ERROR (1.3e-3), not against the CERTIFICATE
     (1.2036).  For a CERTIFICATE the question is whether the extra term is large compared to
     the bound, not compared to the error.
 (3) delta is a-posteriori BOUNDABLE from the residual norm of the sampled ground state
     (Weinstein / Temple), which costs ONE full-sector matvec -- the same cost class as
     Lambda_S itself.  So the hole is patchable, and by a computable quantity.
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
sys.path.insert(0, SRC); sys.path.insert(0, PUERTAS)
import akw_lanczos as AK, c0_lib as C
from scipy.sparse.linalg import eigsh
T0 = time.time()
def log(*a): print("[%6.1fs]" % (time.time()-T0), *a, flush=True)

def lorgrid(poles, eta, per=16.0, pad=None, ntail=800, tailmax=1e7):
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

def main(L, FR=0.85, eta=0.18, FR0s=(0.20, 0.40, 0.60, 0.80)):
    U = 8.0; t = 1.0; nup = nd = L//2
    H0, Du, Dd = AK.build_H_explicit(L, U, nup, nd, t); H0 = H0.tocsr()
    ev, V0 = eigsh(H0, k=2, which='SA')
    idx = np.argsort(ev); E0 = float(ev[idx[0]]); E1 = float(ev[idx[1]])
    psi0 = V0[:, idx[0]]
    gap = E1 - E0
    log("L=%d  N-sector dim=%d  E0=%.12f  E1=%.12f  gap=%.6f" % (L, H0.shape[0], E0, E1, gap))
    # (N+1)-sector objects, exact
    s = C.system(L, U, t); H = s['H']; phi = s['phi']; D = s['D']; n2 = float(phi@phi)
    assert abs(s['E0']-E0) < 1e-8, (s['E0'], E0)
    Em, Um = np.linalg.eigh(H.toarray()); eps = Em - E0; c = Um.T@phi
    mv = lambda x: H@x
    v = (phi/np.linalg.norm(phi)).astype(complex); wc = np.abs(v)**2
    for st in range(1, 18*16+1):
        v = C.krylov_expm(mv, v, 0.5/16, 6)
        if st % 16 == 0: wc += np.abs(v)**2
    order = np.argsort(wc)[::-1]
    k = int(round(FR*D)); Sset = np.sort(order[:k])
    HS = np.asarray(H[Sset][:, Sset].todense()); Es, Us = np.linalg.eigh(HS)
    epst = Es - E0; d = Us.T@phi[Sset]
    w, dw = lorgrid(np.concatenate([eps, epst]), eta)
    A = lor(eps, c*c, w, eta); AS = lor(epst, d*d, w, eta)
    relR0 = float(np.dot(np.abs(A-AS), dw))/n2
    log("  (N+1) FR=%.4f |S|=%d  rel-L1 (exact origin) = %.6e" % (FR, k, relR0))

    # ---------- (1) common shift cancels exactly
    for dl in (1e-3, 1e-2, 1e-1):
        As = lor(eps+dl, c*c, w, eta); ASs = lor(epst+dl, d*d, w, eta)   # BOTH shifted
        r = float(np.dot(np.abs(As-ASs), dw))/n2
        log("  (1) common shift delta=%-8.3g : rel-L1 = %.6e   ratio to unshifted = %.12f" % (dl, r, r/relR0))
    # ---------- (1b) one-sided shift (C6's object): only A_S moves
    print()
    log("  (1b) one-sided shift (A exact origin, A_S shifted) -- C6's object")
    for dl in (1e-4, 3.682e-4, 1e-3, 6.73e-4, 1e-2):
        ASs = lor(epst+dl, d*d, w, eta)
        r = float(np.dot(np.abs(A-ASs), dw))/n2
        law = (4/np.pi)*np.arctan(abs(dl)/(2*eta))
        log("      delta=%-10.4g rel-L1=%.6e  (x%.2f)   arctan law bound=%.4e   sum-bound=%.4e"
            % (dl, r, r/relR0, law, relR0+law))
    # ---------- (2) the right baseline
    return dict(L=L, E0=E0, E1=E1, gap=gap, relR0=relR0, D=D, k=k, H0=H0, psi0=psi0,
                Du=Du, Dd=Dd, eta=eta, n2=n2)

def e0_estimates(st, FR0s):
    H0 = st['H0']; psi0 = st['psi0']; E0 = st['E0']; gap = st['gap']; eta = st['eta']
    p = np.abs(psi0)**2; o = np.argsort(p)[::-1]; D0 = len(p)
    print()
    print("  (3) E0 from a Born-ranked sampled subspace of the N-electron sector, with the")
    print("      A-POSTERIORI residual bound (one full-sector matvec):")
    print("      %-7s %-8s %-12s %-12s %-12s %-12s %-12s" %
          ("FR0", "|S0|", "delta=E0h-E0", "residual r", "Temple bnd", "arctan(delta)", "arctan(bnd)"))
    rows = []
    for fr in FR0s:
        k0 = max(1, int(round(fr*D0))); S0 = np.sort(o[:k0])
        Hb = np.asarray(H0[S0][:, S0].todense())
        e, vv = np.linalg.eigh(Hb)
        E0h = float(e[0]); psih = np.zeros(D0); psih[S0] = vv[:, 0]
        res = float(np.linalg.norm(H0@psih - E0h*psih))
        delta = E0h - E0
        temple = res*res/max(gap - (E0h-E0), 1e-12)        # Temple: E0 >= E0h - r^2/(E1-E0h)
        law = (4/np.pi)*np.arctan(abs(delta)/(2*eta))
        lawb = (4/np.pi)*np.arctan(min(abs(temple), abs(res))/(2*eta))
        print("      %-7.2f %-8d %-12.4e %-12.4e %-12.4e %-12.4e %-12.4e" %
              (fr, k0, delta, res, temple, law, lawb))
        rows.append(dict(FR0=fr, k0=k0, delta=delta, res=res, temple=temple, law=law, lawb=lawb))
    return rows

if __name__ == '__main__':
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    FR = float(sys.argv[2]) if len(sys.argv) > 2 else 0.85
    st = main(L, FR=FR)
    rows = e0_estimates(st, (0.20, 0.40, 0.60, 0.80))
    json.dump(dict(L=L, FR=FR, relR0=st['relR0'], rows=rows),
              open(os.path.join(ADV, 'z_e0_L%d.json' % L), 'w'), indent=1)
