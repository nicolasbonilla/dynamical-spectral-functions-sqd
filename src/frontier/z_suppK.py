# -*- coding: utf-8 -*-
"""ATK2 / the sharp reframing of 'the frontier'.

Lambda_S = 0  <=>  Q H K(H_S,phi_S) = 0  <=>  span(S) contains K(H, phi_S), the Krylov space of
the RETAINED part of the probe; at full captured weight (phi_S = phi) that is K = K(H,phi).
(Corrected 2026-09-25: the earlier statement used K(H,phi) for every S, which fails for w_S<1.)

WARNING, 2026-09-25: the dense projector below is contaminated where eigenvalues of the two
reflection sectors (j -> -j about the seed site) are degenerate: eigh mixes them, and the
reflection-forbidden determinants -- which vanish exactly in every vector of K -- are counted.
At L = 8 it reports 3920 of 3920; the true support is 3902.  The support quoted in the paper
is from z_suppK_sym.py, which projects onto the reflection sector of phi.  Because span(S) is a COORDINATE subspace, this happens iff
S contains the COORDINATE SUPPORT of K,  supp(K) = { i : (P_K)_ii > 0 }.

So  FR_exact := |supp(K)|/D  is a HARD, ranking-independent threshold:
  FR >= FR_exact (with the right determinants) -> A_S = A exactly, certificate == 0
  below it, the certificate is bounded away from 0.
This is the object C6 is estimating by interpolation and calling 'the frontier'.
"""
import os, sys, numpy as np
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
PUERTAS = HERE
sys.path.insert(0, PUERTAS)
import c0_lib as C
for L in [int(x) for x in sys.argv[1:]] or [6, 8]:
    s = C.system(L); H = s['H']; phi = s['phi']; D = s['D']
    Em, Um = np.linalg.eigh(H.toarray())
    c = Um.T @ phi; nph = np.linalg.norm(phi)
    print("L=%d  D=%d  |supp(phi)|=%d (%.4f)" % (L, D, int(np.sum(np.abs(phi) > 0)), np.mean(np.abs(phi) > 0)))
    for tol in (1e-8, 1e-10, 1e-12):
        live = np.abs(c) > tol*nph
        # P_K = projector onto span of the live eigenvectors, but degenerate energies must be
        # collapsed: K only contains the COMPONENT of phi inside each degenerate eigenspace.
        E = Em[live]; Uv = Um[:, live]; cv = c[live]
        # group degenerate energies and form the single Krylov direction in each group
        cols = []
        i = 0
        while i < len(E):
            j = i
            while j+1 < len(E) and E[j+1]-E[i] < 1e-9: j += 1
            v = Uv[:, i:j+1] @ cv[i:j+1]
            n = np.linalg.norm(v)
            if n > tol*nph: cols.append(v/n)
            i = j+1
        K = np.array(cols).T                      # D x dimK, orthonormal columns
        diag = (K*K).sum(1)                       # (P_K)_ii
        for st in (1e-12, 1e-14, 1e-16):
            sup = int(np.sum(diag > st))
            print("    tol=%.0e dimK=%-5d | supp(K) at %.0e : %5d / %d = %.5f" %
                  (tol, K.shape[1], st, sup, D, sup/D))
        print("    tol=%.0e : smallest 8 nonzero (P_K)_ii = %s" %
              (tol, np.array2string(np.sort(diag[diag > 1e-18])[:8], precision=2)))
    print()
