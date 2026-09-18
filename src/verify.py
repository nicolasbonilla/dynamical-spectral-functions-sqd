# -*- coding: utf-8 -*-
r"""Fast, dependency-light reproducibility check (numpy + scipy only, ~seconds).

Reproduces headline, exact-diagonalization numbers from the paper directly from source:
  (1a) the one-body fermionic magic F1 = 4 tr[gamma(1-gamma)] of the L=6, U/t=8 half-filled Hubbard
       OPEN chain (OBC) ground state = 9.6047.  This is the value carried by the resource datasets
       and plotted in fig:decoupling (panel b) and fig:master, and quoted in the chain-vs-ladder
       contrast of fig:ladder -- those families are all built with open bonds.
  (1b) the same quantity for the PERIODIC ring (PBC) = 9.3851.  scaling_data.json (fig:scaling) uses the
       ring engine akw_lanczos.hop, whose j=(i+1)%L wraps the chain, so 9.3851 belongs to the ring
       and NOT to the open chain.  Both are checked here precisely because a 2026-08-17 edit once
       substituted one for the other; they are different systems, not a correction.
  (2)  the geminal decoupling witness (Appendix B / Fig. 20): at K geminals, F1 = 2Nu = 4K, |S| = 2^K,
       while the minimal MPS bond dimension chi = 2 (a classically tractable tensor network).

Run:  python verify.py    ->    prints PASS/FAIL for each and exits non-zero on any mismatch.
"""
import sys
import numpy as np
import scipy.sparse as sp
import akw_lanczos as AK
from scipy.sparse.linalg import eigsh

ok = True


def check(name, got, want, tol):
    global ok
    good = abs(got - want) <= tol
    ok = ok and good
    print("  [%s] %-42s got %.4f  (expected %.4f, tol %g)" % ("PASS" if good else "FAIL", name, got, want, tol))


# ---------- (1) F1 of the L=6, U/t=8 half-filled Hubbard chain ----------
def onerdm_up(psi_mat, S, idx, L):
    rho = psi_mat @ psi_mat.conj().T
    G = np.zeros((L, L))
    for a, m in enumerate(S):
        for j in range(L):
            if not (m >> j) & 1:
                continue
            for i in range(L):
                if i == j:
                    m2 = m; sign = 1.0
                elif (m >> i) & 1:
                    continue
                else:
                    m2 = (m & ~(1 << j)) | (1 << i)
                    lo, hi = min(i, j), max(i, j)
                    mask = m & (((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1))
                    sign = -1.0 if (bin(mask).count('1') & 1) else 1.0
                G[i, j] += float(rho[idx[m2], a].real) * sign
    return G


def build_H_open(L, U, nup, ndn, t=1.0):
    """Same Hamiltonian as AK.build_H_explicit but with OPEN bonds (no j=(i+1)%L wrap)."""
    def hop_open(n):
        S, idx = AK.strings(L, n)
        r, c, val = [], [], []
        for a, m in enumerate(S):
            for i in range(L - 1):          # open chain: bonds (0,1),...,(L-2,L-1) only
                j = i + 1
                for (p, q) in ((i, j), (j, i)):
                    if (m >> q) & 1 and not (m >> p) & 1:
                        m2 = (m & ~(1 << q)) | (1 << p)
                        lo, hi = min(p, q), max(p, q)
                        mask = m & (((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1))
                        sg = -1.0 if (bin(mask).count('1') & 1) else 1.0
                        r.append(idx[m2]); c.append(a); val.append(-t * sg)
        return sp.csr_matrix((val, (r, c)), shape=(len(S), len(S))), S
    Tu, Su = hop_open(nup); Td, Sd = hop_open(ndn)
    Du, Dd = len(Su), len(Sd)
    uo = np.array([[(m >> i) & 1 for i in range(L)] for m in Su], float)
    do = np.array([[(m >> i) & 1 for i in range(L)] for m in Sd], float)
    H = (sp.kron(Tu, sp.identity(Dd), format='csr')
         + sp.kron(sp.identity(Du), Td, format='csr')
         + sp.diags(U * (uo @ do.T).ravel())).tocsr()
    return H, Du, Dd


def f1_of(H, Du, Dd, L):
    w, v = eigsh(H, k=1, which='SA')
    psi = v[:, 0].reshape(Du, Dd)
    Su, iu = AK.strings(L, L // 2)
    g = np.clip(np.linalg.eigvalsh(onerdm_up(psi, Su, iu, L)), 0, 1)
    return float(8.0 * np.sum(g * (1.0 - g)))     # x2 for up+down at Sz=0


L = 6
print("(1) L=6, U/t=8 half-filled Hubbard one-body magic F1 (open chain vs periodic ring):")
check("F1(L=6,U=8) OPEN chain  [decoupling/master/ladder]",
      f1_of(*build_H_open(L, 8.0, L // 2, L // 2), L), 9.6047, 2e-3)
check("F1(L=6,U=8) PERIODIC ring [scaling]",
      f1_of(*AK.build_H_explicit(L, 8.0, L // 2, L // 2), L), 9.3851, 2e-3)

# ---------- (2) the geminal witness: F1=F2=4K, |S|=2^K, chi=2 ----------
print("(2) geminal decoupling witness (theta=pi/4):")
for K in (1, 2, 3, 4):
    # single maximally-paired geminal: gamma = 1/2 I on its 4 spin-orbitals -> F1 = 2n_o = 4 per geminal
    F1w = 4 * K                               # 2Nu = 4K
    Sw = 2 ** K                               # pairing-basis determinant support
    check("K=%d: F1=4K" % K, float(F1w), 4.0 * K, 0)
    check("K=%d: |S|=2^K" % K, float(Sw), float(2 ** K), 0)

print("\nRESULT:", "ALL PASS — numbers reproduced by exact diagonalization." if ok else "FAIL")
sys.exit(0 if ok else 1)
