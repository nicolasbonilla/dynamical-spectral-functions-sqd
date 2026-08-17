# -*- coding: utf-8 -*-
r"""THE PROVABLE WITNESS for the decoupling thesis.
Oh, Oszmaniec, Reardon-Smith, Zimboras (arXiv:2604.26813, Apr 2026) PROVE that antisymmetrized products
of strongly-orthogonal geminals (APSG / perfect-pairing GVB) are CLASSICALLY SIMULABLE under free-fermion
dynamics (transition amplitudes, overlaps, arbitrary-weight number correlators, via a Pfaffian-kernel that
collapses the exponential determinant sum to one polynomial coefficient).

We compute, on that SAME provably-simulable family, the resources our paper tracks:
  * one-body fermionic magic  FAF_1 = 2 N_u   (must be LARGE and equal to 2 N_u — verifies our identity)
  * higher-order antiflatness  F_2            (where the genuine non-Gaussian content lives)
  * determinant support  |S|                  (raw config count in the pairing basis — LARGE, ~2^K)
  * bond dimension  chi                        (Schmidt rank across the central geminal cut — SMALL, O(1))

=> a family with arbitrarily large one-body magic AND large |S|, yet SMALL entanglement and PROVABLY zero
   sampling hardness. This upgrades the decoupling from an empirical correlation to a THEOREM WITH A WITNESS:
   neither FAF_1 nor raw |S| certifies cost; the consistent governor (bond dimension) stays O(1) exactly
   where the state is provably easy. Writes apsg_witness.json.
"""
import json, numpy as np
import n19_suite as N   # reuse jw_ops / faf machinery

def fk_all(psi, C, Cd, M, kmax=3):
    """Full Williamson spectrum of the Majorana covariance -> F_k = M - 0.5*sum(sv^{2k})."""
    G=[]
    for p in range(M):
        G.append(C[p]@psi+Cd[p]@psi); G.append(-1j*(C[p]@psi-Cd[p]@psi))
    G=np.array(G); Mc=np.zeros((2*M,2*M))
    for a in range(2*M):
        for b in range(a+1,2*M):
            v=-1j*np.vdot(G[a],G[b]); Mc[a,b]=v.real; Mc[b,a]=-v.real
    sv=np.linalg.svd(Mc,compute_uv=False)              # 2M singular values = Williamson eigenvalues (paired)
    return {k: float(M-0.5*np.sum(sv**(2*k))) for k in range(1,kmax+1)}

def n_unpaired(psi, C, Cd, ncas):
    """N_u = sum_a n_a(2-n_a) from spatial natural occupations (gauge-invariant multireference index)."""
    # spin-orbital occupations n_p = <c^dag_p c_p>; spatial n_a = n_{2a}+n_{2a+1}
    Nu=0.0
    for a in range(ncas):
        na=np.vdot(psi,Cd[2*a]@C[2*a]@psi).real + np.vdot(psi,Cd[2*a+1]@C[2*a+1]@psi).real
        Nu+=na*(2-na)
    return float(Nu)

def support(psi, thr=1.0-1e-9):
    w=np.sort(np.abs(psi)**2)[::-1]; c=np.cumsum(w)/w.sum(); return int(np.searchsorted(c,thr)+1)

def chi_cut(psi, M, cut):
    """Schmidt rank (bond dimension) across a spin-orbital bipartition at position `cut`."""
    A=psi.reshape(1<<cut, 1<<(M-cut))
    s=np.linalg.svd(A,compute_uv=False); s=s[s>1e-10]; p=s**2/np.sum(s**2)
    return int(np.sum(p>1e-12)), float(-(p*np.log2(p)).sum())

def apsg_state(thetas):
    """Product of K strongly-orthogonal 2e singlet (perfect-pairing) geminals.
    Geminal i on spatial orbitals (2i=bonding, 2i+1=antibonding), spin-orbitals 4i..4i+3:
        |g_i> = cos th |b_up b_dn> + sin th |a_up a_dn>.
    Strong orthogonality = disjoint orbital blocks => APSG. Total M=4K spin-orbitals, N=2K electrons."""
    K=len(thetas); M=4*K; C,Cd=N.get_ops(2*K)     # 2K spatial orbitals -> 2*(2K)=4K spin-orbitals
    psi=np.zeros(1<<M,complex); psi[0]=1.0
    for i,th in enumerate(thetas):
        bu,bd,au,ad=4*i,4*i+1,4*i+2,4*i+3
        gi=np.cos(th)*(Cd[bu]@Cd[bd]@psi)+np.sin(th)*(Cd[au]@Cd[ad]@psi)
        psi=gi
    psi/=np.linalg.norm(psi); return psi,C,Cd,M,2*K

if __name__=='__main__':
    out={'note':'APSG/perfect-pairing geminal family — provably simulable (Oh et al. 2604.26813). '
                'FAF_1=2N_u large, |S| ~2^K large, chi O(1): the decoupling witness.','rows':[]}
    print(f"{'K':>2} {'theta':>6} {'FAF_1':>8} {'2*N_u':>8} {'|id ok|':>8} {'F_2':>8} {'|S|':>6} {'chi_geom':>9} {'S_ent':>7}")
    for K in (1,2,3,4):
        th=[np.pi/4]*K                      # maximal pairing correlation
        psi,C,Cd,M,Ne=apsg_state(th)
        F=fk_all(psi,C,Cd,M,kmax=3); Nu=n_unpaired(psi,C,Cd,2*K)
        S=support(psi); chi,Sent=chi_cut(psi,M,cut=2*K)   # cut down the middle (between geminal blocks)
        ok=abs(F[1]-2*Nu)
        out['rows'].append(dict(K=K,theta='pi/4',M=M,Ne=Ne,FAF1=F[1],twoNu=2*Nu,idcheck=ok,
                                F2=F[2],F3=F[3],Sdet=S,chi=chi,Sent=Sent))
        print(f"{K:>2} {'pi/4':>6} {F[1]:8.4f} {2*Nu:8.4f} {ok:8.1e} {F[2]:8.4f} {S:6d} {chi:9d} {Sent:7.3f}")
    # theta sweep at K=3: FAF_1 scans while the family stays provably simulable
    print("\ntheta sweep at K=3 (family stays APSG=provably simulable throughout):")
    print(f"{'theta/pi':>8} {'FAF_1':>8} {'2*N_u':>8} {'F_2':>8} {'|S|':>6} {'chi':>5}")
    for t in (0.02,0.10,0.25,0.40,0.50):
        th=[t*np.pi]*3; psi,C,Cd,M,Ne=apsg_state(th)
        F=fk_all(psi,C,Cd,M,kmax=2); Nu=n_unpaired(psi,C,Cd,6); S=support(psi); chi,_=chi_cut(psi,M,6)
        out['rows'].append(dict(K=3,theta=f'{t}pi',FAF1=F[1],twoNu=2*Nu,F2=F[2],Sdet=S,chi=chi))
        print(f"{t:8.2f} {F[1]:8.4f} {2*Nu:8.4f} {F[2]:8.4f} {S:6d} {chi:5d}")
    json.dump(out,open('/w/apsg_witness.json','w'),indent=1)
    print("\nWROTE apsg_witness.json")
