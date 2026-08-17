# -*- coding: utf-8 -*-
"""Max-rigor checks: first-moment sum rule, peak dispersion, Haydock convergence."""
import json, numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
trapz=np.trapezoid
import akw_lanczos as ak

d=json.load(open('/w/akw_lanczos_L12.json')); wg=np.array(d['wg'])
ks=np.array(sorted(float(x) for x in d['A'].keys()))
M=np.array([d['A'][f'{k:.5f}'] for k in ks])
L=12; U=8.0; eta=0.18

print("=== FIRST-MOMENT sum rule:  int w A(k,w) dw  vs exact  -2t cos(k)  (m=up) ===")
print("    (exact Hubbard f-sum: M1(k)=eps_k-mu+U<n_dn> = -2cos k at half-filling, mu=U/2)")
maxerr=0
for k in ks:
    A=M[int(np.argmin(np.abs(ks-k)))]
    m1=trapz(wg*A,wg); ref=-2*np.cos(k*np.pi); e=m1-ref; maxerr=max(maxerr,abs(e))
    print(f"  k/pi={k:.3f}  M1_data={m1:+.4f}   -2cos={ref:+.4f}   diff={e:+.4f}")
print(f"  max |diff| = {maxerr:.4f}  (small = f-sum rule satisfied; residual from +-9t truncation)")

print("\n=== dominant-peak dispersion from exact data ===")
for k in ks:
    A=M[int(np.argmin(np.abs(ks-k)))]
    neg=wg<0; pos=wg>0
    wr=wg[neg][np.argmax(A[neg])]; wa=wg[pos][np.argmax(A[pos])]
    print(f"  k/pi={k:.3f}  removal-peak={wr:+.2f}t  addition-peak={wa:+.2f}t")

print("\n=== HAYDOCK convergence at k_F=pi/2 (nl=150 used vs nl=300 reference) ===")
nup=nd=L//2
Hexpl,Du,Dd=ak.build_H_explicit(L,U,nup,nd)
E0,V0=eigsh(Hexpl,k=1,which='SA',maxiter=5000,tol=1e-10); E0=float(E0[0]); psi0=V0[:,0]
Psi=psi0.reshape(Du,Dd)
Ha,_,_,_,_,Dua,_=ak.sector_H(L,U,nup+1,nd)
Hr,_,_,_,_,Dur,_=ak.sector_H(L,U,nup-1,nd)
cdU=ak.cdag_map(L,nup); cU=ak.cdag_map(L,nup-1)
k=np.pi/2; ph=np.exp(1j*k*np.arange(L))/np.sqrt(L)
add=np.zeros((Dua,Dd),dtype=complex)
for j in range(L): add+=ph[j]*(cdU[j]@Psi)
rem=np.zeros((Dur,Dd),dtype=complex)
for j in range(L): rem+=np.conj(ph[j])*(cU[j].T@Psi)
mu=U/2.0; zt=(wg+mu)+1j*eta
def spec(nl):
    Gp=ak.haydock(lambda x:ak._apply(Ha,x)-E0*x, add.ravel(), nl, zt)
    Gm=ak.haydock(lambda x:E0*x-ak._apply(Hr,x), rem.ravel(), nl, zt)
    return -(1/np.pi)*np.imag(Gp+Gm)
A150=spec(150); A300=spec(300)
rel=trapz(np.abs(A150-A300),wg)/trapz(np.abs(A300),wg)
print(f"  rel-L1[A_nl150 - A_nl300] at k_F = {rel:.2e}   (tiny => nl=150 is converged)")
print(f"  sum rule nl=150: {trapz(A150,wg):.4f}   nl=300: {trapz(A300,wg):.4f}")
