# -*- coding: utf-8 -*-
"""Independent adversarial test of Theorem 1(iii) of P2 (dynamical spectral functions).

A_S is defined EXACTLY as in release/src/spectral_validation_max.py :
    H_S = H[S][:,S]  (i.e. P_S H P_S in the determinant basis)
    Em,Um = eigh(H_S);  aS = Um^dag phi[S]
    A_S(w) = sum_m |aS_m|^2 L_eta(w - (Em - E0))
so A_S is Rayleigh-Ritz on the coordinate subspace span(S), NOT a Krylov/Lanczos object.

Quantities computed per subspace S:
  w_S      = ||P_S phi||^2 / ||phi||^2               (captured Born weight)
  LHS      = || A - A_S ||_1                          (absolute, full real line)
  RHS_pap  = 2 ||phi||^2 (1 - w_S)                    (the paper's claim, ignoring C rho^{-K_S})
  K_S      = largest K with span{phi,...,H^K phi} subset span(S)   (-1 if phi not in span(S))
  Lam_S    = ( (eta/pi) int || Q H (z-H_S)^{-1} P_S phi ||^2 dw )^{1/2}   [leakage, closed form]
  Vnorm    = || P_S H (1-P_S) ||_2
  RHS_new  = (||phi||+||P_S phi||) * ( ||(1-P_S)phi|| + Lam_S/eta )   [derived Feshbach/resolvent bound]
  RHS_triv = ||phi||^2 (1 + w_S)                      (trivial cap)
  LB       = ||phi||^2 (1 - w_S)                      (rigorous LOWER bound)
"""
import numpy as np, scipy.sparse as sp, json, time, sys

t0 = time.time()
def log(*a): print("[%6.1fs]" % (time.time()-t0), *a, flush=True)

# ---------------------------------------------------------------- Hubbard builder
def build_hubbard(L=6, U=8.0, t_hop=1.0, pbc=False):
    M = 2*L; dim = 1 << M
    def c_op(p):
        r=[];c=[];d=[]
        for s in range(dim):
            if (s>>p)&1:
                sg=(-1)**bin(s&((1<<p)-1)).count('1')
                r.append(s&~(1<<p)); c.append(s); d.append(float(sg))
        return sp.csr_matrix((d,(r,c)),shape=(dim,dim))
    C=[c_op(p) for p in range(M)]; Cd=[c.T.conj() for c in C]
    Nn=np.array([bin(s).count('1') for s in range(dim)])
    Sz=np.array([sum(((s>>(2*i))&1)-((s>>(2*i+1))&1) for i in range(L)) for s in range(dim)])
    H=sp.csr_matrix((dim,dim))
    links = list(range(L-1)) + ([L-1] if pbc else [])
    for i in links:
        for spn in (0,1):
            a=2*i+spn; b=2*((i+1)%L)+spn
            H = H - t_hop*(Cd[a]@C[b] + Cd[b]@C[a])
    for i in range(L): H = H + U*(Cd[2*i]@C[2*i])@(Cd[2*i+1]@C[2*i+1])
    return H.tocsr(), C, Cd, Nn, Sz, M, dim

def hubbard_problem(L=6, U=8.0, pbc=False, seed_site=0):
    H,C,Cd,Nn,Sz,M,dim = build_hubbard(L,U,1.0,pbc)
    gi = np.where((Nn==L)&(Sz==0))[0]
    w,v = np.linalg.eigh(H[gi][:,gi].toarray()); E0=w[0]
    psi0=np.zeros(dim); psi0[gi]=v[:,0]
    phi_full = Cd[2*seed_site]@psi0
    si = np.where((Nn==L+1)&(Sz==1))[0]
    Hs = H[si][:,si].toarray()
    phi = np.asarray(phi_full[si]).ravel().astype(float)
    return Hs, phi, E0

# ---------------------------------------------------------------- generic machinery
class Spec:
    def __init__(self, Hs, phi, E0, eta):
        self.H=Hs; self.phi=phi; self.E0=E0; self.eta=eta
        self.n=Hs.shape[0]
        self.nrm2=float(phi@phi)                 # ||phi||^2
        self.E,self.V=np.linalg.eigh(Hs)
        self.coef=self.V.T@phi
        self.poles=self.E-E0
        self.wts=self.coef**2
        # integration grid: dense core + log tails, full real line
        lo,hi=self.poles.min()-1.0,self.poles.max()+1.0
        core=np.linspace(lo-2,hi+2,120001)
        tail=np.concatenate([-np.logspace(np.log10(abs(lo-2)+1e-9),5,4000)[::-1],
                              np.logspace(np.log10(abs(hi+2)+1e-9),5,4000)])
        self.grid=np.unique(np.concatenate([core,tail]))
        self.A=self.lorentz(self.poles,self.wts)
        self.normA=np.trapz(self.A,self.grid)
    def lorentz(self,pl,ww):
        A=np.zeros_like(self.grid)
        e=self.eta
        for a,b in zip(pl,ww):
            if b!=0.0: A+= b*(e/np.pi)/((self.grid-a)**2+e**2)
        return A
    def A_of_S(self,S):
        S=np.asarray(sorted(S))
        HS=self.H[np.ix_(S,S)]
        Em,Um=np.linalg.eigh(HS)
        aS=Um.T@self.phi[S]
        return self.lorentz(Em-self.E0, aS**2), Em, Um, aS, S
    def krylov_order(self, S, tol=1e-10):
        """largest K with span{phi,...,H^K phi} subset span(S); -1 if phi not in span(S)."""
        mask=np.zeros(self.n,bool); mask[np.asarray(sorted(S))]=True
        v=self.phi.copy(); K=-1
        for j in range(0, self.n+1):
            nv=np.linalg.norm(v)
            if nv==0: return K            # exhausted
            leak=np.linalg.norm(v[~mask])/nv
            if leak>tol: return K
            K=j
            v=self.H@v
        return K
    def metrics(self,S):
        A_S,Em,Um,aS,S = self.A_of_S(S)
        mask=np.zeros(self.n,bool); mask[S]=True
        phiS=np.where(mask,self.phi,0.0); phiQ=self.phi-phiS
        wS=float(phiS@phiS)/self.nrm2
        LHS=float(np.trapz(np.abs(self.A-A_S),self.grid))
        # leakage Lambda_S in closed form
        # g_m = Q H |m~>  ; c_m = <m~|phi>
        Uf=np.zeros((self.n,len(Em))); Uf[S,:]=Um
        G=self.H@Uf                     # H |m~>
        G[mask,:]=0.0                   # project with Q
        Gram=G.T@G
        c=aS
        dE=Em[None,:]-Em[:,None]        # dE[m',m] = Em'... careful below
        e=self.eta
        # factor F[m,m'] = (4 eta^2 + 2 i eta (E_{m'}-E_m)) / ((E_{m'}-E_m)^2+4 eta^2)
        D=Em[None,:]-Em[:,None]         # D[m,m'] = E_{m'} - E_m
        F=(4*e*e + 2j*e*D)/(D*D+4*e*e)
        # Lam^2 = sum_{m,m'} c_m c_{m'} Gram[m',m] F[m,m']
        M2=(c[:,None]*c[None,:])*Gram.T*F
        Lam2=float(np.real(M2.sum()))
        Lam=np.sqrt(max(Lam2,0.0))
        Vnorm=float(np.linalg.norm(self.H[np.ix_(S,~mask)],2)) if (~mask).any() else 0.0
        nphi=np.sqrt(self.nrm2); nphiS=np.sqrt(wS*self.nrm2); nphiQ=np.linalg.norm(phiQ)
        return dict(nS=int(len(S)), frac=len(S)/self.n, wS=wS, LHS=LHS,
                    RHS_pap=2*self.nrm2*(1-wS), K_S=self.krylov_order(S),
                    Lam=Lam, Vnorm=Vnorm,
                    RHS_new=(nphi+nphiS)*(nphiQ+Lam/e),
                    RHS_sqrt_only=(nphi+nphiS)*nphiQ,
                    RHS_crude=(nphi+nphiS)*(nphiQ+Vnorm*nphiS/e),
                    RHS_triv=self.nrm2*(1+wS), LB=self.nrm2*(1-wS),
                    relL1=LHS/self.normA)
