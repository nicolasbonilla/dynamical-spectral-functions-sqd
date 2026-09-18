# -*- coding: utf-8 -*-
r"""G2 (the programme's namesake): AMORTIZED configuration recovery. The paper's GFlowNet result is a
single-instance null. Here we test the actual open question: train a recovery model ONCE across a family
of instances and generalize to a HELD-OUT instance's under-sampled spectral head, versus the per-instance
classical selector -- leave-one-out, with error bars. Honest first data point.

Family: Hubbard chains (L,U) [exact ED, no pyscf/torch needed]. Task: rank the (N+1)-sector configurations
of the frontier-orbital response seed phi=c^dag_0|psi0>. The model is trained on the OTHER instances only
(instance-transferable features), the classical baseline (Epstein-Nesbet PT1 score) is per-instance and
untrained. Metric: spectral relative-L1 of A(omega) reconstructed in the top-k configs, at matched k.
"""
# --- DEPOSIT PATHS (repaired 2026-09-18, second pass) -----------------------
# Eighteen scripts were repaired earlier today because they hard-coded their output
# under '/w/', the working directory of the Docker container the published runs were
# made in.  THIS FILE WAS NOT AMONG THEM, and it was broken in a quieter way: it wrote
# to the RELATIVE path 'data/...', which lands in whatever directory the reader happens
# to be standing in, and raises FileNotFoundError from anywhere except the repository
# root.  The Data Availability Statement claims that every deposited computation script
# writes into the repository's own data/ directory; that sentence was FALSE for this
# file until now.  Repaired exactly like the other eighteen: the default is computed
# from THIS FILE's location, and --out overrides it.
# Paths only -- no physics and no computational default was changed here.
import os as _os, sys as _sys
_REPO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))


def _flag(name):
    """Value of a `--name VALUE` (or `-n VALUE`) command-line flag, else None."""
    a = _sys.argv[1:]
    for f in ('--' + name, '-' + name[0]):
        if f in a and a.index(f) + 1 < len(a):
            return a[a.index(f) + 1]
    return None


def _outpath(name, sub='data'):
    """Absolute path to write `name` to: --out if given, else <repo>/<sub>/<name>."""
    p = _os.path.abspath(_flag('out') or _os.path.join(_REPO, sub, name))
    _os.makedirs(_os.path.dirname(p), exist_ok=True)
    return p


def _repopath(name, sub='data'):
    """Absolute path <repo>/<sub>/<name>.  Never overridden: for a script's SECOND
    output, which --out (a single flag) cannot address unambiguously."""
    p = _os.path.abspath(_os.path.join(_REPO, sub, name))
    _os.makedirs(_os.path.dirname(p), exist_ok=True)
    return p
# ---------------------------------------------------------------------------
import json, time, itertools, numpy as np
# --- NUMPY_TRAPEZOID_BRIDGE ------------------------------------------------
# numpy 2.0 ADDED np.trapezoid and REMOVED np.trapz.  Files in this repository use
# both names, so without this bridge no single numpy version runs the whole deposit:
# numpy 1.x breaks the files that call trapezoid, numpy 2.x breaks the files that call
# trapz (this guardian included).  requirements.txt asks for numpy>=1.24; with the
# bridge that is true again.
if not hasattr(np, "trapezoid"):
    np.trapezoid = np.trapz          # numpy < 2.0
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid          # numpy >= 2.0
# ---------------------------------------------------------------------------
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
# --- DETERMINISTIC ARPACK START VECTOR (added 2026-09-18, second pass) ------
# eigsh() with no v0= lets ARPACK draw its own start vector from an UNSEEDED
# generator, so E0 converges to a slightly different point on every run (the last
# few digits move) and every quantity derived from it moves with it.  Measured, not
# hypothetical: three consecutive calls on the same matrix gave
# -2.0481308860914536 / ...504 / ...522, and in the leakage certificate a small
# subspace selection moved by 3.2%.  The eigenpair is the same to ARPACK's
# tolerance -- no physics changes -- but a deposit must be bit-reproducible.
# Deliberately NOT a numpy global seed: this touches only the ARPACK start vector.
# The seed 20260918 is the one used by scaling_lanczos.py and the certificate suite.
def _v0(n):
    """Fixed, dimension-dependent ARPACK start vector (never orthogonal to the GS)."""
    return np.random.default_rng(20260918).standard_normal(n)
# ---------------------------------------------------------------------------
from sklearn.ensemble import GradientBoostingRegressor
import gate1_ladder as G
t0=time.time(); log=lambda *a: print(f"[{time.time()-t0:6.1f}s]",*a,flush=True)

def chain_bonds(L): return [(i,i+1,1.0) for i in range(L-1)]
def popcount(x): return bin(x).count('1')

def instance(L,U,eta=0.15,K=20,dt=0.5,p=0):
    """Return per-config features, true importance, PT1 classical score, and A(omega) machinery."""
    bonds=chain_bonds(L); nup=L//2; ndn=L//2
    H,Su,iu,Sd,idd,Du,Dd=G.build_H(L,U,nup,ndn,bonds)
    w,v=eigsh(H,k=1,which='SA',v0=_v0(H.shape[0])); E0=float(w[0]); psi=v[:,0].reshape(Du,Dd)
    # seed phi = c^dag_{p,up}|psi0> in (nup+1,ndn)
    Su1,iu1=G.strings(L,nup+1); r=[];c=[];val=[]
    for a,m in enumerate(Su):
        if not (m>>p)&1:
            m2=m|(1<<p); sg=-1.0 if (popcount(m&((1<<p)-1))&1) else 1.0
            r.append(iu1[m2]); c.append(a); val.append(sg)
    Cd=sp.csr_matrix((val,(r,c)),shape=(len(Su1),Du)); phi=(Cd@psi).reshape(-1)
    H1,Su1b,iu1b,Sd1,idd1,Du1,Dd1=G.build_H(L,U,nup+1,ndn,bonds); nS=H1.shape[0]
    Harr=H1.toarray(); Em,Um=np.linalg.eigh(Harr); coef=Um.conj().T@phi
    poles=Em-E0; grid=np.linspace(poles.min()-1,poles.max()+1,600)
    def spec(pw,ww):
        A=np.zeros_like(grid)
        for a,b in zip(pw,ww): A+=b*(eta/np.pi)/((grid-a)**2+eta**2)
        return A
    A_ex=spec(poles,np.abs(coef)**2); nrm=np.trapz(A_ex,grid)
    # true importance: time-evolved-seed union weight (the response-relevant weight)
    wc=np.zeros(nS)
    for k in range(K+1):
        vk=Um@(np.exp(-1j*Em*k*dt)*coef); wc+=np.abs(vk)**2
    # seed config index (the (nup+1,ndn) det with max |phi|)
    seed_idx=int(np.argmax(np.abs(phi)))
    Hdiag=np.diag(Harr); Hs=Hdiag[seed_idx]
    conn=np.abs(Harr[:,seed_idx]); ediff=np.abs(Hdiag-Hs)
    pt1=conn**2/np.maximum(ediff,0.1)                 # Epstein-Nesbet PT1 (classical, per-instance)
    # instance-transferable features per config
    su_of=Su1b; nUp=Du1; nDn=Dd1
    feats=np.zeros((nS,5))
    # decode (up,dn) strings per flat index a = iu*Dd1 + idd (row-major of Du1 x Dd1)
    for a in range(nS):
        iu_=a//Dd1; idn_=a%Dd1; mu=Su1b[iu_]; md=Sd1[idn_]
        # doublons
        dbl=popcount(mu&md)
        # hamming from seed (in both spin strings)
        su_s=Su1b[seed_idx//Dd1]; sd_s=Sd1[seed_idx%Dd1]
        ham=popcount(mu^su_s)+popcount(md^sd_s)
        feats[a]=[ham, dbl, conn[a], ediff[a]/max(U,1.0), np.log10(pt1[a]+1e-12)]
    def relL1_topk(idxs):
        Ss=np.sort(np.array(idxs)); e,u=np.linalg.eigh(Harr[np.ix_(Ss,Ss)]); aa=u.conj().T@phi[Ss]
        return float(np.trapz(np.abs(spec(e-E0,np.abs(aa)**2)-A_ex),grid)/nrm)
    return dict(L=L,U=U,nS=nS,feats=feats,y=np.log10(wc+1e-14),pt1=pt1,true=wc,relL1=relL1_topk)

if __name__=='__main__':
    FAM=[(6,4.0),(6,8.0),(6,12.0),(8,4.0),(8,8.0),(8,12.0)]
    log("building family (exact ED per instance)...")
    inst=[instance(L,U) for (L,U) in FAM]
    for I in inst: log(f"  L={I['L']} U={I['U']:.0f}: (N+1) sector nS={I['nS']}")
    # choose a common evaluation subspace size: fraction f of each held-out sector
    fracs=[0.05,0.10,0.20,0.30]
    results=[]
    for ho in range(len(FAM)):
        tr=[i for i in range(len(FAM)) if i!=ho]
        Xtr=np.vstack([inst[i]['feats'] for i in tr]); ytr=np.concatenate([inst[i]['y'] for i in tr])
        # amortized model uses ONLY transferable raw features (cols 0..3), NOT the classical pt1 (col 4)
        model=GradientBoostingRegressor(n_estimators=200,max_depth=3,learning_rate=0.05,subsample=0.8,random_state=0)
        model.fit(Xtr[:,:4],ytr)
        I=inst[ho]; nS=I['nS']
        pred=model.predict(I['feats'][:,:4])          # amortized ranking (trained on OTHER instances)
        order_amort=np.argsort(pred)[::-1]
        order_class=np.argsort(I['pt1'])[::-1]         # per-instance classical (untrained)
        order_oracle=np.argsort(I['true'])[::-1]
        row={'held_out':f"L{I['L']}U{int(I['U'])}",'nS':nS}
        for f in fracs:
            k=max(2,int(f*nS))
            row[f'amort_{f}']=I['relL1'](order_amort[:k])
            row[f'class_{f}']=I['relL1'](order_class[:k])
            row[f'oracle_{f}']=I['relL1'](order_oracle[:k])
        results.append(row)
        log(f"held-out {row['held_out']}: "+" | ".join(
            f"f={f}: amort {row[f'amort_{f}']:.3f} vs class {row[f'class_{f}']:.3f} (oracle {row[f'oracle_{f}']:.3f})" for f in fracs))
    # leave-one-out summary: mean +- std of (amort - class) per fraction; positive = amortized WORSE
    log("=== LEAVE-ONE-OUT SUMMARY (rel-L1; lower is better) ===")
    summary={}
    for f in fracs:
        am=np.array([r[f'amort_{f}'] for r in results]); cl=np.array([r[f'class_{f}'] for r in results])
        d=am-cl
        summary[f]={'amort_mean':float(am.mean()),'amort_std':float(am.std()),
                    'class_mean':float(cl.mean()),'class_std':float(cl.std()),
                    'delta_mean':float(d.mean()),'delta_std':float(d.std())}
        verdict='amortized BEATS classical' if d.mean()<-0.005 else ('amortized ~ classical (null)' if abs(d.mean())<=0.02 else 'amortized WORSE')
        log(f"  f={f}: amortized {am.mean():.3f}+-{am.std():.3f}  classical {cl.mean():.3f}+-{cl.std():.3f}  "
            f"delta {d.mean():+.3f}+-{d.std():.3f}  -> {verdict}")
    json.dump({'note':'G2 amortized recovery: train-once-across-family vs per-instance classical (Epstein-Nesbet), leave-one-out over 6 Hubbard instances. metric=spectral rel-L1 at matched subspace fraction. delta=amort-class (positive => amortized worse).',
               'family':FAM,'per_fold':results,'summary':summary},open(_outpath('amortized_recovery.json'),'w'),indent=1)
    log("WROTE amortized_recovery.json")
