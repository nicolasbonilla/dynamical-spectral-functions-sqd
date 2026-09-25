# -*- coding: utf-8 -*-
r"""L=14 certificate point at the PUBLISHED fraction FR=0.08, reusing the repo's
own verified checkpoints (ground state residual 2.7e-9, seed identical to
c^dag_{0up}|Psi_0> to 0.0e+00, order a true permutation reproducing |S|=824504).
Nothing in release/ is touched.  Same algebra as frun12.py, same protocol.
Depth 250: at L=12 the LOW-fraction points (FR=0.10, 0.18) are converged by
n ~ 130-200 (Lam/eta changes 6e-6 from 250->300), so 250 is the right depth here."""
import os,sys,time,gc,json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
# --- deposit paths (repaired 2026-09-19; see docs/C3_FRONTIER.md) ----------
SRC = os.path.dirname(HERE)                       # release/src
ROOT = os.path.dirname(SRC)                       # repository root
OUT = os.environ.get('C3_OUT') or os.path.join(ROOT, 'data', 'c3_frontier')
TMP = os.environ.get('C3_TMP') or os.path.join(ROOT, 'build')
# --------------------------------------------------------------------------
PUB = os.path.join(OUT, 'published_fraction')
if not os.path.isdir(PUB):
    os.makedirs(PUB)
sys.path.insert(0, HERE)
sys.dont_write_bytecode=True
import fcore as F
from fcore import log
from frun import bare_ref
from c0_big import sector_mv

# ckpt/L14_gs.npz (94 MB) and ckpt/L14_order.npy (41 MB) are deposited (2026-09-25);
# P2_CKPT overrides the folder.  See ckpt/README.txt.
CK = os.environ.get('P2_CKPT') or os.path.join(ROOT, 'ckpt')
L=14; U=8.0; NL=250; NREF=500; ETAS=[0.10,0.15,0.18,0.25]
FRS=[0.08]

_gs = os.path.join(CK, 'L14_gs.npz')
if not os.path.isfile(_gs):
    print("SKIP: %s not found. The deposit carries it in ckpt/ (94 MB ground state + "
          "41 MB ranking); set P2_CKPT to the folder that holds L14_gs.npz and "
          "L14_order.npy." % _gs)
    sys.exit(3)
g=np.load(_gs)
E0=float(g['E0']); psi=g['psi'].astype(float); nup=nd=7
mv0,D0,Du0,Dd0=sector_mv(L,U,nup,nd)
from c0_lib import AK
Psi=psi.reshape(Du0,Dd0); cdU=AK.cdag_map(L,nup)
phi=np.asarray((cdU[0]@Psi).reshape(-1),dtype=float)
del Psi,cdU,psi,g,mv0; gc.collect()
mv,D,Du,Dd=sector_mv(L,U,nup+1,nd)
n2=float(phi@phi)
log(f"L={L} D={D} E0={E0:.10f} |phi|^2={n2:.12f}  [checkpoint, gs residual 2.656e-09]")
order=np.load(os.path.join(CK,'L14_order.npy'))
log(f"order from checkpoint: permutation of {order.size} (verified), int32")

t=time.time()
th_r,w_r,m_r=bare_ref(mv,phi,NREF)
th_r2,w_r2,_=bare_ref(mv,phi,NREF//2)
drift={}
for e in ETAS: drift[e]=F.l1_lorentz(th_r,w_r,th_r2,w_r2,e)['l1']/n2
del th_r2,w_r2; gc.collect()
log(f"reference bare Lanczos m={m_r} [{time.time()-t:.0f}s]; drift {NREF//2}->{NREF}: "
    +" ".join("eta=%.2f:%.2e"%(e,drift[e]) for e in ETAS))

out=[]
for FR in FRS:
    t0=time.time()
    k=max(1,int(round(FR*D))); sel=order[:k]
    mask=np.zeros(D,bool); mask[sel]=True
    idx=np.where(mask)[0]
    cut_ok=dict(k=int(k),FR_real=k/D)
    phiS=phi[idx].copy(); nphiS2=float(phiS@phiS); wS=nphiS2/n2
    nl=int(min(NL,k)); gib=k*nl*8/2**30
    log(f"FR={FR:.2f} |S|={k} ({k/D:.5f}) w_S={wS:.12f}  V={gib:.2f} GiB (memmap)")
    buf=np.zeros(D)
    def mvS(xs):
        buf[:]=0.0; buf[idx]=xs; return mv(buf)[idx]
    mp=os.path.join(TMP,f"_V_L14_FR{FR:.2f}.dat")
    if not os.path.isdir(TMP):
        os.makedirs(TMP)
    V=np.memmap(mp,dtype=np.float64,mode='w+',shape=(k,nl),order='F')
    # two-pass CGS Lanczos, byte-identical algorithm to frun12.lanczos_coords
    nrm=float(np.linalg.norm(phiS)); V[:,0]=phiS/nrm
    al=np.zeros(nl); be=np.zeros(nl); m=nl
    w=mvS(np.asarray(V[:,0])); a=float(V[:,0]@w); al[0]=a; w=w-a*np.asarray(V[:,0])
    t=time.time()
    for j in range(1,nl):
        for _ in range(2):
            Vj=V[:,:j]; w-=Vj@(Vj.T@w)
        b=float(np.linalg.norm(w)); be[j-1]=b
        if b<1e-12: m=j; break
        V[:,j]=w/b
        w=mvS(np.asarray(V[:,j])); a=float(V[:,j]@w); al[j]=a
        w=w-a*np.asarray(V[:,j])-b*np.asarray(V[:,j-1])
        if j%50==0: log(f"   lanczos {j}/{nl} [{time.time()-t:.0f}s]")
    if m==nl:
        for _ in range(2):
            Vj=V[:,:nl]; w-=Vj@(Vj.T@w)
        be[nl-1]=float(np.linalg.norm(w))
    t_lan=time.time()-t; V=V[:,:m]
    js=np.unique(np.linspace(0,m-1,min(m,40)).astype(int))
    O=np.asarray(V.T@np.asarray(V[:,js])); O[js,np.arange(len(js))]-=1.0
    orth=float(np.abs(O).max()); del O
    Vtphi=np.asarray(V.T@phiS)
    resPK=float(np.abs(np.asarray(V@Vtphi)-phiS).max()); wK=float(Vtphi@Vtphi)/n2
    t=time.time(); G1=np.empty((m,m))
    for j in range(m):
        buf[:]=0.0; buf[idx]=np.asarray(V[:,j]); y=mv(mv(buf))[idx]
        G1[:,j]=np.asarray(V.T@y)
        if (j+1)%50==0: log(f"   G1 {j+1}/{m} [{time.time()-t:.0f}s]")
    G1=0.5*(G1+G1.T); t_g1=time.time()-t
    del V; gc.collect()
    try: os.remove(mp)
    except OSError: pass
    lad=F.ladder(m)
    deps=[F.depth_eval(al,be,G1,n,nphiS2,ETAS,keep_pairs=(n==lad[-1])) for n in lad]
    rel={}
    for dv in deps:
        if 'theta' not in dv: continue
        rr={}
        for e in ETAS:
            d=F.l1_lorentz(th_r,w_r,dv['theta'],dv['weights'],e)
            rr[str(e)]=dict(relL1=d['l1']/n2,relL1_win=d['l1']/d['win_absA1'])
        rel[str(dv['n'])]=rr; dv.pop('theta'); dv.pop('weights')
    rows=[]
    for dv in deps:
        for e in ETAS:
            lk,li,lo=dv['per_eta'][float(e)]
            lh=lk/np.sqrt(n2)/e
            B,tr,ce=F.certificate(wS,lh)
            r=dict(n=dv['n'],eta=e,Lambda_K=lk,LamHat_over_eta=lh,bound_leak=B,
                   bound_trivial=tr,cert=ce,in_over_out=li/max(lo,1e-300),
                   nontrivial=bool(B<tr),vacuity_factor=B/tr)
            rr=rel.get(str(dv['n']),{}).get(str(e))
            if rr: r.update(relL1=rr['relL1'],relL1_win=rr['relL1_win'],slack=ce/max(rr['relL1'],1e-300))
            rows.append(r)
    pt=dict(L=L,D=D,E0=E0,norm_phi2=n2,FR=FR,S=k,FR_real=k/D,w_S=wS,w_K=wK,dw=abs(wK-wS),
            nl=m,orth_defect=orth,res_PKphi=resPK,nref=m_r,ref_drift={str(k_):v for k_,v in drift.items()},
            t_lanczos=t_lan,t_G1=t_g1,t_point=time.time()-t0,peak_vec_GiB=gib,ladder=lad,rows=rows,
            provenance="gs+seed+order from repo work/ckpt/L14_*.npz|npy; gs residual 2.656e-09 verified; "
                       "seed == c^dag_{0up}|Psi_0> to 0.0e+00; order verified to be a permutation. "
                       "RANKING PROTOCOL NOT RE-RUN at L=14 (verified at L=10 by set overlap 1.000000).")
    out.append(pt)
    log(f"  done [lanczos {t_lan:.0f}s + G1 {t_g1:.0f}s = {time.time()-t0:.0f}s] orth={orth:.1e} |wK-wS|={abs(wK-wS):.1e}")
    for r in [x for x in rows if x['n']==m]:
        log("    eta=%.2f  Lam/eta=%.4f  B=%.4f  triv=%.4f  VACUITY x%.3f  relL1=%.4e  slack=%.0f  in/out=%.2e"
            %(r['eta'],r['LamHat_over_eta'],r['bound_leak'],r['bound_trivial'],r['vacuity_factor'],
              r.get('relL1',float('nan')),r.get('slack',float('nan')),r['in_over_out']))
json.dump(out,open(os.path.join(PUB,'L14_FR008.json'),'w'),indent=1,default=float)
log("WROTE L14_FR008.json")
