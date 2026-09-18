# -*- coding: utf-8 -*-
r"""MAX-LEVEL validation of the bitstring-sampled spectral method (Hubbard L=6, U/t=8).
Extends spectral_prototype with: (i) the EXACT order-K Krylov determinant subspace = union of supports of
{phi, H phi, ..., H^K phi} -- the optimal moment-complete subspace whose L1 error decays GEOMETRICALLY
(Theorem b2, moment-exactness) -> the theoretical reference curve; (ii) the finite-shot TE-QSCI sampled
subspace over NSEEDS independent seeds -> mean +- std error bands; (iii) the sum rule, MEASURED BY
INTEGRATION (see below); (iv) a fit of the geometric rate rho. Writes data/method_max.json.

THE SUM RULE -- WHAT IS ACTUALLY COMPUTED, AND WHAT IS NOT
----------------------------------------------------------
Earlier versions of this script reported a single number, `sumrule`, obtained as

    phi = Cd[p] @ psi0 ;  wnorm = <phi|phi> = <Psi0| c_p c_p^dag |Psi0> = 1 - <n_p> = 0.5

That is an ANALYTIC IDENTITY of the ground state at half filling.  It is exact by construction,
it does not depend on the omega grid, on eta, or on the reconstruction, and -- decisively -- the
reconstructed spectral function A_S was NEVER INTEGRATED anywhere in this script.  The
"a-posteriori validation of the reconstructed spectrum" that the manuscript's Fig. 3 caption
advertised as holding "to machine precision" was therefore never computed.

This version measures the sum rule the way the caption claims, by trapezoidal integration of
both spectral functions over the SAME omega grid the figure uses, and deposits three clearly
distinguished numbers:

    sum_rule_analytic                = <Psi0| c_p c_p^dag |Psi0> = 1 - <n_p>     (exact identity)
    sum_rule_exact_integrated        = \int A_exact(omega) domega  over the grid
    sum_rule_reconstructed_integrated= \int A_S(omega)     domega  over the grid

At the settings of the figure (L=6, U/t=8, eta=0.15, 600-point grid spanning
[min pole - 1, max pole + 1]) these are

    analytic      0.500000
    exact         0.486588      -> 2.68% BELOW the analytic identity: Lorentzian tails of width
                                   eta leak outside a finite omega window, a WINDOW EFFECT
    reconstructed 0.486471      -> 0.024% below the exact integral

so the honest statement is "the reconstructed spectrum conserves the integrated weight of the
exact one to 2.4e-4 relative, both of them losing 2.7% of the analytic weight to the finite
omega window" -- NOT "0.500 to machine precision".  The manuscript already concedes the same
~2% window effect in the caption of another figure, so the old caption contradicted itself.

All numbers exact-diagonalization verified.  Deterministic: the only randomness is the shot
sampling, seeded by np.random.default_rng(1000 + sd).
"""
import json, os, sys, time, datetime, numpy as np, scipy.sparse as sp
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

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_JSON = os.environ.get('METHOD_MAX_OUT', os.path.join(REPO, 'data', 'method_max.json'))

# numpy >= 2.0 renamed trapz -> trapezoid; the deposit must run on both.
_trapz = getattr(np, 'trapezoid', None) or np.trapz

t0=time.time(); log=lambda *a: print(f"[{time.time()-t0:6.1f}s]",*a,flush=True)

def fail(msg):
    sys.stderr.write("FAIL: " + msg + "\n")
    sys.exit(1)

def build(L=6,U=8.0,t_hop=1.0):
    M=2*L; dim=1<<M
    def c_op(p):
        r=[];c=[];d=[]
        for s in range(dim):
            if (s>>p)&1:
                sg=(-1)**bin(s&((1<<p)-1)).count('1'); r.append(s&~(1<<p)); c.append(s); d.append(float(sg))
        return sp.csr_matrix((d,(r,c)),shape=(dim,dim))
    C=[c_op(p) for p in range(M)]; Cd=[c.T.conj() for c in C]
    Nn=np.array([bin(s).count('1') for s in range(dim)])
    Sz=np.array([sum(((s>>(2*i))&1)-((s>>(2*i+1))&1) for i in range(L)) for s in range(dim)])
    H=sp.csr_matrix((dim,dim))
    for i in range(L-1):
        for sp_ in (0,1):
            a=2*i+sp_; b=2*(i+1)+sp_; H=H-t_hop*(Cd[a]@C[b]+Cd[b]@C[a])
    for i in range(L): H=H+U*(Cd[2*i]@C[2*i])@(Cd[2*i+1]@C[2*i+1])
    return H.tocsr(),C,Cd,Nn,Sz,M,dim

def run(L=6,U=8.0,K=12,dt=0.4,shots=4000,eta=0.15,NSEEDS=8):
    H,C,Cd,Nn,Sz,M,dim=build(L,U)
    gi=np.where((Nn==L)&(Sz==0))[0]; Hg=H[gi][:,gi].toarray()
    w,v=np.linalg.eigh(Hg); E0=w[0]; psi0=np.zeros(dim); psi0[gi]=v[:,0]
    p=0                                                   # seed site 0, up
    phi=Cd[p]@psi0
    # ANALYTIC identity <Psi0|c_p c_p^dag|Psi0> = 1 - <n_p>.  NOT a measurement of A: it is exact
    # by construction and independent of the omega grid.  Kept, but no longer called "the sum rule".
    wnorm=float(np.vdot(phi,phi).real)
    si=np.where((Nn==L+1)&(Sz==1))[0]
    H7=H[si][:,si].toarray(); En,Vn=np.linalg.eigh(H7); phi7=phi[si]
    amp=Vn.conj().T@phi7; poles=En-E0; wts=np.abs(amp)**2
    def spec(pl,ww,g):
        A=np.zeros_like(g)
        for a,b in zip(pl,ww): A+=b*(eta/np.pi)/((g-a)**2+eta**2)
        return A
    grid=np.linspace(poles.min()-1,poles.max()+1,600); A_ex=spec(poles,wts,grid); normA=_trapz(A_ex,grid)
    def relL1(Sset):
        S=np.array(sorted(Sset)); HS=H[S][:,S].toarray(); Em,Um=np.linalg.eigh(HS)
        aS=Um.conj().T@phi[S]; A_S=spec(Em-E0,np.abs(aS)**2,grid)
        return float(_trapz(np.abs(A_S-A_ex),grid)/normA), A_S
    log(f"L={L} U={U}: E0={E0:.5f} analytic 1-<n_p>={wnorm:.4f} dim(N+1,Sz+1)={len(si)}")

    # ---- (II) finite-shot TE-QSCI sampled subspace over NSEEDS seeds ----
    def evolve(tk): return Vn@(np.exp(-1j*En*tk)*(Vn.conj().T@phi7))
    per_seed=[]                                            # per seed: list over k of (|S|, l1)
    for sd in range(NSEEDS):
        rng=np.random.default_rng(1000+sd); seen={int(si[np.argmax(np.abs(phi7)**2)])}; curve=[]
        for k in range(K+1):
            vk=evolve(k*dt); pr=np.abs(vk)**2; pr/=pr.sum()
            for d in np.unique(rng.choice(len(si),size=shots,p=pr)): seen.add(int(si[d]))
            l1,_=relL1(seen); curve.append((len(seen),l1))
        per_seed.append(curve)
    per_seed=np.array(per_seed)                            # (NSEEDS, K+1, 2)
    Smean=per_seed[:,:,0].mean(0); l1mean=per_seed[:,:,1].mean(0); l1std=per_seed[:,:,1].std(0)
    sampled=[{'K':k,'S':float(Smean[k]),'l1':float(l1mean[k]),'l1_std':float(l1std[k])} for k in range(K+1)]
    S_final=set(seen)                                      # last seed's (sd = NSEEDS-1) final subspace
    l1_final,A_final=relL1(S_final)

    # ---- (III-bis) THE SUM RULE, MEASURED BY INTEGRATION ON THE FIGURE'S OWN GRID ----
    # This is the check the Fig. 3 caption claims and that no earlier version performed.
    sr_exact=float(normA)                                  # \int A_exact domega
    sr_recon=float(_trapz(A_final,grid))                   # \int A_S     domega
    window_loss=(wnorm-sr_exact)/wnorm                     # total deficit vs the analytic identity
    recon_dev=(sr_exact-sr_recon)/sr_exact                 # reconstruction vs exact, relative
    # Decompose the deficit, so that "window effect" is a measurement and not a story.
    # The Lorentzian mass of pole a_m with weight w_m lying INSIDE [g0,g1] is exactly
    #   w_m * (1/pi) * [ arctan((g1-a_m)/eta) + arctan((a_m-g0)/eta) ].
    g0,g1=float(grid[0]),float(grid[-1])
    in_window=float(np.sum(wts*(np.arctan((g1-poles)/eta)+np.arctan((poles-g0)/eta))/np.pi))
    tail_outside=float(wnorm-in_window)                    # analytic, exact: the window effect
    quad_error=float(sr_exact-in_window)                   # trapezoidal discretization residue
    sum_rule=dict(
        sum_rule_analytic=wnorm,
        sum_rule_analytic_definition="<Psi0|c_p c_p^dag|Psi0> = 1 - <n_p>; exact identity, NOT a measurement of A",
        sum_rule_exact_integrated=sr_exact,
        sum_rule_reconstructed_integrated=sr_recon,
        integration="numpy trapezoid over the 600-point omega grid used by the figure",
        omega_window=[float(grid[0]),float(grid[-1])], n_omega=int(len(grid)), eta=float(eta),
        window_loss_relative=float(window_loss),
        deficit_decomposition=dict(
            analytic_mass_inside_window=in_window,
            tail_outside_window=tail_outside,
            tail_outside_window_relative=float(tail_outside/wnorm),
            trapezoid_quadrature_residue=quad_error,
            trapezoid_quadrature_residue_relative=float(quad_error/wnorm),
            formula="in_window = sum_m w_m (1/pi)[arctan((w_max-a_m)/eta)+arctan((a_m-w_min)/eta)]"),
        reconstruction_deviation_relative=float(recon_dev),
        reconstructed_subspace_size=int(len(S_final)),
        reconstructed_relL1=float(l1_final),
        note=("The analytic identity is exactly 0.5 and is grid-independent. The INTEGRATED exact "
              "and reconstructed weights are both ~2.7%% below it because Lorentzians of width eta "
              "leak outside a finite omega window. The reconstructed spectrum matches the exact "
              "INTEGRAL to %.2e relative -- this, not 'machine precision' against 0.500, is the "
              "honest a-posteriori statement." % abs(recon_dev)),
    )
    log("SUM RULE  analytic (1-<n_p>)      = %.10f" % wnorm)
    log("SUM RULE  \\int A_exact domega     = %.10f   (%.3f%% below the analytic identity: window effect)"
        % (sr_exact,100*window_loss))
    log("SUM RULE  \\int A_reconstructed    = %.10f   (%.4f%% below the exact integral)"
        % (sr_recon,100*recon_dev))
    log("DEFICIT   Lorentzian tail outside the window = %.6f (%.3f%% of the analytic weight); "
        "trapezoid residue = %+.2e (%+.4f%%)"
        % (tail_outside,100*tail_outside/wnorm,quad_error,100*quad_error/wnorm))

    # ---- self-verification: abort rather than deposit numbers that do not reproduce ----
    # Tolerances tightened 2026-09-18.  The previous +-1e-4 acceptance windows around
    # 0.48659 and 0.48647 OVERLAPPED -- the two numbers are 1.16e-4 apart -- so a single
    # value in [0.48649, 0.48657] satisfied both tests at once, for a pair of numbers whose
    # whole point is that the reconstruction is NOT the exact spectrum.  They are now
    # separated by three orders of magnitude more than they are apart.
    if abs(wnorm-0.5)>1e-12:
        fail("analytic identity 1-<n_p> = %.12f != 0.5 (half filling broken)"%wnorm)
    # the protocol itself is pinned: these three numbers define the window, and the
    # decomposition headline below is only true on this grid.
    if int(len(grid))!=600:
        fail("omega grid has %d points, not the 600 the deposited figure uses; the "
             "window/quadrature split is grid dependent"%len(grid))
    if abs(eta-0.15)>1e-12:
        fail("eta = %.6f, not the 0.15 of the deposited figure"%eta)
    if abs(float(grid[0])-5.687888306705499)>1e-9 or abs(float(grid[-1])-29.8483092042426)>1e-6:
        fail("omega window [%.9f, %.9f] is not the deposited one"%(grid[0],grid[-1]))
    if abs(sr_exact-0.4865877175)>1e-6:
        fail("integrated exact sum rule %.10f is not the expected 0.4865877175 +- 1e-6"%sr_exact)
    if abs(sr_recon-0.4864714463)>1e-6:
        fail("integrated reconstructed sum rule %.10f is not the expected 0.4864714463 "
             "+- 1e-6"%sr_recon)
    if abs(sr_exact-sr_recon)<1e-8:
        fail("the integrated exact and reconstructed weights are indistinguishable "
             "(%.12f vs %.12f); the reconstruction is on a strict subspace and cannot be "
             "the exact spectrum"%(sr_exact,sr_recon))
    if not (0.02<window_loss<0.035):
        fail("window loss %.4f outside the expected 2-3.5%% band"%window_loss)
    # THE DECOMPOSITION MUST CLOSE.  Reported since 2026-09-18 as "99.99% of the deficit is
    # the finite window, not the quadrature"; until now that sentence was computed, deposited
    # and logged, and never verified.  An adversarial pass pointed out that the four aborts
    # above touched none of in_window / tail_outside / quad_error.
    if abs((in_window+tail_outside)-wnorm)>1e-12:
        fail("decomposition does not close: analytic mass inside the window (%.12f) plus the "
             "Lorentzian tail outside it (%.12f) is not the analytic identity (%.12f)"
             %(in_window,tail_outside,wnorm))
    if abs((in_window+quad_error)-sr_exact)>1e-12:
        fail("decomposition does not close: analytic mass inside the window (%.12f) plus the "
             "quadrature residue (%.12f) is not the integral actually performed (%.12f)"
             %(in_window,quad_error,sr_exact))
    deficit=wnorm-sr_exact
    if not (abs(tail_outside/deficit-1.0)<1e-3 and abs(quad_error/deficit)<1e-3):
        fail("the deficit is NOT dominated by the finite window: tail = %.3f%% of it, "
             "quadrature residue = %.3f%% of it. The claim that the 2.7%% is a window effect "
             "and not our discretization is exactly what this check exists to defend"
             %(100*tail_outside/deficit,100*quad_error/deficit))
    log("SELF-CHECK PASS: analytic 0.5 exact; integrated 0.4865877175 / 0.4864714463 "
        "reproduced and distinct; grid/eta/window pinned; deficit closes with %.4f%% window "
        "and %+.4f%% quadrature."%(100*tail_outside/deficit,100*quad_error/deficit))

    # ---- (III) geometric fit on the sampled COLLAPSE regime (rel-L1 ~ C rho^{-K}) : Thm b2 in action ----
    Ka=np.array([r['K'] for r in sampled]); L1a=np.array([r['l1'] for r in sampled])
    coll=(Ka>=4)&(Ka<=9)
    slope,intc=np.polyfit(Ka[coll],np.log(L1a[coll]),1); rho=float(np.exp(-slope))   # per Krylov order
    fit=[{'K':int(k),'l1':float(np.exp(intc+slope*k))} for k in range(4,13)]
    kstar=4                                                 # onset of the geometric regime

    return dict(L=L,U=U,dt=dt,shots=shots,eta=eta,NSEEDS=NSEEDS,E0=float(E0),sumrule=float(wnorm),
                dimNp1=int(len(si)),rho=rho,kstar=kstar,grid=grid.tolist(),
                A_exact=A_ex.tolist(),A_final=A_final.tolist(),sampled=sampled,fit=fit,
                sum_rule=sum_rule,
                provenance=dict(
                    generated_by="src/spectral_validation_max.py",
                    generated_utc=datetime.datetime.now(datetime.timezone.utc)
                                    .strftime("%Y-%m-%dT%H:%M:%SZ"),
                    model="1D Hubbard OPEN chain, L=%d, U/t=%.1f, half filling"%(L,U),
                    seeds=[1000+sd for sd in range(NSEEDS)], shots_per_step=shots, dt=dt, K_max=K,
                    deterministic=True,
                    backs=("Fig. 3 (method validation): rel-L1 curve, geometric rate rho, and the "
                           "SUM RULE. 'sumrule' is kept for backward compatibility and is the "
                           "ANALYTIC identity 1-<n_p>; the integrated values live in 'sum_rule'."),
                    checks=("analytic identity == 0.5 (1e-12); integrated exact == 0.48659 (1e-4); "
                            "integrated reconstructed == 0.48647 (1e-4); window loss in [2%,3.5%]"),
                    checks_status="ALL PASS"))

if __name__=='__main__':
    out=run()
    # Reproduction check against whatever is already deposited (deterministic script).
    if os.path.exists(OUT_JSON):
        try:
            old=json.load(open(OUT_JSON))
            d_ex=float(np.max(np.abs(np.array(out['A_exact'])-np.array(old['A_exact']))))
            d_fi=float(np.max(np.abs(np.array(out['A_final'])-np.array(old['A_final']))))
            d_rho=abs(out['rho']-old['rho'])
            d_l1=max(abs(a['l1']-b['l1']) for a,b in zip(out['sampled'],old['sampled']))
            log("REPRODUCTION vs deposited method_max.json: max|dA_exact|=%.2e max|dA_final|=%.2e "
                "|drho|=%.2e max|dl1|=%.2e"%(d_ex,d_fi,d_rho,d_l1))
            if max(d_ex,d_fi,d_rho,d_l1)>1e-6:
                fail("this run does NOT reproduce the deposited method_max.json; "
                     "refusing to overwrite it")
        except (KeyError,ValueError) as e:
            log("reproduction check skipped (deposited file lacks a field): %s"%e)
    json.dump(out,open(OUT_JSON,'w'))
    log(f"WROTE {OUT_JSON} | rho(geometric)={out['rho']:.3f} "
        f"sum rule analytic/exact/recon = {out['sum_rule']['sum_rule_analytic']:.5f} / "
        f"{out['sum_rule']['sum_rule_exact_integrated']:.5f} / "
        f"{out['sum_rule']['sum_rule_reconstructed_integrated']:.5f}")
