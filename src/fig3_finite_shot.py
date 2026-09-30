# -*- coding: utf-8 -*-
r"""Finite-shot run of the literal Fig. 3(a,b) configuration (fig:akwsampled, "The shot budget of
(a)--(b), measured").

L=8, U=8, eta=0.18t, K=16 (17 Born slices), dt=0.5, all 8 momenta x 2 branches = 16 channels.
Selection = determinants actually OBSERVED in T shots (multinomial per slice), not argsort.
A_S by Haydock on the restricted operator (verified identical to dense Rayleigh-Ritz to
printed precision in fig5_cost.py, a 2026-09-19 probe that is not deposited).

PROVENANCE OF THIS FILE
-----------------------
This is the generator of the finite-shot numbers of the Fig. 3 caption.  It was written on
2026-09-19 as `fig5_literal.py` (T = 2.6e6) and `fig5_literal_T52.py` (T = 5.2e6, the same file
after a two-line sed: T_SHOTS and the output name), run in a session scratch folder that was
later purged, and never deposited.  On 2026-09-28 it was recovered byte for byte from the session
transcript (subagent agent-a07c2e0927cf28823 of workflow wf_bb98610c-331, transcript lines
240/241 = the heredoc that wrote the file, 625/626 = the sed, 959/960 = the post-processing
step); sha256 of the recovered file (LF endings, 3515 bytes):
e2b3793a2267953f4fab6adfb293e891033bc492c918069c2d27117a237873b5.

The ONLY changes against the recovered file are:
  1. T (was the constant T_SHOTS) and the frequency window (was wg=np.linspace(-9,9,600)) are
     command-line options;
  2. the import path is repository-relative and the output goes to data/ (was: next to the
     script); the numpy trapz/trapezoid bridge used elsewhere in the repository is added;
  3. a provenance block, the post-processing step, and the reproduction check are added.
Every algorithmic line -- the Hamiltonians, the ARPACK start vector, the dense eigh of the two
sectors, the multinomial draws of T//(K+1) shots per slice, the per-(seed, momentum) generator
default_rng(seed*1000+n), the Haydock depth NL=260, the Lorentzian reference, the trapezoid
rel-L1 -- is the recovered one, unchanged.

POST-PROCESSING (the 2026-09-19 compare step, now part of the generator)
------------------------------------------------------------------------
Over the 32 (momentum, seed) rows (each row sums the addition and removal branches; the log
line below keeps the 2026-09-19 wording "(channel,seed)" because the reproduction check compares
it character by character): mean rel-L1, sample s.d. (ddof=1), max, mean frac_add, mean
frac_rem; the ratio of the mean to the mean of Fig. 3(b) read from data/sampled_akw_L8.json; and
the UNDERSTATEMENT of the transferred penalty, defined here explicitly because no 2026-09-19
script printed it:

    understatement = 1 - (transferred penalty) / (measured ratio at T = 2.6e6),

where the transferred penalty is Table S6's L=8 factor at |S| = 2180, read from
data/sampled_honest.json (`penalty_factor_topm` of the matched_subspace row with S_matched = 2180,
2.052136..., printed in Sec. S6 as 2.05) and the measured ratio is the T = 2.6e6 mean divided by
the Fig. 3(b) mean.  (On the former window this was 1 - 2.05/2.80 = 27%.)
The ratios are written only when the run is scored on the window of Fig. 3(b) (the window and
grid of data/sampled_akw_L8.json); on any other window they would divide two different integrals
and are stored as null.  The understatement is written only for T = 2.6e6, the budget whose
observed fraction (0.839) is the one closest to the 0.850 of Fig. 3(a,b); for other T it is null.

PROVENANCE HASH
---------------
provenance.generator_sha256 is the sha256 of THIS file with line endings normalised to LF, so
it does not depend on git's core.autocrlf; src/verify.py checks it against the file.

REPRODUCTION CHECK
------------------
Run on the ORIGINAL window [-9t, 9t] with 600 points (--wmin -9 --wmax 9 --npts 600) this file
must reproduce the 2026-09-19 output exactly: the per-channel relL1 printed then (all 32 lines at
T = 2.6e6, the 28 lines the transcript still shows at T = 5.2e6), the |S| sizes, the MEAN/MAX/frac
line and the compare step's s.d. and frac_rem.  `--repro OLD.json` compares such an old-window
output against those printed values and records the result in the new output JSON under
`reproduction_check`; the run fails (exit 1) if anything differs.

Run (from anywhere; about 5-6 min per run single-threaded on the reference machine, longer when
several runs share the cores):
    python src/fig3_finite_shot.py --T 2600000 --wmin -9 --wmax 9 --npts 600 --out SCRATCH/old_T2.6e6.json
    python src/fig3_finite_shot.py --T 2600000 --repro SCRATCH/old_T2.6e6.json   # -> data/fig3_finite_shot_T2.6e6.json
    python src/fig3_finite_shot.py --T 5200000 --wmin -9 --wmax 9 --npts 600 --out SCRATCH/old_T5.2e6.json
    python src/fig3_finite_shot.py --T 5200000 --repro SCRATCH/old_T5.2e6.json   # -> data/fig3_finite_shot_T5.2e6.json
    python src/fig3_finite_shot.py --summary      # prints every caption number from the two data files
Defaults: the window of Fig. 3(a,b), omega-E0 in [-9t, 17t] on 866 points, as in src/sampled_akw.py.
"""
import os, sys, time, json, argparse, hashlib, platform, datetime, numpy as np
sys.dont_write_bytecode = True
_SRC = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_SRC)
sys.path.insert(0, _SRC)
# --- NUMPY_TRAPEZOID_BRIDGE (as in src/sampled_akw.py) ---------------------
if not hasattr(np, "trapezoid"):
    np.trapezoid = np.trapz          # numpy < 2.0
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid          # numpy >= 2.0
# ---------------------------------------------------------------------------

K_CHANNELS = 16                      # 8 momenta x 2 branches

# ---- What the 2026-09-19 run printed on the window [-9t, 9t] (600 points) ----------------------
# Copied from the transcript's tool_result records; used ONLY by --repro.
# (k/pi, seed): (relL1 as printed '%.4e', |S|add, |S|rem)
PUBLISHED_OLD_WINDOW = {
    2_600_000: dict(
        rows={
            (0.000, 7001): ('6.2585e-03', 3292, 3084), (0.000, 7002): ('7.1255e-03', 3285, 3099),
            (0.000, 7003): ('5.9450e-03', 3293, 3101), (0.000, 7004): ('6.0166e-03', 3287, 3101),
            (0.250, 7001): ('6.4591e-03', 3406, 3220), (0.250, 7002): ('6.3171e-03', 3407, 3204),
            (0.250, 7003): ('6.1366e-03', 3410, 3228), (0.250, 7004): ('6.9682e-03', 3408, 3216),
            (0.500, 7001): ('5.9357e-03', 3335, 3341), (0.500, 7002): ('5.9557e-03', 3336, 3345),
            (0.500, 7003): ('5.4468e-03', 3367, 3350), (0.500, 7004): ('5.9866e-03', 3353, 3351),
            (0.750, 7001): ('5.9041e-03', 3199, 3407), (0.750, 7002): ('5.9451e-03', 3201, 3417),
            (0.750, 7003): ('5.5947e-03', 3229, 3387), (0.750, 7004): ('5.7134e-03', 3204, 3430),
            (1.000, 7001): ('6.9140e-03', 3074, 3290), (1.000, 7002): ('6.1987e-03', 3098, 3272),
            (1.000, 7003): ('6.2126e-03', 3095, 3292), (1.000, 7004): ('6.2063e-03', 3088, 3309),
            (1.250, 7001): ('5.6452e-03', 3196, 3447), (1.250, 7002): ('5.5441e-03', 3217, 3404),
            (1.250, 7003): ('5.6158e-03', 3198, 3408), (1.250, 7004): ('5.8918e-03', 3209, 3405),
            (1.500, 7001): ('5.6615e-03', 3369, 3330), (1.500, 7002): ('5.8148e-03', 3342, 3356),
            (1.500, 7003): ('5.7657e-03', 3374, 3329), (1.500, 7004): ('5.6341e-03', 3357, 3349),
            (1.750, 7001): ('7.6383e-03', 3391, 3179), (1.750, 7002): ('6.3715e-03', 3403, 3211),
            (1.750, 7003): ('6.5130e-03', 3432, 3214), (1.750, 7004): ('6.5659e-03', 3396, 3222)},
        summary_line="MEAN over 32 (channel,seed) = 6.1219e-03   MAX = 7.6383e-03   mean frac_add=0.8391",
        sd='5.05e-04', frac_rem='0.8394', ratio_mean='2.80', ratio_max='2.33'),
    5_200_000: dict(
        rows={  # the transcript still shows 28 of the 32 lines
            (0.000, 7001): ('2.5625e-03', 3455, 3289), (0.000, 7002): ('2.5691e-03', 3460, 3302),
            (0.000, 7003): ('2.3772e-03', 3443, 3298), (0.000, 7004): ('2.4815e-03', 3444, 3299),
            (0.250, 7001): ('2.2024e-03', 3561, 3414), (0.250, 7002): ('1.8860e-03', 3568, 3423),
            (0.250, 7003): ('2.2439e-03', 3535, 3415), (0.250, 7004): ('2.4119e-03', 3551, 3405),
            (0.500, 7001): ('1.6965e-03', 3512, 3524), (0.500, 7002): ('1.9347e-03', 3512, 3505),
            (0.500, 7003): ('1.7912e-03', 3496, 3494), (0.750, 7004): ('1.9961e-03', 3394, 3552),
            (1.000, 7001): ('2.1901e-03', 3304, 3453), (1.000, 7002): ('2.4135e-03', 3301, 3441),
            (1.000, 7003): ('2.3540e-03', 3290, 3450), (1.000, 7004): ('2.4260e-03', 3300, 3435),
            (1.250, 7001): ('1.9603e-03', 3404, 3565), (1.250, 7002): ('2.1146e-03', 3392, 3560),
            (1.250, 7003): ('1.6951e-03', 3413, 3563), (1.250, 7004): ('1.9197e-03', 3406, 3544),
            (1.500, 7001): ('1.9777e-03', 3500, 3491), (1.500, 7002): ('1.7864e-03', 3509, 3498),
            (1.500, 7003): ('1.9311e-03', 3494, 3497), (1.500, 7004): ('1.8403e-03', 3516, 3511),
            (1.750, 7001): ('2.1875e-03', 3553, 3421), (1.750, 7002): ('2.2885e-03', 3563, 3416),
            (1.750, 7003): ('2.3975e-03', 3557, 3413), (1.750, 7004): ('2.0308e-03', 3550, 3426)},
        summary_line="MEAN over 32 (channel,seed) = 2.1134e-03   MAX = 2.5691e-03   mean frac_add=0.8822",
        sd='2.58e-04', frac_rem='0.8833', ratio_mean='0.97', ratio_max='0.78'),
}
# The Fig. 3(b) mean and max the 2026-09-19 compare step divided by: data/sampled_akw_L8.json
# as committed then (window [-9t, 9t]); used ONLY by --repro, to reproduce its printed ratios.
OLD_WINDOW_PANEL_B = dict(mean_relL1=0.0021858473952153625, max_relL1=0.0032826276455945177)


def tag(T):
    """2600000 -> '2.6e6' (the data-file suffix)."""
    return "%ge6" % (T / 1e6)


def default_out(T):
    return os.path.join(_REPO, "data", "fig3_finite_shot_T%s.json" % tag(T))


def sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def sha256_lf(path):
    """sha256 of a text file with CRLF normalised to LF (independent of git core.autocrlf)."""
    with open(path, "rb") as f:
        return hashlib.sha256(f.read().replace(b"\r\n", b"\n")).hexdigest()


def committed_references():
    """Fig. 3(b) mean/max (data/sampled_akw_L8.json) and Table S6's transferred L=8 factor at
    |S| = 2180 (data/sampled_honest.json)."""
    p_b = os.path.join(_REPO, "data", "sampled_akw_L8.json")
    p_h = os.path.join(_REPO, "data", "sampled_honest.json")
    b = json.load(open(p_b))
    h = json.load(open(p_h))
    tr = None
    for r in h["results"]:
        for m in (r.get("matched_subspace") or []):
            if isinstance(m, dict) and m.get("S_matched") == 2180:
                tr = (r["L"], m["penalty_factor_topm"], m["penalty_factor_union"])
    assert tr is not None and tr[0] == 8, "Table S6 L=8 |S|=2180 row not found in sampled_honest.json"
    assert tr[1] == tr[2], "top-m and union penalty factors differ; the definition would be ambiguous"
    return dict(panel_b_file="data/sampled_akw_L8.json", panel_b_sha256=sha256(p_b),
                panel_b_window=b.get("window"), panel_b_npts=len(b["wg"]), panel_b_frac=b["frac"],
                panel_b_mean_relL1=b["mean_relL1"], panel_b_max_relL1=b["max_relL1"],
                transferred_file="data/sampled_honest.json", transferred_sha256=sha256(p_h),
                transferred_penalty_L8_S2180=tr[1])


def same_window_as_panel_b(refs, window, npts):
    """True only if the run is scored on the window and grid of Fig. 3(b)."""
    return (refs is not None and window is not None
            and refs.get('panel_b_window') == [float(window[0]), float(window[1])]
            and refs.get('panel_b_npts') == int(npts))


def postprocess(rows, refs=None, window=None, npts=None, T=None):
    """The 2026-09-19 compare step (mean, ddof=1 s.d., max, fractions, ratios) plus the derived
    figures the caption quotes.  Ratios only on the window of Fig. 3(b); the understatement only
    at T = 2.6e6 (null otherwise, see the module docstring)."""
    a = np.array([r['relL1'] for r in rows])
    fa = np.array([r['frac_add'] for r in rows]); fr = np.array([r['frac_rem'] for r in rows])
    ks = sorted(set(r['k_over_pi'] for r in rows))
    out = dict(n_pairs=int(a.size), mean_relL1=float(a.mean()), sd_relL1_ddof1=float(a.std(ddof=1)),
               max_relL1=float(a.max()), min_relL1=float(a.min()),
               mean_frac_add=float(fa.mean()), mean_frac_rem=float(fr.mean()),
               per_k_mean_relL1={"%.3f" % k: float(np.mean([r['relL1'] for r in rows if r['k_over_pi'] == k]))
                                 for k in ks})
    if refs is not None:
        if same_window_as_panel_b(refs, window, npts):
            out['ratio_mean_to_panel_b'] = out['mean_relL1'] / refs['panel_b_mean_relL1']
            out['ratio_max_to_panel_b'] = out['max_relL1'] / refs['panel_b_max_relL1']
        else:
            out['ratio_mean_to_panel_b'] = out['ratio_max_to_panel_b'] = None
            out['ratio_note'] = ("not computed: this run is scored on %s x %s, Fig. 3(b) on %s x %s"
                                 % (window, npts, refs.get('panel_b_window'), refs.get('panel_b_npts')))
        if T == 2_600_000 and out['ratio_mean_to_panel_b'] is not None:
            out['understatement_of_transferred_penalty'] = \
                1.0 - refs['transferred_penalty_L8_S2180'] / out['ratio_mean_to_panel_b']
        else:
            out['understatement_of_transferred_penalty'] = None
    return out


def reproduction_check(old_path, T):
    """Compare an old-window ([-9t,9t], 600 points) output of THIS file with what the 2026-09-19
    run printed.  Returns a dict; 'exact' is True only if every printed value is identical."""
    d = json.load(open(old_path))
    pub = PUBLISHED_OLD_WINDOW[T]
    cfg = d['config']
    assert cfg['T'] == T and cfg['window'] == [-9.0, 9.0] and cfg['npts'] == 600, \
        "--repro needs an output of this file at the same T on the original window [-9,9] x 600"
    got = {("%.3f" % r['k_over_pi'], r['seed']): r for r in d['rows']}
    mism = []
    for (k, sd), (rel, sa, sr) in pub['rows'].items():
        r = got[("%.3f" % k, sd)]
        if ("%.4e" % r['relL1'], r['S_add'], r['S_rem']) != (rel, sa, sr):
            mism.append(dict(k_over_pi=k, seed=sd, printed=[rel, sa, sr],
                             reproduced=["%.4e" % r['relL1'], r['S_add'], r['S_rem']]))
    s = postprocess(d['rows'])
    line = "MEAN over %d (channel,seed) = %.4e   MAX = %.4e   mean frac_add=%.4f" % (
        s['n_pairs'], s['mean_relL1'], s['max_relL1'], s['mean_frac_add'])
    checks = {
        'summary_line': [pub['summary_line'], line],
        'sd_ddof1': [pub['sd'], "%.2e" % s['sd_relL1_ddof1']],
        'mean_frac_rem': [pub['frac_rem'], "%.4f" % s['mean_frac_rem']],
        'ratio_mean_to_old_panel_b': [pub['ratio_mean'], "%.2f" % (s['mean_relL1'] / OLD_WINDOW_PANEL_B['mean_relL1'])],
        'ratio_max_to_old_panel_b': [pub['ratio_max'], "%.2f" % (s['max_relL1'] / OLD_WINDOW_PANEL_B['max_relL1'])],
    }
    exact = not mism and all(a == b for a, b in checks.values())
    return dict(what="this generator re-run on the ORIGINAL window omega-E0 in [-9t, 9t] (600 points) at the "
                     "same T, compared with the values the 2026-09-19 run printed (transcript tool_result "
                     "records); ratios against the Fig. 3(b) mean/max of that date on that window",
                window=[-9.0, 9.0], npts=600, T=T, old_window_output_sha256=sha256(old_path),
                old_window_output_generated_utc=d.get('provenance', {}).get('generated_utc'),
                n_rows_compared=len(pub['rows']), n_rows_mismatched=len(mism), row_mismatches=mism,
                printed_vs_reproduced=checks, old_window_statistics=s, exact=bool(exact))


def summary():
    refs = committed_references()
    D = {}
    for T in (2_600_000, 5_200_000):
        p = default_out(T)
        D[T] = json.load(open(p))
        w = D[T]['config']['window']
        s = postprocess(D[T]['rows'], refs, w, D[T]['config']['npts'], T)
        print("T=%.2e per channel (%.3g over the %d channels), window [%g,%g] x %d, n=%d" % (
            T, K_CHANNELS * T, K_CHANNELS, w[0], w[1], D[T]['config']['npts'], s['n_pairs']))
        print("   rel-L1 mean %.4e +- %.2e (ddof=1)  max %.4e  min %.4e" % (
            s['mean_relL1'], s['sd_relL1_ddof1'], s['max_relL1'], s['min_relL1']))
        print("   mean frac_add %.4f  mean frac_rem %.4f" % (s['mean_frac_add'], s['mean_frac_rem']))
        print("   ratio to Fig. 3(b): mean %.4f (%.2f)  max %.4f" % (
            s['ratio_mean_to_panel_b'], s['ratio_mean_to_panel_b'], s['ratio_max_to_panel_b']))
        print("   per k:", " ".join("%s:%.3e" % kv for kv in s['per_k_mean_relL1'].items()))
        rc = D[T].get('reproduction_check')
        print("   reproduction check on [-9t,9t]: %s" % (("EXACT (%d rows)" % rc['n_rows_compared'])
                                                        if rc and rc['exact'] else rc))
    c26 = D[2_600_000]['config']
    s26 = postprocess(D[2_600_000]['rows'], refs, c26['window'], c26['npts'], 2_600_000)
    print("Fig. 3(b): mean %.4e  max %.4e on window %s, fraction %.3f" % (
        refs['panel_b_mean_relL1'], refs['panel_b_max_relL1'], refs['panel_b_window'], refs['panel_b_frac']))
    print("transferred penalty (Table S6, L=8, |S|=2180): %.4f" % refs['transferred_penalty_L8_S2180'])
    print("understatement = 1 - %.4f / %.4f = %.4f  -> %.0f%%" % (
        refs['transferred_penalty_L8_S2180'], s26['ratio_mean_to_panel_b'],
        s26['understatement_of_transferred_penalty'], 100 * s26['understatement_of_transferred_penalty']))


ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("--T", type=int, help="shots per channel (2600000 or 5200000 in the paper)")
ap.add_argument("--wmin", type=float, default=-9.0, help="window lower edge, omega-E0 in units of t")
ap.add_argument("--wmax", type=float, default=17.0, help="window upper edge, omega-E0 in units of t")
ap.add_argument("--npts", type=int, default=866, help="grid points on the window")
ap.add_argument("--out", default=None, help="output JSON (default data/fig3_finite_shot_T<T>.json)")
ap.add_argument("--repro", default=None, help="old-window output of this file to check and record")
ap.add_argument("--summary", action="store_true", help="print the caption numbers from the two data files")
args = ap.parse_args()
if args.summary:
    summary(); sys.exit(0)
if args.T is None:
    ap.error("--T is required unless --summary")

# ===================== the recovered 2026-09-19 script (algorithm unchanged) =====================
import akw_lanczos as AK
from scipy.sparse.linalg import eigsh
T0=time.time()
def log(*a): print("[%7.1fs]"%(time.time()-T0), *a, flush=True)
L,U,t,eta,K,dt = 8,8.0,1.0,0.18,16,0.5
NL=260; SEEDS=[7001,7002,7003,7004]; T_SHOTS=args.T
nup=nd=L//2
H0,Du,Dd=AK.build_H_explicit(L,U,nup,nd,t)
E0,V0=eigsh(H0,k=1,which='SA',v0=np.random.default_rng(20260918).standard_normal(H0.shape[0]),tol=0)
E0=float(E0[0]); Psi=V0[:,0].reshape(Du,Dd); log("E0=%.9f"%E0)
Hadd,Dua,_=AK.build_H_explicit(L,U,nup+1,nd,t); Hrem,Dur,_=AK.build_H_explicit(L,U,nup-1,nd,t)
cdU=AK.cdag_map(L,nup); cUr=AK.cdag_map(L,nup-1)
wg=np.linspace(args.wmin,args.wmax,args.npts)
def spec_poles(poles,wts):
    A=np.zeros_like(wg)
    for a,b in zip(poles,wts): A+=b*(eta/np.pi)/((wg-a)**2+eta**2)
    return A
ta=time.time(); Ea,Ua=np.linalg.eigh(Hadd.toarray()); log("eigh add %.1fs"%(time.time()-ta))
ta=time.time(); Er,Ur=np.linalg.eigh(Hrem.toarray()); log("eigh rem %.1fs"%(time.time()-ta))
Hadd=Hadd.tocsr(); Hrem=Hrem.tocsr()
def channel(seed,Em,Um,Hc,sgn,rng):
    nS=Hc.shape[0]; coef=Um.conj().T@seed
    A_ex=spec_poles(sgn*(Em-E0),np.abs(coef)**2); nrm=np.trapz(np.abs(A_ex),wg)
    seen=np.zeros(nS,bool); per=T_SHOTS//(K+1)
    for kk in range(K+1):
        vk=Um@(np.exp(-1j*Em*kk*dt)*coef); q=np.abs(vk)**2
        seen |= (rng.multinomial(per,q/q.sum())>0)
    S=np.flatnonzero(seen)
    def rmv(xs):
        xf=np.zeros(nS,dtype=complex); xf[S]=xs; return Hc.dot(xf)[S]
    G=AK.haydock(lambda x: sgn*(rmv(x)-E0*x), seed[S], min(NL,S.size), wg+1j*eta)
    A_s=-G.imag/np.pi
    return A_ex,A_s,nrm,S.size,nS
res=[]
for n in range(L):
    k=2*np.pi*n/L; ph=np.exp(1j*k*np.arange(L))/np.sqrt(L)
    add=np.zeros((Dua,Dd),dtype=complex); rem=np.zeros((Dur,Dd),dtype=complex)
    for j in range(L):
        add+=ph[j]*(cdU[j]@Psi); rem+=np.conj(ph[j])*(cUr[j].T@Psi)
    sa=add.ravel(); sr=rem.ravel()
    for sd in SEEDS:
        rng=np.random.default_rng(sd*1000+n)
        Aa_e,Aa_s,_,ma,nSa=channel(sa,Ea,Ua,Hadd,+1.0,rng)
        Ar_e,Ar_s,_,mr,nSr=channel(sr,Er,Ur,Hrem,-1.0,rng)
        A_e=Aa_e+Ar_e; A_s=Aa_s+Ar_s
        r=float(np.trapz(np.abs(A_s-A_e),wg)/np.trapz(np.abs(A_e),wg))
        res.append(dict(k_over_pi=2*n/L,seed=sd,relL1=r,S_add=int(ma),S_rem=int(mr),
                        frac_add=ma/nSa,frac_rem=mr/nSr))
        log("k/pi=%.3f seed=%d  relL1=%.4e  |S|add=%d(%.4f) rem=%d(%.4f)"%(2*n/L,sd,r,ma,ma/nSa,mr,mr/nSr))
arr=np.array([d['relL1'] for d in res]); fr=np.array([d['frac_add'] for d in res])
log("MEAN over %d (channel,seed) = %.4e   MAX = %.4e   mean frac_add=%.4f"%(len(res),arr.mean(),arr.max(),fr.mean()))
# ==================================================================================================

refs = committed_references()
out = dict(config=dict(L=L,U=U,eta=eta,K=K,dt=dt,T=T_SHOTS,seeds=SEEDS,NL=NL,
                       window=[float(args.wmin),float(args.wmax)],npts=int(args.npts),
                       n_channels=K_CHANNELS,shots_all_channels=K_CHANNELS*T_SHOTS,
                       shots_per_slice=T_SHOTS//(K+1),rng="numpy default_rng(seed*1000+n), n = momentum index",
                       selection="union of determinants observed in multinomial draws of T//(K+1) shots per Born slice",
                       reconstruction="Haydock on the restricted operator, depth min(NL,|S|)",
                       reference="dense eigh of each sector + Lorentzian Lehmann sum on the grid"),
           mean_relL1=float(arr.mean()),max_relL1=float(arr.max()),rows=res)
out['summary'] = postprocess(res, refs, [args.wmin, args.wmax], args.npts, T_SHOTS)
out['references'] = refs
out['definitions'] = dict(
    ratio_mean_to_panel_b="summary.mean_relL1 / references.panel_b_mean_relL1 (Fig. 3(b) mean, same window); "
        "null if this run is not scored on the window and grid of Fig. 3(b)",
    understatement_of_transferred_penalty="1 - references.transferred_penalty_L8_S2180 / summary.ratio_mean_to_panel_b; "
        "written only at T=2.6e6 (null otherwise), whose observed subspace covers 0.839, just under the 0.850 of Fig. 3(a,b)")
out['provenance'] = dict(
    script="src/fig3_finite_shot.py", generator_sha256=sha256_lf(os.path.abspath(__file__)),
    generator_sha256_convention="sha256 of the file with line endings normalised to LF",
    recovered_from="2026-09-19 session transcript, subagent agent-a07c2e0927cf28823 (workflow wf_bb98610c-331), "
                   "lines 240/241 (fig5_literal.py), 625/626 (sed to T=5.2e6), 959/960 (compare step)",
    recovered_file_sha256_lf="e2b3793a2267953f4fab6adfb293e891033bc492c918069c2d27117a237873b5",
    generated_utc=datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    wall_seconds=round(time.time()-T0,1),
    argv=[os.path.basename(a) if os.path.isabs(a) else a for a in sys.argv[1:]], E0=E0,
    python=platform.python_version(), numpy=np.__version__, scipy=__import__('scipy').__version__,
    platform=platform.platform())
if args.repro:
    out['reproduction_check'] = reproduction_check(args.repro, T_SHOTS)
    log("reproduction check on [-9t,9t]:", "EXACT" if out['reproduction_check']['exact'] else "FAILED")
s = out['summary']
log("mean %.4e +- %.2e (ddof=1)  max %.4e | frac_add %.4f frac_rem %.4f" % (
    s['mean_relL1'], s['sd_relL1_ddof1'], s['max_relL1'], s['mean_frac_add'], s['mean_frac_rem']))
if s['ratio_mean_to_panel_b'] is None:
    log("ratio to Fig. 3(b): not computed (" + s['ratio_note'] + ")")
else:
    log("ratio to Fig. 3(b) mean %.2f max %.2f" % (s['ratio_mean_to_panel_b'], s['ratio_max_to_panel_b']))
if s['understatement_of_transferred_penalty'] is not None:
    log("understatement 1-%.4f/%.4f = %.4f" % (refs['transferred_penalty_L8_S2180'], s['ratio_mean_to_panel_b'],
                                               s['understatement_of_transferred_penalty']))
p = os.path.abspath(args.out or default_out(T_SHOTS))
os.makedirs(os.path.dirname(p), exist_ok=True)
json.dump(out, open(p, "w"), indent=1)
log("wrote", p)
if args.repro and not out['reproduction_check']['exact']:
    sys.exit(1)
log("done")
