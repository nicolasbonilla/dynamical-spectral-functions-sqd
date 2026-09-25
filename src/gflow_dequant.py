# -*- coding: utf-8 -*-
r"""gflow_dequant.py -- the dequantization control of Sec. IX / App. dequant, per seed.

WHAT IT COMPUTES
----------------
N2, cc-pVDZ, R = 2.5 A, CAS(10e,12o): 792 alpha (and 792 beta) strings.  A device head is
emulated by drawing `shots` configurations from the exact FCI string marginal and corrupting
them with a Heron-class noise model (5% depolarizing replacement plus asymmetric readout,
1->0 and 0->1 rates read from qiskit's FakeTorino calibration when available).  From the
SAME corrupted histogram five subspaces of D = 120 strings per spin are built and scored on
the exact FCI energy (selected-CI diagonalization in the subspace):

    raw          valid-particle-number shots only, ranked by count
    ibm          IBM-style occupancy recovery of every shot, ranked by count
    ibm+cheap    ibm, then filled to D with the top Epstein--Nesbet PT1 strings from HF
                 (one application of H, no FCI, no further shots)       -- classical control
    gfn-nocheap  a GFlowNet (trajectory balance) trained on the recovered histogram only
    gfn-fused    a GFlowNet trained on the histogram fused with the EN-PT1 prior

The shot model, the recovery, the classical fill and the fused GFlowNet are exactly those of the
companion review's published noise ladder (calculations/n2_ladder_crossover.py: symmetry-pinned
gauge, 800 training iterations, 15000 draws, seeds 0-4 and 5-9), so at R = 2.5 A and 1x noise the
ibm+cheap minus gfn-fused gap reproduces that paper's two runs seed for seed; the raw, ibm and
gfn-nocheap arms are those of validate_realnoise.py of the same review (Bonilla Vargas,
arXiv:2608.05314; repository github.com/nicolasbonilla/ml-for-sqd-review).  Every seed is
WRITTEN, so the paired comparisons the paper makes are computed, not bounded.  Seeds 0-4
reproduce the companion's first ladder run to four decimals (+5.2314 +- 3.1891 mHa); seeds 5-9
give -1.62 against its -1.83.  (v2 of this paper printed five means from an unpinned-gauge run
of validate_realnoise.py whose per-seed values were never stored; superseded 2026-09-25.)

DETERMINISM
-----------
Run single-threaded (OMP/OPENBLAS/MKL threads = 1, torch.set_num_threads(1)).  With default
threading, PySCF's RHF picks a thread-dependent orientation inside the degenerate pi shells
of N2 (PySCF issue #1717), which changes the string basis and every number downstream.  The
per-seed values of the original run were not stored and are not reproduced bit for bit; the
numbers the paper prints are the ones written by THIS script.

USAGE (pyscf and torch are needed; the repository's Docker image sqd-nb has both)
    python src/gflow_dequant.py --shots 1000 --seed 0      # one seed -> data/gflow_runs/
    python src/gflow_dequant.py --merge                    # -> data/gflow.json
About 15 minutes per run single-threaded; do not launch all twenty at once (memory).  No QPU.
"""
import argparse
import datetime
import json
import os
import platform
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
RUNS = os.path.join(REPO, "data", "gflow_runs")
OUT = os.path.join(REPO, "data", "gflow.json")
METHODS = ["raw", "ibm", "ibm+cheap", "gfn-nocheap", "gfn-fused"]
SEEDS = list(range(10))      # 0-4 and 5-9: the two runs of the companion's ladder
D, BETA_T, R_BOND = 120, 0.5, 2.5
ITERS, NDRAW = 800, 15000    # the companion's published ladder (n2_ladder_crossover.py)
t0 = time.time()


def log(*a):
    print("[%7.1fs]" % (time.time() - t0), *a, flush=True)


def run(shots, seed):
    import torch
    import torch.nn as nn
    import pyscf
    from pyscf import gto, scf, mcscf, ao2mo
    from pyscf.fci import cistring, selected_ci, direct_spin1
    torch.set_num_threads(1)

    # noise rates from FakeTorino (fallback Heron-typical), exactly as validate_realnoise.py
    LAMBDA = 0.05
    P10, P01 = 0.0229, 0.0200
    noise_source = "fallback constants"
    try:
        from qiskit_ibm_runtime.fake_provider import FakeTorino
        props = FakeTorino().properties()
        nq = FakeTorino().num_qubits
        P10 = float(np.median([props.qubit_property(q, "prob_meas0_prep1")[0] for q in range(nq)]))
        P01 = float(np.median([props.qubit_property(q, "prob_meas1_prep0")[0] for q in range(nq)]))
        noise_source = "FakeTorino median readout"
    except Exception:
        pass
    log("noise: 1->0=%.4f 0->1=%.4f depol=%.2f (%s)" % (P10, P01, LAMBDA, noise_source))

    NCAS, NELECAS = 12, (5, 5)
    na = NELECAS[0]
    # symmetry=True pins the orbital gauge inside the degenerate pi shells (D-infinity-h), as in
    # the companion's ladder; without it the string basis is not defined run to run
    mol = gto.M(atom="N 0 0 0; N 0 0 %s" % R_BOND, basis="cc-pvdz", symmetry=True, verbose=0)
    mf = scf.RHF(mol).run()
    cas = mcscf.CASCI(mf, NCAS, NELECAS)
    h1, ecore = cas.get_h1cas()
    h2 = ao2mo.restore(1, cas.get_h2cas(), NCAS)
    strs_a = cistring.make_strings(range(NCAS), na)
    dim_a = len(strs_a)
    str_to_idx = {int(s): i for i, s in enumerate(strs_a)}
    set_to_idx = {frozenset(p for p in range(NCAS) if (int(s) >> p) & 1): i
                  for i, s in enumerate(strs_a)}
    HF = (1 << na) - 1
    hf_idx = str_to_idx[HF]
    civ_hf = np.zeros((dim_a, dim_a))
    civ_hf[hf_idx, hf_idx] = 1.0
    h2e = direct_spin1.absorb_h1e(h1, h2, NCAS, NELECAS, 0.5)
    Hc = direct_spin1.contract_2e(h2e, civ_hf, NCAS, NELECAS).reshape(dim_a, dim_a)
    hdiag = direct_spin1.make_hdiag(h1, h2, NCAS, NELECAS).reshape(dim_a, dim_a)
    den = Hc[hf_idx, hf_idx] - hdiag
    den[hf_idx, hf_idx] = 1.0
    c1 = Hc / den
    c1[hf_idx, hf_idx] = 1.0
    w_cheap = (c1 ** 2).sum(1)
    w_cheap = w_cheap / w_cheap.sum()
    cheap_order = np.argsort(w_cheap)[::-1]
    e_fci, civec = pyscf.fci.direct_spin1.FCI().kernel(h1, h2, NCAS, NELECAS, ecore=ecore)
    civec = civec.reshape(dim_a, dim_a)
    w_true = (civec ** 2).sum(1)
    w_true /= w_true.sum()
    log("E_FCI = %.8f Ha  (%d strings per spin)" % (e_fci, dim_a))
    _sci = selected_ci.SelectedCI()

    def E(idxs):
        s = np.asarray(sorted(set(int(strs_a[i]) for i in idxs)), dtype=np.int64)
        out = selected_ci.kernel_fixed_space(_sci, h1, h2, NCAS, NELECAS, (s, s), ecore=ecore)
        return (float(out[0] if isinstance(out, (tuple, list)) else out) - e_fci) * 1000

    class Policy(nn.Module):
        def __init__(s, n):
            super().__init__()
            s.net = nn.Sequential(nn.Linear(n, 256), nn.ReLU(), nn.Linear(256, 256), nn.ReLU(),
                                  nn.Linear(256, n))
            s.logZ = nn.Parameter(torch.zeros(1))

        def forward(s, x):
            return s.net(x)

    def sample_batch(net, B):
        state = torch.zeros(B, NCAS)
        logPF = torch.zeros(B)
        for _ in range(na):
            lp = torch.log_softmax(net(state).masked_fill(state.bool(), -1e9), 1)
            a = torch.distributions.Categorical(logits=lp).sample()
            logPF += lp.gather(1, a[:, None]).squeeze(1)
            state = state.scatter(1, a[:, None], 1.0)
        return np.array([set_to_idx[frozenset(np.where(state[b].numpy() > 0)[0].tolist())]
                         for b in range(B)]), logPF

    def train_and_draw(reward, ndraw=NDRAW, iters=ITERS):
        net = Policy(NCAS)
        opt = torch.optim.Adam([{"params": net.net.parameters(), "lr": 1e-3},
                                {"params": [net.logZ], "lr": 1e-1}])
        Rt = torch.tensor(reward / reward.sum(), dtype=torch.float32)
        for it in range(iters):
            idxs, logPF = sample_batch(net, 256)
            loss = ((net.logZ + logPF - torch.log(Rt[idxs] + 1e-12)) ** 2).mean()
            opt.zero_grad()
            loss.backward()
            opt.step()
        o = []
        while len(o) < ndraw:
            idxs, _ = sample_batch(net, 512)
            o += idxs.tolist()
        c = {}
        for i in o:
            c[i] = c.get(i, 0) + 1
        return c

    def b2s(b):
        return int(sum(int(v) << p for p, v in enumerate(b)))

    def fill_to_D(ranked):
        sel = list(dict.fromkeys(ranked))[:D]
        if len(sel) < D:
            for i in cheap_order:
                if i not in sel:
                    sel.append(int(i))
                if len(sel) >= D:
                    break
        return sel[:D]

    rng = np.random.default_rng(seed)
    torch.manual_seed(seed)
    ideal = rng.choice(dim_a, size=shots, p=w_true)
    bits = np.array([[(int(strs_a[i]) >> p) & 1 for p in range(NCAS)] for i in ideal], dtype=np.int8)
    dep = rng.random(shots) < LAMBDA
    bits[dep] = (rng.random((dep.sum(), NCAS)) < 0.5).astype(np.int8)
    o, z = bits == 1, bits == 0
    bits[o & (rng.random(bits.shape) < P10)] = 0
    bits[z & (rng.random(bits.shape) < P01)] = 1
    good = np.array([int(b.sum()) == na for b in bits])
    craw = {}
    for b in bits[good]:
        i = str_to_idx[b2s(b)]
        craw[i] = craw.get(i, 0) + 1
    occ = bits[good].mean(0) if good.any() else np.full(NCAS, na / NCAS)
    rec = []
    for row in bits:
        s = row.copy()
        m = int(s.sum())
        if m == na:
            rec.append(b2s(s))
            continue
        if m > na:
            od = np.where(s == 1)[0]
            s[od[np.argsort(occ[od])[:m - na]]] = 0
        else:
            em = np.where(s == 0)[0]
            s[em[np.argsort(occ[em])[::-1][:na - m]]] = 1
        rec.append(b2s(s))
    cibm = {}
    for s in rec:
        i = str_to_idx[s]
        cibm[i] = cibm.get(i, 0) + 1
    f_emp = np.zeros(dim_a)
    for i, cc in cibm.items():
        f_emp[i] = cc
    f_emp = f_emp / max(f_emp.sum(), 1)
    r_fused = np.maximum(f_emp + 0.1 * w_cheap / w_cheap.max() * max(f_emp.max(), 1e-9),
                         1e-9) ** BETA_T
    r_noch = np.maximum(f_emp, 1e-9) ** BETA_T
    log("seed %d: %d of %d shots have the right particle number; training two GFlowNets"
        % (seed, int(good.sum()), shots))
    gfn_f = train_and_draw(r_fused)
    gfn_n = train_and_draw(r_noch)

    def rank(c):
        return [i for i, _ in sorted(c.items(), key=lambda kv: -kv[1])]

    res = {
        "raw": E(rank(craw)[:D]) if craw else float("nan"),
        "ibm": E(rank(cibm)[:D]),
        "ibm+cheap": E(fill_to_D(rank(cibm))),
        "gfn-nocheap": E(rank(gfn_n)[:D]),
        "gfn-fused": E(rank(gfn_f)[:D]),
    }
    log("seed %d: " % seed + "  ".join("%s %.2f" % (k, v) for k, v in res.items()))
    return dict(shots=shots, seed=seed, error_mHa=res, E_FCI_Ha=float(e_fci),
                valid_shots=int(good.sum()), noise=dict(p10=P10, p01=P01, depol=LAMBDA,
                                                        source=noise_source),
                versions=dict(python=platform.python_version(), numpy=np.__version__,
                              pyscf=pyscf.__version__, torch=torch.__version__),
                threads=1, seconds=round(time.time() - t0, 1))


def _t_paired(a, b):
    d = np.asarray(a, float) - np.asarray(b, float)
    n = d.size
    sd = d.std(ddof=1)
    t = float(d.mean() / (sd / np.sqrt(n))) if sd > 0 else float("inf")
    from math import comb
    k = int((d > 0).sum())
    tail = sum(comb(n, j) for j in range(0, min(k, n - k) + 1)) / 2 ** n
    return dict(mean_difference=float(d.mean()), sd_difference=float(sd), t=t, df=n - 1,
                n_positive=k, sign_test_two_sided_p=float(min(1.0, 2 * tail)))


def merge():
    blocks = {}
    for f in sorted(os.listdir(RUNS)):
        if f.endswith(".json"):
            r = json.load(open(os.path.join(RUNS, f), encoding="utf-8"))
            blocks.setdefault(r["shots"], {})[r["seed"]] = r
    out = dict(generated_by="src/gflow_dequant.py --merge",
               generated_utc=datetime.datetime.now(datetime.timezone.utc)
               .strftime("%Y-%m-%dT%H:%M:%SZ"),
               system="N2 cc-pVDZ R=2.5 A, CAS(10e,12o), 792 strings per spin; subspace D=120 "
                      "strings per spin (14400 of 627264 determinants)",
               error_definition="E_subspace - E_FCI in mHa (lower is better)",
               source_of_algorithm="validate_realnoise.py of the companion review, unchanged",
               runs={})
    for shots, seeds in sorted(blocks.items()):
        got = sorted(seeds)
        per = {m: [seeds[s]["error_mHa"][m] for s in got] for m in METHODS}
        summ = {m: dict(mean=float(np.mean(v)), std_pop=float(np.std(v)),
                        std_sample=float(np.std(v, ddof=1)), per_seed=v) for m, v in per.items()}
        out["runs"][str(shots)] = dict(
            seeds=got, methods=summ,
            paired=dict(
                ibm_minus_ibmcheap=_t_paired(per["ibm"], per["ibm+cheap"]),
                ibmcheap_minus_gfnfused=_t_paired(per["ibm+cheap"], per["gfn-fused"]),
                ibmcheap_minus_gfnnocheap=_t_paired(per["ibm+cheap"], per["gfn-nocheap"])),
            noise=seeds[got[0]]["noise"], versions=seeds[got[0]]["versions"],
            E_FCI_Ha=seeds[got[0]]["E_FCI_Ha"],
            valid_shots=[seeds[s]["valid_shots"] for s in got])
        log("shots=%d  seeds %s" % (shots, got))
        for m in METHODS:
            s = summ[m]
            log("  %-12s %7.2f +- %5.2f (pop)  per seed %s"
                % (m, s["mean"], s["std_pop"], " ".join("%.2f" % x for x in s["per_seed"])))
        for k, p in out["runs"][str(shots)]["paired"].items():
            log("  paired %-26s d=%.2f sd=%.2f t4=%.2f sign-test p=%.4f"
                % (k, p["mean_difference"], p["sd_difference"], p["t"], p["sign_test_two_sided_p"]))
    json.dump(out, open(OUT, "w", encoding="utf-8"), indent=1)
    log("WROTE", os.path.relpath(OUT, REPO))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--shots", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--merge", action="store_true")
    a = ap.parse_args()
    if a.merge:
        return merge()
    os.makedirs(RUNS, exist_ok=True)
    for sd in ([a.seed] if a.seed is not None else SEEDS):
        r = run(a.shots, sd)
        json.dump(r, open(os.path.join(RUNS, "shots%d_seed%d.json" % (a.shots, sd)), "w",
                          encoding="utf-8"), indent=1)


if __name__ == "__main__":
    main()
