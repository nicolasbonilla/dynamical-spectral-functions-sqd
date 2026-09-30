# -*- coding: utf-8 -*-
r"""Deviation of the plotted spin-peak dispersion of Fig. 8 (fig:spin) from the des Cloizeaux-Pearson
lower boundary (pi J / 2)|sin q|.  Added 2026-09-28: the caption's percentages (about 3% for
q <= pi/2, 10% at q = 2pi/3, 29% at q = 5pi/6) had no deposited script that printed them.

The peak is the one the cyan markers of Fig. 8 draw: the argmax of the broadened S^zz(q, omega) of
data/spinqw_L12.json on its own frequency grid (src/spin_lanczos.py, column `peak` of
paper/figs/spinqw_edges.dat), on the 0.006 t grid of that file; the caption quotes these
plotted (grid-argmax) values, not the positions of the unbroadened poles.

Run:  python src/spin_peak_deviation.py      -> data/spin_peak_deviation.json  (seconds)
"""
import os, sys, json, hashlib, datetime
import numpy as np

sys.dont_write_bytecode = True
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    d = json.load(open(os.path.join(_REPO, "data", "spinqw_L12.json")))
    wg = np.array(d["wg"])
    J = float(d["J"])
    rows = []
    for key in sorted(d["S"], key=float):
        qpi = float(key)
        A = np.array(d["S"][key])
        pk = float(wg[int(np.argmax(A))])
        dcp = 0.5 * np.pi * J * abs(np.sin(np.pi * qpi))
        rows.append(dict(q_over_pi=qpi, peak=pk, dCP_lower=float(dcp),
                         rel_dev=(pk - dcp) / dcp if dcp > 1e-12 else None))
    le = [abs(r["rel_dev"]) for r in rows
          if r["rel_dev"] is not None and (r["q_over_pi"] <= 0.5 + 1e-9 or r["q_over_pi"] >= 1.5 - 1e-9)]
    at = {("%.4f" % r["q_over_pi"]): r["rel_dev"] for r in rows}
    summ = dict(max_abs_rel_dev_q_le_half_pi=max(le),
                rel_dev_2pi_3=at["0.6667"], rel_dev_5pi_6=at["0.8333"],
                grid_spacing=float(wg[1] - wg[0]), J=J)
    for r in rows:
        print("q/pi=%.4f  peak %.4f  dCP %.4f  dev %s%s" % (
            r["q_over_pi"], r["peak"], r["dCP_lower"],
            "%+.2f%%" % (100 * r["rel_dev"]) if r["rel_dev"] is not None else "   n/a", ""))
    print("max |dev| for q <= pi/2 (and q >= 3pi/2): %.2f%%   2pi/3: %+.1f%%   5pi/6: %+.1f%%   grid %.4f t" % (
        100 * summ["max_abs_rel_dev_q_le_half_pi"], 100 * summ["rel_dev_2pi_3"], 100 * summ["rel_dev_5pi_6"],
        summ["grid_spacing"]))
    with open(os.path.abspath(__file__), "rb") as f:
        h = hashlib.sha256(f.read().replace(b"\r\n", b"\n")).hexdigest()
    out = dict(what="plotted spin-peak dispersion of Fig. 8 vs the des Cloizeaux-Pearson lower boundary",
               source="data/spinqw_L12.json (argmax on its grid, as the markers of paper/figs/spinqw_edges.dat)",
               rows=rows, summary=summ,
               provenance=dict(script="src/spin_peak_deviation.py", generator_sha256_lf=h,
                               generated_utc=datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")))
    json.dump(out, open(os.path.join(_REPO, "data", "spin_peak_deviation.json"), "w"), indent=1)
    print("wrote data/spin_peak_deviation.json")


if __name__ == "__main__":
    main()
