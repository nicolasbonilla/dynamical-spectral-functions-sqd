# -*- coding: utf-8 -*-
"""Caption repairs for fig:thm1iii-violation.

(a) The header claimed "every number below is printed by build_n3.py".  Two numbers in the
    caption, 156 and 4000, are NOT emitted by that script (the provenance reviewer checked
    by grep and was right).  They are traceable, and this script traces them itself from the
    JSONs rather than leaving the claim standing.
(b) The Fig.-5 overlay is now a min/median/max whisker, so the caption must say so and must
    quote the medians it now draws.
"""
import io
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SCR = os.path.dirname(HERE)
CAP = os.path.join(HERE, "caption_thm1iii_violation.tex")
SUMF = os.path.join(HERE, "n3_fig5_summary.json")
OUT = os.path.join(SCR, "p2_calc", "out")

# ---- trace 156 and 4000 from the deposited JSONs, do not take them on trust -------------
teq = json.load(io.open(os.path.join(OUT, "cert_teqsci.json"), encoding="utf-8"))
shots = int(teq["shots"])
rows_teq = [r for r in teq["rows"] if r.get("relL1", 1.0) > 1e-10]

thm3 = None
for base in (SCR, os.path.join(SCR, "p2_calc"), OUT, HERE):
    p = os.path.join(base, "thm3_results.json")
    if os.path.exists(p):
        thm3 = json.load(io.open(p, encoding="utf-8"))
        break
assert thm3 is not None, "thm3_results.json not found"
rows3 = thm3["rows"] if isinstance(thm3, dict) and "rows" in thm3 else thm3
n_thm3_teq = sum(1 for r in rows3
                 if r.get("kind") == "teqsci" and r.get("relL1", 1.0) > 1e-10)
assert len(rows_teq) + n_thm3_teq == 156, (len(rows_teq), n_thm3_teq)
assert shots == 4000, shots

sm = json.load(io.open(SUMF, encoding="utf-8"))
med = [g["median"] for g in sm["frac_groups"]]
fr = [g["x"] for g in sm["frac_groups"]]

s = io.open(CAP, encoding="utf-8").read()

old_hdr = ("% Every number below is printed by build_n3.py from the JSON sources; "
           "none is typed from memory.")
L = []
L.append("% Every number below is traceable to the JSON sources.  The counts, factors and ratios")
L.append("% are printed by build_n3.py; the two exceptions are 156 and 4000, which that script does")
L.append("% NOT emit and which fix_n3_caption.py re-derives here from the deposited files:")
L.append("%%   156  = %d live rows of p2_calc/out/cert_teqsci.json  +  %d rows of thm3_results.json"
         % (len(rows_teq), n_thm3_teq))
L.append("%          carrying kind = teqsci")
L.append("%%   4000 = cert_teqsci.json /shots (= %d), the same value as thm3_run.py:28" % shots)
L.append("% The min/median/max whiskers of the Fig.-5 series come from n3_fig5_summary.json,")
L.append("% written by fix_n3_fig5_series.py from n3_fig5_frac.dat / n3_fig5_eta.dat.")
new_hdr = "\n".join(L)
assert old_hdr in s
s = s.replace(old_hdr, new_hdr)

old = "at the published fraction $0.85$."
assert old in s
new = ("at the published fraction $0.85$. Those $64$ rows sit at only four abscissae, so they are "
       "drawn as a min--median--max whisker at each one rather than as $64$ overlapping markers "
       "(medians $%s$ at fractions $%s$); each of them is also plotted individually inside the "
       "filled orange series, which is where their scatter can be read."
       % (", ".join("%.1f" % m for m in med), ", ".join("%.2f" % x for x in fr)))
s = s.replace(old, new)
io.open(CAP, "w", encoding="utf-8").write(s)
print("156 traced as %d + %d ; shots = %d" % (len(rows_teq), n_thm3_teq, shots))
print("medians", ["%.1f" % m for m in med])
print("patched", CAP)
