# -*- coding: utf-8 -*-
"""Replace the illegible 64-diamond overlay of fig_thm1iii_violation_native.tex by a
min / median / max summary at each distinct abscissa.

Reviewer defect D3: the 64 rows of the paper's Fig.-5 configuration sit at only four
subspace fractions (and at a single eta), so 2.2pt diamonds fuse into four solid black
blobs that (i) resolve no individual marker and (ii) occlude the orange points beneath.

No data is hidden by this change: all 64 rows remain plotted as part of the orange
w_S >= 0.99 series (n3_pub_frac_hi.dat / n3_pub_eta_hi.dat).  The black overlay is a
HIGHLIGHT of a subset, and it is now drawn as what the caption actually quotes: the
range of the violation factor at each fraction.

Every number written here is read from the .dat files that build_n3.py emitted from the
JSONs; nothing is typed.  The script asserts the row count and the caption's quoted
range before writing.
"""
import io
import os
import statistics

HERE = os.path.dirname(os.path.abspath(__file__))
TEX = os.path.join(HERE, "fig_thm1iii_violation_native.tex")
SUM = os.path.join(HERE, "n3_fig5_summary.json")


def read(name):
    rows = []
    with io.open(os.path.join(HERE, name), encoding="utf-8") as fh:
        for i, ln in enumerate(fh):
            ln = ln.strip()
            if not ln or i == 0:
                continue
            a, b = ln.split()
            rows.append((float(a), float(b)))
    return rows


def summarise(rows):
    groups = {}
    for x, y in rows:
        groups.setdefault(x, []).append(y)
    out = []
    for x in sorted(groups):
        v = sorted(groups[x])
        out.append({"x": x, "n": len(v), "min": v[0], "max": v[-1],
                    "median": statistics.median(v)})
    return out


frac = read("n3_fig5_frac.dat")
eta = read("n3_fig5_eta.dat")
assert len(frac) == 64 and len(eta) == 64, (len(frac), len(eta))
sf, se = summarise(frac), summarise(eta)
assert [g["x"] for g in sf] == [0.5, 0.7, 0.85, 0.95], [g["x"] for g in sf]
assert len(se) == 1 and se[0]["x"] == 0.18, se

allmin = min(y for _, y in frac)
allmax = max(y for _, y in frac)
at85 = [g for g in sf if g["x"] == 0.85][0]
# the caption quotes 15.1-57.3 overall and 15.1-33.7 at the published fraction 0.85
assert "%.1f" % allmin == "15.1" and "%.1f" % allmax == "57.3", (allmin, allmax)
assert "%.1f" % at85["min"] == "15.1" and "%.1f" % at85["max"] == "33.7", at85


def whisker(groups, cap):
    """Vertical min-max whisker with horizontal caps, plus an open median marker."""
    L = []
    for g in groups:
        x, lo, hi = g["x"], g["min"], g["max"]
        L.append(r"\draw[vlInk,line width=0.5pt] (axis cs:%.6g,%.6g) -- (axis cs:%.6g,%.6g);"
                 % (x, lo, x, hi))
        for yy in (lo, hi):
            L.append(r"\draw[vlInk,line width=0.5pt] (axis cs:%.6g,%.6g) -- (axis cs:%.6g,%.6g);"
                     % (x * (1 - cap), yy, x * (1 + cap), yy))
    pts = " ".join("(%.6g,%.6g)" % (g["x"], g["median"]) for g in groups)
    L.append(r"\addplot[forget plot,only marks,mark=diamond,mark size=3.0pt,draw=vlInk,"
             r"line width=0.6pt,mark options={fill=white,fill opacity=1}] coordinates {%s};" % pts)
    return L


s = io.open(TEX, encoding="utf-8").read()

old_a = (r"\addplot[forget plot,only marks,mark=diamond,mark size=2.2pt,draw=vlInk,line width=0.45pt]"
         "\n"
         r"  table {n3_fig5_frac.dat};")
old_b = (r"\addplot[only marks,mark=diamond,mark size=2.2pt,draw=vlInk,line width=0.45pt]"
         "\n"
         r"  table {n3_fig5_eta.dat};")
assert old_a in s and old_b in s

head = (r"%% --- paper's Fig.5 configuration (64 rows): drawn as min / median / max at each" "\n"
        r"%% abscissa.  The 64 individual rows are NOT removed from the figure -- they are part of" "\n"
        r"%% the orange w_S>=0.99 series.  As 64 separate diamonds they fused into solid blobs that" "\n"
        r"%% resolved no marker and hid the orange points underneath (reviewer defect D3)." "\n"
        r"%% Caps are drawn at +-%d%% of the abscissa; that is layout, the y values are the data." "\n")

new_a = (head % 3) + "\n".join(whisker(sf, 0.03))
new_b = (head % 6) + "\n".join(whisker(se, 0.06))
s = s.replace(old_a, new_a).replace(old_b, new_b)

# legend proxy must describe what is now drawn
old_leg = (r"\addlegendimage{only marks,mark=diamond,mark size=3.0pt,draw=vlInk,line width=0.7pt}"
           "\n"
           r"\addlegendentry{paper's Fig.~5 setting (64)}")
assert old_leg in s
new_leg = (r"\addlegendimage{only marks,mark=diamond,mark size=3.0pt,draw=vlInk,line width=0.7pt}"
           "\n"
           r"\addlegendentry{paper's Fig.~5 setting (64): min--median--max}")
s = s.replace(old_leg, new_leg)

io.open(TEX, "w", encoding="utf-8").write(s)

import json
json.dump({"frac_groups": sf, "eta_groups": se, "n_rows": len(frac),
           "overall_min": allmin, "overall_max": allmax,
           "at_published_fraction_0.85": at85,
           "source_dat": ["n3_fig5_frac.dat", "n3_fig5_eta.dat"]},
          io.open(SUM, "w", encoding="utf-8"), indent=1)

for g in sf:
    print("frac %.2f  n=%2d  min=%8.3f  med=%8.3f  max=%8.3f" % (g["x"], g["n"], g["min"], g["median"], g["max"]))
print("eta  %.2f  n=%2d  min=%8.3f  med=%8.3f  max=%8.3f"
      % (se[0]["x"], se[0]["n"], se[0]["min"], se[0]["median"], se[0]["max"]))
print("patched", TEX)
