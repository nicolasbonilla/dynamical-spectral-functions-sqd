# -*- coding: utf-8 -*-
"""Min / median / max whiskers of the Fig. 3-setting rows in Fig. S5 (fig:thm1iii-violation).

Added 2026-09-28, when Fig. S5 was rebuilt from the deposited certificate scans (build_n3.py).
It does what src/recovered/fix_n3_fig5_series.py did once, on 2026-09-18, as a one-off patch of
the fragment: it reads the two tables build_n3.py writes for the 64 rows of the Fig. 3 setting
(paper/figs/n3_fig5_frac.dat, n3_fig5_eta.dat), writes their min / median / max at each abscissa
to data/thm1iii_violation/n3_fig5_summary.json, and rewrites the two whisker blocks of
paper/figs/fig_thm1iii_violation_native.tex from them (same drawing code, caps at +-3% and +-6%
of the abscissa).  The medians it prints are the ones the Fig. S5 caption quotes.

Run after build_n3.py:   python src/recovered/n3_fig5_whiskers.py
"""
import io, os, json, statistics

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
FIGS = os.path.join(REPO, "paper", "figs")
TEX = os.path.join(FIGS, "fig_thm1iii_violation_native.tex")
SUM = os.path.join(REPO, "data", "thm1iii_violation", "n3_fig5_summary.json")


def read(name):
    rows = []
    with io.open(os.path.join(FIGS, name), encoding="utf-8") as fh:
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
    return [{"x": x, "n": len(v), "min": min(v), "max": max(v), "median": statistics.median(v)}
            for x, v in sorted(groups.items())]


def whisker(groups, cap):
    L = []
    for g in groups:
        x, lo, hi = g["x"], g["min"], g["max"]
        L.append(r"\draw[vlInk,line width=0.5pt] (axis cs:%.6g,%.6g) -- (axis cs:%.6g,%.6g);" % (x, lo, x, hi))
        for yy in (lo, hi):
            L.append(r"\draw[vlInk,line width=0.5pt] (axis cs:%.6g,%.6g) -- (axis cs:%.6g,%.6g);"
                     % (x * (1 - cap), yy, x * (1 + cap), yy))
    pts = " ".join("(%.6g,%.6g)" % (g["x"], g["median"]) for g in groups)
    L.append(r"\addplot[forget plot,only marks,mark=diamond,mark size=3.0pt,draw=vlInk,"
             r"line width=0.6pt,mark options={fill=white,fill opacity=1}] coordinates {%s};" % pts)
    return L


frac, eta = read("n3_fig5_frac.dat"), read("n3_fig5_eta.dat")
assert len(frac) == 64 and len(eta) == 64, (len(frac), len(eta))
sf, se = summarise(frac), summarise(eta)
assert [g["x"] for g in sf] == [0.5, 0.7, 0.85, 0.95] and [g["x"] for g in se] == [0.18]

# rewrite the two whisker blocks: maximal runs of whisker \draw lines closed by the median addplot
lines = io.open(TEX, encoding="utf-8").read().split("\n")
IS_W = lambda s: s.startswith(r"\draw[vlInk,line width=0.5pt] (axis cs:")
IS_D = lambda s: s.startswith(r"\addplot[forget plot,only marks,mark=diamond,mark size=3.0pt,draw=vlInk,")
blocks, i = [], 0
while i < len(lines):
    if IS_W(lines[i]):
        j = i
        while j < len(lines) and IS_W(lines[j]):
            j += 1
        assert j < len(lines) and IS_D(lines[j]), "whisker block at line %d not closed by the median plot" % (i + 1)
        blocks.append((i, j + 1)); i = j + 1
    else:
        i += 1
assert len(blocks) == 2, blocks
new = lines[:blocks[0][0]] + whisker(sf, 0.03) + lines[blocks[0][1]:blocks[1][0]] + whisker(se, 0.06) + lines[blocks[1][1]:]
io.open(TEX, "w", encoding="utf-8", newline="\n").write("\n".join(new))

at85 = [g for g in sf if g["x"] == 0.85][0]
out = {"frac_groups": sf, "eta_groups": se, "n_rows": len(frac),
       "overall_min": min(y for _, y in frac), "overall_max": max(y for _, y in frac),
       "at_published_fraction_0.85": at85,
       "source_dat": ["n3_fig5_frac.dat", "n3_fig5_eta.dat"],
       "written_by": "src/recovered/n3_fig5_whiskers.py (2026-09-28; formerly fix_n3_fig5_series.py)"}
json.dump(out, io.open(SUM, "w", encoding="utf-8"), indent=1)
print("medians at 0.50/0.70/0.85/0.95: %s" % ", ".join("%.1f" % g["median"] for g in sf))
print("overall %.1f-%.1f; at 0.85 %.1f-%.1f; eta group median %.1f" % (
    out["overall_min"], out["overall_max"], at85["min"], at85["max"], se[0]["median"]))
