# -*- coding: utf-8 -*-
r"""Derive the figure -> source -> generator -> data map from the MANUSCRIPT, and check
docs/FIGURE_PROVENANCE.md against it.

WHY THIS FILE EXISTS
--------------------
docs/FIGURE_PROVENANCE.md is hand-written, and on 2026-09-19 it was measured to document
FIVE figures the manuscript no longer typesets (fig:amort, fig:bench, fig:gallery,
fig:gflow, fig:molsuite) and to omit TWO it does (fig:gapscaling, fig:thm1iii-violation).
docs/FILE_INDEX.md described the v2 file set in its entirety.  Neither document was wrong
when it was written; both went stale silently, because nothing derived them from the
manuscript.  This file derives them.

WHAT IT DOES
------------
1. Expands `\input` from paper/main.tex INLINE and with LaTeX comment rules, producing the
   set of source files the document actually reaches.
2. Censuses every float in that stream: environment, starred or not, its \label, the
   files it is assembled from, its \includegraphics, and the pgfplots `table{...}` data
   tables its fragments read.
3. Checks the float census against the count printed in the main.tex preamble comment
   ("N of the F+T floats are starred"), which was false until this file was written.
4. Checks docs/FIGURE_PROVENANCE.md: the set of `fig:` labels it documents must be
   EXACTLY the set the manuscript typesets, and every file it names in the source,
   generator and data columns must exist unless it is marked `NOT DEPOSITED`.
5. Reports the orphan inventory in three tiers (orphan / dead-end / unsourced).

WHAT IT DOES NOT DO
-------------------
It does not check numbers.  `src/verify.py` does that.  It does not re-run generators.
`src/check_figures.py` does that -- see the note about its one blind spot in
docs/KNOWN_DISCREPANCIES.md.  This file checks only that the DOCUMENTATION still
describes the manuscript that is in the tree.

SELF-TEST FIRST
---------------
A PASS is worth nothing without a control that fires.  Nine synthetic controls run
before anything else and this file exits 2, reporting nothing, if any of them does not
behave as stated.  (The nine are a list in `selftest()`; the number printed is its
length, not a typed constant.)  Two of them are the exact mistakes made while writing it:
a `\begin{figure}` inside a `%` comment (paper/figs/fig_gapscaling_caption.tex really
has one, and counting it gives 18 figures instead of 17), and an `\input` sharing a line
with `\begin{figure*}` (paper/carried/*.tex are all written that way, and a line-level
expander drops the float entirely and finds 2 floats instead of 33).

    python src/check_provenance.py            # exit 0 iff the documentation still fits
    python src/check_provenance.py --markdown # emit the provenance table, ready to paste
    python src/check_provenance.py -v         # also print the orphan inventory in full
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.join(ROOT, "paper")

RE_INPUT = re.compile(r"\\(?:input|include)\s*\{([^}]*)\}")
RE_FLOAT_BEGIN = re.compile(r"\\begin\{(figure\*?|table\*?)\}")
RE_FLOAT_END = re.compile(r"\\end\{(figure\*?|table\*?)\}")
RE_LABEL = re.compile(r"\\label\{([^}]*)\}")
RE_GRAPHIC = re.compile(r"\\includegraphics(?:\[[^\]]*\])?\s*\{([^}]*)\}")
RE_PGFTABLE = re.compile(r"\btable\s*(?:\[[^\]]*\])?\s*\{([^}]*\.dat)\}")


# ---------------------------------------------------------------------------
# 1.  The expander.  Two rules that a naive version gets wrong, both of which
#     occur in this manuscript and both of which are controlled for below.
# ---------------------------------------------------------------------------
def strip_comment(line):
    r"""Everything from an unescaped % to end of line is not source. `\%` is a percent."""
    out = []
    i = 0
    while i < len(line):
        c = line[i]
        if c == "\\" and i + 1 < len(line):
            out.append(line[i:i + 2])
            i += 2
            continue
        if c == "%":
            break
        out.append(c)
        i += 1
    return "".join(out)


def expand(read_file, start, order=None, missing=None):
    r"""Return the document as [(file, line, text)] chunks, `\input` spliced INLINE.

    `read_file(rel)` returns the text of a path relative to paper/, or None.
    Splicing inline (not line by line) is what keeps `\begin{figure*}...\input{...}`
    on one line from losing its float.
    """
    order = [] if order is None else order
    missing = [] if missing is None else missing
    chunks = []

    def resolve(name):
        for cand in (name, name + ".tex"):
            if read_file(cand) is not None:
                return cand
        return None

    def walk(rel, stack):
        if rel in stack:                       # an \input cycle; do not recurse forever
            return
        if rel not in order:
            order.append(rel)
        raw = read_file(rel).replace("\r\n", "\n").replace("\r", "\n")
        for ln, line in enumerate(raw.split("\n"), 1):
            code = strip_comment(line)
            pos = 0
            while True:
                m = RE_INPUT.search(code, pos)
                if not m:
                    chunks.append((rel, ln, code[pos:]))
                    break
                chunks.append((rel, ln, code[pos:m.start()]))
                tgt = resolve(m.group(1))
                if tgt is None:
                    missing.append((rel, ln, m.group(1)))
                else:
                    walk(tgt, stack + [rel])
                pos = m.end()
            chunks.append((rel, ln, "\n"))

    walk(start, [])
    return chunks, order, missing


def census(chunks):
    """Every top-level float in the expanded stream, in order of appearance."""
    floats, depth, cur = [], 0, None
    for rel, ln, text in chunks:
        for tok in re.split(r"(?=\\begin\{|\\end\{)", text):
            mb = RE_FLOAT_BEGIN.match(tok)
            if mb:
                if depth == 0:
                    cur = dict(env=mb.group(1), file=rel, line=ln,
                               labels=[], graphics=[], tables=[], files=[])
                depth += 1
            if cur is not None and depth > 0:
                if rel not in cur["files"]:
                    cur["files"].append(rel)
                cur["labels"] += RE_LABEL.findall(tok)
                cur["graphics"] += RE_GRAPHIC.findall(tok)
                cur["tables"] += [t.split("/")[-1] for t in RE_PGFTABLE.findall(tok)]
            if RE_FLOAT_END.match(tok):
                depth -= 1
                if depth == 0 and cur is not None:
                    floats.append(cur)
                    cur = None
    return floats


# ---------------------------------------------------------------------------
# 2.  The self-test.  Runs before anything else; exit 2 if a control misbehaves.
# ---------------------------------------------------------------------------
SELFTEST = {
    # (b) an \input sharing a line with \begin{figure*}: the real shape of paper/carried/*.tex
    "main": (r"\begin{figure*}[tp]\centering\input{frag}" "\n"
             r"\caption{c}\label{fig:spliced}\end{figure*}" "\n"
             r"\input{commented}" "\n"
             r"\input{unstarred}" "\n"
             r"\input{nolabel}" "\n"
             r"100\% of the sector \begin{figure*}\caption{p}\label{fig:pct}\end{figure*}" "\n"),
    "frag": r"\addplot table[x=a,y=b]{some_table.dat};" "\n",
    # (a) a float declared inside a comment: the real shape of fig_gapscaling_caption.tex
    "commented": ("% Use as: " + r"\begin{figure}[p]\centering\input{x}\end{figure}" "\n"),
    # (c) an unstarred float must be seen AND reported as unstarred
    "unstarred": (r"\begin{figure}\caption{u}\label{fig:unstarred}\end{figure}" "\n"),
    # (d) a float with no \label at all
    "nolabel": (r"\begin{figure*}\caption{n}\end{figure*}" "\n"),
}


def selftest():
    """Return (problems, n_controls).  The count is derived, not typed."""
    def rd(rel):
        return SELFTEST.get(rel.replace(".tex", ""), SELFTEST.get(rel))

    chunks, order, missing = expand(rd, "main")
    fl = census(chunks)
    labels = [l for f in fl for l in f["labels"]]
    starred = {l: f["env"].endswith("*") for f in fl for l in f["labels"]}

    controls = [
        (len(fl) == 4,
         "the fixture's four floats were not all found (got %d: %s)" % (len(fl), labels)),
        ("fig:spliced" in labels,
         "(b) an \\input sharing a line with \\begin{figure*} lost its float; a line-level "
         "expander finds 2 floats in this manuscript instead of 33"),
        (not any(f["labels"] for f in fl if f["file"] == "commented"),
         "(a) a \\begin{figure} inside a %-comment was counted; counting the one in "
         "fig_gapscaling_caption.tex gives 18 figures instead of 17"),
        (starred.get("fig:unstarred") is False,
         "(c) an unstarred float was not reported as unstarred"),
        (starred.get("fig:spliced") is True,
         "(c) a starred float was not reported as starred"),
        (any(not f["labels"] for f in fl),
         "(d) a float with no \\label was not noticed"),
        ("fig:pct" in labels,
         "(e) an escaped \\%% was treated as the start of a comment"),
        (any("some_table.dat" in f["tables"] for f in fl),
         "(f) a pgfplots table{...} inside an \\input'ed fragment was not attributed to its float"),
        (not missing,
         "the self-test fixture is broken: unresolved inputs %s" % (missing,)),
    ]
    return [msg for ok, msg in controls if not ok], len(controls)


# ---------------------------------------------------------------------------
# 3.  The real run
# ---------------------------------------------------------------------------
def read_paper(rel):
    p = os.path.join(PAPER, rel.replace("/", os.sep))
    if not os.path.isfile(p):
        return None
    return open(p, "r", encoding="utf-8", errors="replace", newline="").read()


def rel_files(subdir, exts):
    d = os.path.join(ROOT, subdir.replace("/", os.sep))
    if not os.path.isdir(d):
        return []
    return sorted(n for n in os.listdir(d)
                  if os.path.isfile(os.path.join(d, n)) and os.path.splitext(n)[1] in exts)


# Files whose only job is to check the deposit.  A checker naming an artefact proves it
# is GUARDED, not that it is USED, so the tier annotation says so instead of letting
# them make every orphan look employed.
CHECKERS = {"src/verify.py", "src/check_figures.py", "src/check_provenance.py",
            "src/check_tarball.py", "src/build_arxiv_bundle.py",
            "src/build_repro_notebook.py"}

# Deliberately NOT a producer/consumer heuristic.  Writes in this repository happen
# through variables (`json.dump(d, open(_dest,'w'))`, `DEFAULT_OUT = os.path.join(...)`,
# `args.out`), so a line-local "is this a write?" test is wrong often enough to be
# worse than useless -- and a provenance tool that guesses is exactly how the documents
# this file exists to fix went stale.  The tiers below use only two facts that can be
# established without guessing: does any deposited FILE name the artefact, and does the
# MANUSCRIPT reach it.


def named_by(name, texts):
    return {src for src, txt in texts.items() if name in txt}


# A cell may name a file that is deliberately not here -- a standalone source kept
# outside the deposit, a generator that was never written down, a table that was
# deleted.  The cell must SAY so, and the exemption is granted per CELL, not per row:
# one honest "NOT DEPOSITED" must not license every other filename in the same row.
RE_EXEMPT = re.compile(r"NOT DEPOSITED|NOT IN THE DEPOSIT|DELETED|SUPERSEDED|RETIRED"
                       r"|NOT DEPOSITED HERE", re.I)
RE_FNAME = re.compile(r"`([A-Za-z0-9_./*-]+\.(?:tex|py|pdf|json|dat|png))`")


def parse_provenance_doc(path):
    """{label: [(filename, exempt?)]} over the cells of each table row."""
    rows = {}
    if not os.path.isfile(path):
        return None
    for line in open(path, encoding="utf-8", errors="replace"):
        if not line.lstrip().startswith("|"):
            continue
        m = re.search(r"\bfig:[A-Za-z0-9:_-]+", line)
        if not m:
            continue
        files = []
        for cell in line.split("|"):
            ex = bool(RE_EXEMPT.search(cell))
            files += [(fn, ex) for fn in RE_FNAME.findall(cell)]
        rows[m.group(0)] = dict(files=files, line=line.rstrip())
    return rows


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    verbose = "-v" in argv or "--verbose" in argv
    as_md = "--markdown" in argv

    problems, ncontrols = selftest()
    if problems:
        print("SELF-TEST FAILED (%d of %d controls) -- this file's own parser is wrong, so "
              "nothing it would report can be trusted:" % (len(problems), ncontrols))
        for p in problems:
            print("  * %s" % p)
        return 2
    if verbose:
        print("self-test: %d controls, all behaved as stated\n" % ncontrols)

    chunks, order, missing = expand(read_paper, "main.tex")
    floats = census(chunks)
    figs = [f for f in floats if f["env"].startswith("figure")]
    tabs = [f for f in floats if f["env"].startswith("table")]
    nstar = sum(1 for f in floats if f["env"].endswith("*"))

    reached_tex = {"paper/figs/" + o.split("/")[-1] for o in order if o.startswith("figs/")}
    reached_tex |= {"paper/" + o for o in order if not o.startswith("figs/")}
    reached_dat = {"paper/figs/" + t for f in floats for t in f["tables"]}
    reached_pdf = {"paper/" + g for f in floats for g in f["graphics"]}
    reached = reached_tex | reached_dat | reached_pdf

    fail = []
    if missing:
        fail.append("paper/main.tex reaches \\input files that do not exist: %s" % missing)
    nolabel = [f for f in floats if not f["labels"]]
    if nolabel:
        fail.append("floats with no \\label (nothing can \\ref them): %s"
                    % [(f["file"], f["line"]) for f in nolabel])

    print("MANUSCRIPT (derived from paper/main.tex by inline \\input expansion)")
    print("  source files reached      : %d" % len(order))
    print("  floats                    : %d  (%d figure, %d table)"
          % (len(floats), len(figs), len(tabs)))
    print("  double-column (starred)   : %d of %d" % (nstar, len(floats)))

    # -- 3a. the preamble comment must state the measured census -----------------
    mt = read_paper("main.tex") or ""
    m = re.search(r"(\d+)\s+of\s+the\s+(\d+)\+(\d+)\s+floats are starred", mt)
    if m is None:
        fail.append("the float-census sentence of the main.tex preamble could not be found; "
                    "it read '35 of the 19+16 floats are starred' and was false")
    else:
        got = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
        want = (nstar, len(figs), len(tabs))
        if got != want:
            fail.append("paper/main.tex says '%d of the %d+%d floats are starred'; the "
                        "manuscript has %d of %d+%d" % (got + want))
        else:
            print("  main.tex preamble census  : agrees (%d of %d+%d)" % want)

    # -- 3b. FIGURE_PROVENANCE.md must cover exactly the figures that exist -------
    derived = []
    for i, f in enumerate(figs, 1):
        derived.append(dict(n=i, label=f["labels"][0] if f["labels"] else "(none)",
                            decl="%s:%d" % (f["file"], f["line"]),
                            files=[x for x in f["files"] if x != f["file"]] or [],
                            wrapper=f["file"],
                            graphics=f["graphics"], tables=sorted(set(f["tables"]))))
    doc = parse_provenance_doc(os.path.join(ROOT, "docs", "FIGURE_PROVENANCE.md"))
    if doc is None:
        fail.append("docs/FIGURE_PROVENANCE.md is not in the deposit")
    else:
        have = {d["label"] for d in derived}
        documented = set(doc)
        ghosts = sorted(documented - have)
        absent = sorted(have - documented)
        if ghosts:
            fail.append("docs/FIGURE_PROVENANCE.md documents %d figure(s) the manuscript "
                        "does not typeset: %s" % (len(ghosts), ", ".join(ghosts)))
        if absent:
            fail.append("docs/FIGURE_PROVENANCE.md is missing %d figure(s) the manuscript "
                        "typesets: %s" % (len(absent), ", ".join(absent)))
        if not ghosts and not absent:
            print("  FIGURE_PROVENANCE.md      : covers exactly the %d figures" % len(have))
        # every file a row names must exist, unless the row says NOT DEPOSITED
        for label, row in sorted(doc.items()):
            for fn, exempt in row["files"]:
                if exempt or "*" in fn:
                    continue
                hits = [p for p in (os.path.join(ROOT, "paper", "figs", fn),
                                    os.path.join(ROOT, "paper", fn),
                                    os.path.join(ROOT, "src", fn),
                                    os.path.join(ROOT, "data", fn),
                                    os.path.join(ROOT, "docs", "img", fn),
                                    os.path.join(ROOT, fn.replace("/", os.sep)))
                        if os.path.isfile(p)]
                if not hits:
                    fail.append("FIGURE_PROVENANCE.md row %s names `%s`, which is in no "
                                "deposited directory, and the row does not say so (mark it "
                                "NOT DEPOSITED / DELETED / SUPERSEDED)" % (label, fn))

    # -- 3b'. every artefact check_figures.py guards must still be in the manuscript ----
    cf = os.path.join(ROOT, "src", "check_figures.py")
    if os.path.isfile(cf):
        txt = open(cf, encoding="utf-8", errors="replace").read()
        blk = re.search(r"ARTEFACTS\s*=\s*\[(.*?)\]", txt, re.S)
        if blk:
            guarded = re.findall(r"\"([^\"]+)\"", blk.group(1))
            stray = [g for g in guarded if g.endswith(".tex") and g not in reached_tex]
            if stray:
                fail.append("src/check_figures.py guards fragment(s) the manuscript no "
                            "longer \\input's: %s" % stray)

    # -- 3c. the orphan inventory -------------------------------------------------
    def load(rels):
        out = {}
        for r in rels:
            p = os.path.join(ROOT, r.replace("/", os.sep))
            if os.path.isfile(p):
                out[r] = open(p, encoding="utf-8", errors="replace").read()
        return out

    # src/frontier/ is deposited code too (the C3 frontier sweep, 2026-09-19).  A checker
    # that reads only the top level of src/ has a blind spot the size of a subdirectory,
    # and blind spots are how the documents this file exists to fix went stale.  Measured
    # when it was added: it moves nothing today -- no data/*.json and no paper/figs
    # artefact is named by any frontier script -- so the tiers below are unchanged and
    # the counts quoted in FILE_INDEX.md still hold.
    code = load(["src/" + n for n in rel_files("src", {".py"})]
                + ["src/frontier/" + n for n in rel_files("src/frontier", {".py"})]
                + ["Makefile", ".github/workflows/ci.yml",
                   "notebooks/00_Reproduce_Everything.ipynb",
                   "notebooks/Spectral_Heron.ipynb",
                   "notebooks/HW_LUCJ_N2_Heron_READY.ipynb"])
    docs = load(["README.md", "CITATION.cff"]
                + ["docs/" + n for n in rel_files("docs", {".md"})])

    # Only artefacts that are SUPPOSED to be produced by something are scanned: data
    # tables, datasets, compiled figure bodies, README thumbnails, and the .tex
    # fragments that declare themselves machine-written.  Hand-written LaTeX (the float
    # wrappers, the captions, the section files) has no producer by design and listing
    # it would bury the real orphans.
    GEN_MARK = re.compile(r"AUTO-GENERATED|Emitted by|Generated by|no hand-typed numbers",
                          re.I)
    scanned = []
    for sub, exts in (("paper/figs", {".dat", ".pdf", ".json", ".tex"}), ("data", {".json"}),
                      ("docs/img", {".png"})):
        for n in rel_files(sub, exts):
            rel = sub + "/" + n
            if n.endswith(".tex"):
                head = open(os.path.join(ROOT, sub.replace("/", os.sep), n),
                            encoding="utf-8", errors="replace").read(1200)
                if not GEN_MARK.search(head):
                    continue          # hand-written source: not an artefact
            scanned.append((rel, n))

    # The tiers use CODE only.  Documenting an orphan does not give it a producer, and
    # if a doc mention counted the orphan list would empty itself the moment it was
    # written down -- which is the opposite of the point.
    ORDER = ["orphan", "off-paper", "unsourced"]
    tiers = {t: [] for t in ORDER}
    for rel, n in scanned:
        c, d = named_by(n, code), named_by(n, docs)
        hit = rel in reached
        if hit and not c:
            tiers["unsourced"].append((rel, ["the manuscript reaches it; no code names it"]))
        elif hit:
            continue
        elif not c:
            tiers["orphan"].append((rel, ["documented in " + ", ".join(sorted(d))] if d
                                    else ["UNDOCUMENTED"]))
        else:
            tag = "checkers only" if c <= CHECKERS else ""
            tiers["off-paper"].append((rel, sorted(c) + ([tag] if tag else [])))

    print("\nORPHAN INVENTORY  (%d produced artefacts scanned; hand-written LaTeX excluded)"
          % len(scanned))
    print("  tier 1  ORPHAN     no deposited code names it and the manuscript does not use"
          " it : %d" % len(tiers["orphan"]))
    print("  tier 2  OFF-PAPER  code names it, the manuscript does not use it             "
          "      : %d" % len(tiers["off-paper"]))
    print("  tier 3  UNSOURCED  the manuscript uses it, no deposited code names it        "
          "      : %d" % len(tiers["unsourced"]))
    undoc = [r for r, w in tiers["orphan"] if w == ["UNDOCUMENTED"]]
    if undoc:
        print("           of which UNDOCUMENTED (no .md mentions them either)          "
              "          : %d" % len(undoc))
    if verbose:
        for t in ORDER:
            for rel, who in sorted(tiers[t]):
                print("    [%-9s] %-46s %s" % (t, rel, ", ".join(who) if who else "--"))

    if as_md:
        print("\n<!-- derived by src/check_provenance.py; do not hand-edit the label column -->")
        print("| # | Body ref | Assembled in | Figure body | pgfplots tables |")
        print("|---|----------|--------------|-------------|-----------------|")
        for d in derived:
            body = ", ".join("`%s`" % x for x in (d["graphics"] or d["files"])) or "--"
            tb = ", ".join("`%s`" % x for x in d["tables"]) or "--"
            print("| %d | %s | `paper/%s` | %s | %s |"
                  % (d["n"], d["label"], d["wrapper"], body, tb))

    if fail:
        print("\nFAILURES:")
        for f in fail:
            print("  * %s" % f)
        print("\nRESULT: FAIL -- %d disagreement(s) between the deposit's documentation and "
              "the manuscript in the tree." % len(fail))
        return 1
    print("\nRESULT: PASS -- the documentation describes the manuscript that is in the tree.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
