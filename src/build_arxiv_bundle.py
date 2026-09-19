# -*- coding: utf-8 -*-
"""Build the arXiv **v3** source package from the working tree, and prove it.

    python src/build_arxiv_bundle.py                 # build + audit + clean-room compile
    python src/build_arxiv_bundle.py --no-compile    # closure + audit only (fast, ~1 s)
    python src/build_arxiv_bundle.py --selftest      # run the five negative controls
    python src/build_arxiv_bundle.py --list          # print the derived file list
    python src/build_arxiv_bundle.py --passes 3 --out build/arxiv-v3.tar.gz

WHY THIS FILE EXISTS
--------------------
``paper/arxiv-submission.tar.gz`` is the FROZEN arXiv **v2** package.  It carries the
files this manuscript no longer has (``introduction/method/results/discussion/conclusion/
hardware/theorem_b2``), five withdrawn figures, the two fabricated rows of
``figs/sqw_edges.dat`` and the superseded ``19184``/``824504``.  It is kept as a
historical record and **must never be uploaded again**.  This script builds the v3
package instead, into ``build/`` (a scratch directory, gitignored: the package is a
build product, not a deposited artefact -- it is one command away at any time).

THE FILE LIST IS NEVER WRITTEN BY HAND
--------------------------------------
It is derived by recursive expansion of ``\\input`` / ``\\include`` /
``\\includegraphics`` / ``table[...]{...}`` from ``paper/main.tex``.  A previous
attempt wrote the pattern as ``table\\[[^]]*\\]\\{`` -- it REQUIRED the optional
bracket group -- and so silently dropped the ``.dat`` files this manuscript loads as
``table {n3_pub_frac_lo.dat}`` (no brackets, and the keyword sits on its own line,
separated from its ``\\addplot``).  Nothing noticed until a clean-room rehearsal died
inside pgfplots.  So the parser here is written as a scanner, not as a line regex: it
reads the keyword, then an OPTIONAL bracket group, then the brace group, across
newlines.  The exact damage is not quoted from memory -- ``--selftest`` control 1
re-runs the broken derivation and prints it: **6 of the 10 .dat files lost**, 60 files
instead of 66.

AND THE LIST IS CHECKED THREE INDEPENDENT WAYS
----------------------------------------------
A check that can pass for the wrong reason is worth nothing, so the closure is derived
once and confirmed twice more, by mechanisms that share no code with it:

  1. the scanner above (keyword-driven);
  2. a filename-shaped-literal sweep: every ``{...}`` group in every shipped .tex whose
     content looks like ``name.dat|csv|tsv|txt|pdf|png|jpg|eps`` must be in the closure
     (catches a keyword form the scanner does not know about);
  3. TeX's own testimony: the clean-room compile runs with ``-recorder``, and every
     INPUT line of the resulting ``main.fls`` that points inside the extraction
     directory must be in the closure, and vice versa (catches everything the two
     static passes can miss, including files pulled in by macro expansion).

and the equivalence test compares the clean-room PDF against BOTH a same-pass-count
recompilation of the full tree AND the committed ``paper/main.pdf``.

``--selftest`` mutates a copy of the package five ways and asserts that the five checks
FAIL.  A PASS with no negative control is not evidence.
"""
from __future__ import print_function

import argparse
import gzip
import hashlib
import io
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.join(ROOT, "paper")
MAIN = "main.tex"

# ----------------------------------------------------------------------------------
# Search paths TeX itself will use inside the extracted package.  Both come from
# paper/main.tex and are re-read from it below, never assumed:
#     \graphicspath{{figs/}}                  -> \includegraphics
#     \pgfplotsset{table/search path={figs}}  -> table{...}
SEARCH_DIRS_DEFAULT = ["", "figs"]

EXT_BY_KIND = {
    "input": ["", ".tex"],
    "graphic": ["", ".pdf", ".png", ".jpg", ".jpeg", ".eps"],
    "table": ["", ".dat", ".txt", ".csv", ".tsv"],
}

# keyword -> (kind, how many brace groups the macro takes before the useful one)
INPUT_MACROS = {
    "input": "input",
    "include": "input",
    "subfile": "input",
    "InputIfFileExists": "input",
    "lstinputlisting": "input",
    "verbatiminput": "input",
}
GRAPHIC_MACROS = {"includegraphics": "graphic"}
TABLE_KEYWORDS = {"table": "table", "pgfplotstableread": "table"}

DATA_EXT_RE = re.compile(r"^[\w./+-]+\.(?:dat|csv|tsv|txt|pdf|png|jpg|jpeg|eps)$", re.I)

# ---- what must NOT be in the package ---------------------------------------------
FORBIDDEN_BASENAMES = [
    # v2 body files, superseded by sec_*.tex
    "introduction.tex", "method.tex", "results.tex", "discussion.tex",
    "conclusion.tex", "hardware.tex", "theorem_b2.tex", "resource.tex",
    "table_molecules.tex", "fig_akw_sampled_native.tex",
    # the five withdrawn figures
    "fig_amort_native.tex", "fig_gflownet.pdf", "fig_benchmark.pdf",
    "fig_molecular_gallery.pdf", "fig_molecular_suite.pdf",
    # editor / build by-products
    "mainNotes.bib", "main.pdf", "arxiv-submission.tar.gz",
    "_method_standalone.tex", "_method_standalone.pdf",
]
FORBIDDEN_EXTS = [
    ".aux", ".log", ".out", ".toc", ".fls", ".fdb_latexmk", ".bbl", ".blg",
    ".synctex", ".nav", ".snm", ".vrb", ".bib", ".py", ".pyc", ".ipynb",
    ".gz", ".zip", ".tgz", ".bak", ".orig", ".rej",
]
FORBIDDEN_PATH_PARTS = ["_superseded", "__pycache__", ".git", "rebuild", "ckpt",
                        "notebooks", "data"]

# the two rows that were fabricated in figs/sqw_edges.dat and removed on 2026-09-18
FABRICATED_ROWS = [
    re.compile(r"^\s*0\.0000\s+4\.969\s+5\.000\s+5\.000\s*$", re.M),
    re.compile(r"^\s*2\.0000\s+4\.969\s+5\.000\s+5\.000\s*$", re.M),
]
# the superseded |S| values of Fig. 12 (correct: 19183 / 824503)
SUPERSEDED_NUMBERS = [
    ("19184", re.compile(r"(?<![\d.])19184(?![\d])")),
    ("824504", re.compile(r"(?<![\d.])824504(?![\d])")),
]
# A Windows drive letter, or a unix build path.  The negative look-behind is not
# cosmetic: without it `https://` matches as drive `s:` and the check drowns in 57
# false positives from the bibliography's DOIs -- a check that cries wolf gets
# switched off, which is how a real absolute path would then slip through.
ABSOLUTE_PATH_RE = re.compile(
    r"(?<![A-Za-z0-9])(?:[A-Za-z]:[\\/](?![/])|/home/|/Users/|/mnt/[a-z]/|/tmp/|/w/)")


# ==================================================================================
# 1.  TeX scanning
# ==================================================================================
def strip_tex_comments(text):
    """Remove TeX comments, honouring every line terminator TeX honours.

    MiKTeX/pdfTeX break lines on CR, LF *and* CRLF.  A lone CR inside a comment
    therefore ENDS the comment and the rest of the line is typeset -- which is exactly
    the defect that printed a build note on page 43 of this manuscript.  Normalising
    CR -> LF before stripping makes this parser see what TeX sees; reading the file
    through a newline-translating text handle would instead make the parser see what
    the author INTENDED, and that is how the defect survived a scan once already.
    """
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    out = []
    for line in text.split("\n"):
        i, cut = 0, None
        while i < len(line):
            c = line[i]
            if c == "\\":          # \% is a literal percent; \\ is a break, then rescan
                i += 2
                continue
            if c == "%":
                cut = i
                break
            i += 1
        out.append(line if cut is None else line[:cut])
    return "\n".join(out)


def _skip_space(s, i):
    while i < len(s) and s[i] in " \t\r\n":
        i += 1
    return i


def _read_group(s, i, opener, closer):
    """s[i] must be `opener`.  Returns (content, index_after_closer) or (None, i)."""
    if i >= len(s) or s[i] != opener:
        return None, i
    depth, j = 0, i
    while j < len(s):
        c = s[j]
        if c == "\\":
            j += 2
            continue
        if c == opener:
            depth += 1
        elif c == closer:
            depth -= 1
            if depth == 0:
                return s[i + 1:j], j + 1
        j += 1
    return None, i


def scan_references(text):
    """Return [(kind, target, char_offset)] for every file reference in `text`.

    `text` must already be comment-stripped.  The scanner walks the string; it does
    NOT work line by line, because `\\addplot[...]` and its `table {file.dat}` sit on
    different lines in figs/fig_thm1iii_violation_native.tex.
    """
    refs = []
    n = len(text)
    i = 0
    word = re.compile(r"[A-Za-z]+")
    while i < n:
        c = text[i]
        if c == "\\":
            m = word.match(text, i + 1)
            if not m:
                i += 2
                continue
            name = m.group(0)
            j = m.end()
            if name in INPUT_MACROS or name in GRAPHIC_MACROS or name in TABLE_KEYWORDS:
                kind = (INPUT_MACROS.get(name) or GRAPHIC_MACROS.get(name)
                        or TABLE_KEYWORDS.get(name))
                got = _consume_args(text, j)
                if got is not None:
                    refs.append((kind, got[0], i))
                    i = got[1]
                    continue
            i = j
            continue
        if c.isalpha():
            m = word.match(text, i)
            name = m.group(0)
            j = m.end()
            # bare keyword, e.g. pgfplots' `table {file.dat}` after \addplot[...]
            if name in TABLE_KEYWORDS and (i == 0 or not (text[i - 1].isalpha()
                                                          or text[i - 1] == "\\")):
                got = _consume_args(text, j)
                if got is not None:
                    refs.append((TABLE_KEYWORDS[name], got[0], i))
                    i = got[1]
                    continue
            i = j
            continue
        i += 1
    return refs


def _consume_args(text, j):
    """After a keyword: skip whitespace, an OPTIONAL [..] group, whitespace, then read
    the mandatory {..} group.  Returns (content, index_after) or None.

    The optionality of the bracket group is the whole point: requiring it is what lost
    ten .dat files last time.
    """
    j = _skip_space(text, j)
    if j < len(text) and text[j] == "[":
        opt, j2 = _read_group(text, j, "[", "]")
        if opt is None:
            return None
        j = _skip_space(text, j2)
    if j < len(text) and text[j] == "{":
        arg, j2 = _read_group(text, j, "{", "}")
        if arg is None:
            return None
        return arg.strip(), j2
    return None


def read_search_dirs(main_text):
    """Read \\graphicspath and pgfplots' table/search path OUT OF main.tex."""
    dirs = [""]
    for m in re.finditer(r"\\graphicspath\s*\{(.*?)\}\s*(?:\n|%|$)", main_text, re.S):
        for d in re.findall(r"\{([^{}]*)\}", m.group(1)):
            dirs.append(d.rstrip("/"))
    for m in re.finditer(r"table/search\s+path\s*=\s*\{?([\w./-]+)\}?", main_text):
        dirs.append(m.group(1).rstrip("/"))
    seen, out = set(), []
    for d in dirs:
        if d not in seen:
            seen.add(d)
            out.append(d)
    return out


# ==================================================================================
# 2.  the closure
# ==================================================================================
class Closure(object):
    def __init__(self):
        self.files = []          # relative-to-paper/ POSIX paths, in discovery order
        self.seen = set()
        self.unresolved = []     # (kind, target, referencing file)
        self.edges = []          # (referencing file, kind, resolved file)

    def add(self, rel):
        if rel not in self.seen:
            self.seen.add(rel)
            self.files.append(rel)
            return True
        return False


def _resolve(target, kind, search_dirs, from_rel):
    cands = []
    for d in search_dirs:
        base = os.path.join(d, target) if d else target
        for ext in EXT_BY_KIND[kind]:
            cands.append(os.path.normpath(base + ext).replace("\\", "/"))
    # also, last resort, next to the referencing file (MiKTeX will do this too)
    here = os.path.dirname(from_rel)
    if here:
        for ext in EXT_BY_KIND[kind]:
            cands.append(os.path.normpath(os.path.join(here, target + ext)).replace("\\", "/"))
    for c in cands:
        if c.startswith(".."):
            continue
        if os.path.isfile(os.path.join(PAPER, c.replace("/", os.sep))):
            return c
    return None


def build_closure(paper=PAPER, main=MAIN):
    cl = Closure()
    main_text = open(os.path.join(paper, main), "rb").read().decode("utf-8", "replace")
    search_dirs = read_search_dirs(strip_tex_comments(main_text))
    cl.add(main)
    queue = [main]
    while queue:
        rel = queue.pop(0)
        path = os.path.join(paper, rel.replace("/", os.sep))
        if os.path.splitext(rel)[1].lower() != ".tex":
            continue
        raw = open(path, "rb").read().decode("utf-8", "replace")
        body = strip_tex_comments(raw)
        for kind, target, _off in scan_references(body):
            if kind == "table" and not DATA_EXT_RE.match(target) and "." not in target:
                continue                       # `table` in a tabular spec, not a file
            got = _resolve(target, kind, search_dirs, rel)
            if got is None:
                if kind == "table" and not DATA_EXT_RE.match(target):
                    continue                   # e.g. \begin{table}{...}: not a file
                cl.unresolved.append((kind, target, rel))
                continue
            cl.edges.append((rel, kind, got))
            if cl.add(got):
                queue.append(got)
    return cl, search_dirs


def literal_sweep(closure, paper=PAPER):
    """Second, INDEPENDENT detector: any brace group in any shipped .tex whose content
    is filename-shaped must already be in the closure."""
    missed = []
    for rel in closure.files:
        if os.path.splitext(rel)[1].lower() != ".tex":
            continue
        body = strip_tex_comments(
            open(os.path.join(paper, rel.replace("/", os.sep)), "rb")
            .read().decode("utf-8", "replace"))
        for grp in re.findall(r"\{([^{}]{1,120})\}", body):
            t = grp.strip()
            if not DATA_EXT_RE.match(t):
                continue
            norm = os.path.normpath(t).replace("\\", "/")
            hit = norm in closure.seen or any(
                os.path.normpath(os.path.join(d, t)).replace("\\", "/") in closure.seen
                for d in ("", "figs"))
            if not hit:
                missed.append((rel, t))
    return missed


# ==================================================================================
# 3.  packaging (deterministic: two runs give byte-identical tarballs)
# ==================================================================================
def make_tarball(closure, out_path, paper=PAPER):
    members = sorted(closure.files)
    buf = io.BytesIO()
    tf = tarfile.open(fileobj=buf, mode="w", format=tarfile.USTAR_FORMAT)
    for rel in members:
        src = os.path.join(paper, rel.replace("/", os.sep))
        data = open(src, "rb").read()
        ti = tarfile.TarInfo(rel)
        ti.size = len(data)
        ti.mtime = 0
        ti.mode = 0o644
        ti.uid = ti.gid = 0
        ti.uname = ti.gname = ""
        ti.type = tarfile.REGTYPE
        tf.addfile(ti, io.BytesIO(data))
    tf.close()
    raw = buf.getvalue()
    gz = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=gz, mtime=0) as g:
        g.write(raw)
    blob = gz.getvalue()
    d = os.path.dirname(out_path)
    if d and not os.path.isdir(d):
        os.makedirs(d)
    open(out_path, "wb").write(blob)
    return blob, members


# ==================================================================================
# 4.  the audit -- run against the TARBALL, not against the tree
# ==================================================================================
def audit_bundle(tar_path, paper=PAPER):
    """Returns (failures, info).  Reads the bytes that would actually be uploaded."""
    fails, info = [], {}
    with tarfile.open(tar_path, "r:gz") as tf:
        members = [m for m in tf.getmembers()]
        names = [m.name for m in members if m.isfile()]
        blobs = {}
        for m in members:
            if m.isfile():
                blobs[m.name] = tf.extractfile(m).read()
    info["n_files"] = len(names)
    info["bytes_uncompressed"] = sum(len(b) for b in blobs.values())
    info["bytes_gz"] = os.path.getsize(tar_path)

    for m_name in names:
        base = os.path.basename(m_name)
        low = base.lower()
        if base in FORBIDDEN_BASENAMES:
            fails.append("RETIRED FILE SHIPPED: %s" % m_name)
        for e in FORBIDDEN_EXTS:
            if low.endswith(e):
                fails.append("FORBIDDEN EXTENSION %s: %s" % (e, m_name))
        parts = m_name.replace("\\", "/").split("/")
        for p in FORBIDDEN_PATH_PARTS:
            if p in parts[:-1]:
                fails.append("FORBIDDEN DIRECTORY %r in path: %s" % (p, m_name))
        if m_name.startswith("/") or re.match(r"^[A-Za-z]:", m_name) or ".." in parts:
            fails.append("UNSAFE MEMBER PATH: %s" % m_name)
        if any(ord(ch) > 126 for ch in m_name) or " " in m_name:
            fails.append("NON-ASCII OR SPACE IN FILENAME: %s" % m_name)

    # fabricated rows + superseded numbers, in any shipped file
    for name, blob in sorted(blobs.items()):
        try:
            txt = blob.decode("utf-8")
        except UnicodeDecodeError:
            continue                                   # binary (the figure PDFs)
        for k, rx in enumerate(FABRICATED_ROWS):
            if rx.search(txt):
                fails.append("FABRICATED sqw_edges ROW #%d present in %s" % (k + 1, name))
        for label, rx in SUPERSEDED_NUMBERS:
            if rx.search(txt):
                for mm in rx.finditer(txt):
                    line = txt.count("\n", 0, mm.start()) + 1
                    fails.append("SUPERSEDED NUMBER %s in %s line %d" % (label, name, line))
        if name.lower().endswith(".tex"):
            for mm in ABSOLUTE_PATH_RE.finditer(strip_tex_comments(txt)):
                fails.append("ABSOLUTE PATH %r in %s" % (mm.group(0), name))

    # the two withdrawn floats must not be pulled in even by a live \input
    main_body = strip_tex_comments(blobs["main.tex"].decode("utf-8", "replace"))
    for tag in ("fig_amort", "fig_gflow", "fig_gflownet", "fig_benchmark",
                "fig_molecular_gallery", "fig_molecular_suite"):
        if tag in main_body:
            fails.append("WITHDRAWN FIGURE %s is referenced by live (uncommented) "
                         "code in main.tex" % tag)

    # stray CR: the page-43 defect class, checked at byte level in what we ship
    for name, blob in sorted(blobs.items()):
        if not name.lower().endswith((".tex", ".dat")):
            continue
        i = 0
        while True:
            i = blob.find(b"\r", i)
            if i < 0:
                break
            if blob[i + 1:i + 2] != b"\n":
                fails.append("STRAY CR (0x0D) at byte %d of %s -- this is the defect "
                             "class that typeset a build note on page 43" % (i, name))
            i += 1

    # every shipped byte must be the working tree's byte
    for name, blob in sorted(blobs.items()):
        src = os.path.join(paper, name.replace("/", os.sep))
        if not os.path.isfile(src):
            fails.append("MEMBER HAS NO SOURCE IN THE TREE: %s" % name)
        elif open(src, "rb").read() != blob:
            fails.append("MEMBER DIFFERS FROM THE TREE: %s" % name)
    return fails, info


# ==================================================================================
# 5.  arXiv's own requirements
# ==================================================================================
def arxiv_requirements(tar_path):
    """Returns [(id, verdict, detail)] -- verdict in {'PASS','FAIL','N/A','INFO'}.

    Every rule below was re-read from info.arxiv.org on 2026-09-19, NOT recalled.  That
    mattered: arXiv retired AutoTeX in April 2025 (Submission System 1.5) and the page
    now says, verbatim, "You should not use \\pdfoutput to change the output format."
    The processor is chosen by the submitter at upload instead.  This repository's own
    notes still carried the pre-2025 rule "\\pdfoutput=1 within the first five lines",
    so that rule is checked here in its CURRENT form -- inverted.
    """
    out = []
    with tarfile.open(tar_path, "r:gz") as tf:
        names = [m.name for m in tf.getmembers() if m.isfile()]
        main = tf.extractfile("main.tex").read().decode("utf-8", "replace")
    main_lf = main.replace("\r\n", "\n").replace("\r", "\n")
    lines = main_lf.split("\n")
    body = strip_tex_comments(main_lf)

    raster = [n for n in names if n.lower().endswith((".pdf", ".png", ".jpg", ".jpeg"))]
    eps = [n for n in names if n.lower().endswith((".eps", ".ps"))]

    # --- 1. the processor, and the rule that replaced \pdfoutput ---------------------
    po = re.search(r"\\pdfoutput\s*=", body)
    out.append(("no \\pdfoutput (retired rule, now discouraged)",
                "PASS" if po is None else "FAIL",
                'info.arxiv.org/help/submit_tex: "You should not use \\pdfoutput to '
                'change the output format." %s'
                % ("absent from main.tex -- correct" if po is None
                   else "main.tex sets it; remove it")))
    out.append(("submission processor must be pdflatex",
                "ACTION" if raster else "N/A",
                "%d PDF figures are embedded, so 'LaTeX with PDFLaTeX' must be chosen "
                "in the upload form (00README.json compiler: pdflatex). Measured "
                "consequence of the wrong choice is printed by the DVI probe below."
                % len(raster)))
    out.append(("figure formats are unified (no PS/EPS mixed with PDF)",
                "PASS" if not (eps and raster) else "FAIL",
                'arXiv "requires that you use a unified figure format"; %d eps/ps, '
                "%d pdf/png/jpg" % (len(eps), len(raster))))

    # --- 2. paper size --------------------------------------------------------------
    m = re.search(r"\\documentclass\s*(\[[^\]]*\])?\s*\{([^}]*)\}", body)
    opts = re.sub(r"\s+", "", (m.group(1) or "")) if m else ""
    cls = m.group(2) if m else "?"
    explicit = bool(re.search(r"(a4paper|letterpaper|a5paper|legalpaper)", opts)) or \
        bool(re.search(r"\\(special\s*\{papersize|pdfpagewidth|geometry)", body))
    out.append(("explicit paper size",
                "N/A" if not explicit else "PASS",
                "no paper-size rule on info.arxiv.org for the 1.5 system (it was a "
                "dvips-era rule); %s{%s} sets the size itself -- the MediaBox measured "
                "below is the evidence" % (cls, opts)))
    out.append(("exactly one \\documentclass in the package",
                "PASS" if sum(
                    1 for n in names if n.endswith(".tex") and "\\documentclass" in
                    strip_tex_comments(open_member(tar_path, n))) == 1 else "FAIL",
                "the top-level file must be unambiguous"))

    # --- 3. \today ------------------------------------------------------------------
    has_today = "\\today" in body
    out.append(("no \\today", "WARN" if has_today else "PASS",
                'info.arxiv.org/help/submit_tex: "arXiv recommends against using the '
                '\\today macro in the standard \\date field." %s'
                % ("main.tex line %d has it; REVTeX prints no date at all if the "
                   "\\date line is deleted (\\def\\@date{} at revtex4-2.cls:2521)"
                   % next((i + 1 for i, l in enumerate(lines) if "\\today" in l), 0)
                   if has_today else "absent")))

    # --- 4. bibliography ------------------------------------------------------------
    bibs = [n for n in names if n.endswith(".bib")]
    bbls = [n for n in names if n.endswith(".bbl")]
    needs_bbl = re.search(r"\\(?:bibliography|addbibresource)\s*\{", body) is not None
    out.append(("bibliography is self-contained",
                "PASS" if (not needs_bbl and not bibs and not bbls) else "FAIL",
                "inline thebibliography; %d .bib, %d .bbl shipped, \\bibliography{} "
                "%s -- arXiv blocks a submission that needs a .bbl and has neither"
                % (len(bibs), len(bbls), "used" if needs_bbl else "unused")))

    # --- 5. no auxiliary or output files, no hidden files ---------------------------
    auxbad = [n for n in names
              if os.path.splitext(n)[1].lower() in (".aux", ".log", ".toc", ".lot",
                                                    ".lof", ".dvi", ".pdf", ".ps")
              and os.path.splitext(os.path.basename(n))[0] == "main"]
    out.append(("no auxiliary/output file of main.tex",
                "PASS" if not auxbad else "FAIL",
                'arXiv: "do not include any associated auxiliary file or intermediate '
                'or resulting output file"; %s' % (auxbad or "none present")))
    hidden = [n for n in names if any(p.startswith(".") for p in n.split("/"))]
    out.append(("no hidden files", "PASS" if not hidden else "FAIL",
                'arXiv: "Hidden files will be deleted upon announcement"; %s'
                % (hidden or "none present")))

    # --- 6. embedded JavaScript in the figure PDFs is an automatic rejection ---------
    js = []
    with tarfile.open(tar_path, "r:gz") as tf:
        for n in names:
            if not n.lower().endswith(".pdf"):
                continue
            blob = tf.extractfile(n).read()
            for tok in (b"/JavaScript", b"/JS ", b"/RichMedia", b"/Movie", b"/Launch"):
                if tok in blob:
                    js.append("%s carries %s" % (n, tok.decode()))
    out.append(("no embedded JavaScript / rich media in the figure PDFs",
                "PASS" if not js else "FAIL",
                'arXiv: "Submissions with embedded JavaScript are automatically '
                'rejected"; %s' % (js or "%d figure PDFs scanned, clean" % len(
                    [n for n in names if n.lower().endswith(".pdf")]))))

    size_mb = os.path.getsize(tar_path) / 1048576.0
    out.append(("package size", "PASS" if size_mb < 50 else "FAIL",
                "%.2f MB gzipped" % size_mb))
    return out


def open_member(tar_path, name):
    with tarfile.open(tar_path, "r:gz") as tf:
        return tf.extractfile(name).read().decode("utf-8", "replace")


# ==================================================================================
# 6.  clean-room compile and PDF equivalence
# ==================================================================================
def run_latex(workdir, passes, recorder=True):
    cmd = ["pdflatex", "-interaction=nonstopmode", "-file-line-error"]
    if recorder:
        cmd.append("-recorder")
    cmd.append("main.tex")
    rc = None
    for _ in range(passes):
        p = subprocess.Popen(cmd, cwd=workdir, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        p.communicate()
        rc = p.returncode
    log = ""
    lp = os.path.join(workdir, "main.log")
    if os.path.isfile(lp):
        log = open(lp, "rb").read().decode("utf-8", "replace")
    return rc, log


# pdfTeX prints an error in TWO shapes and this script asks for -file-line-error, which
# selects the SECOND one.  A detector that only knows `! ...` therefore reports "0
# errors" on a run that failed -- measured here on 2026-09-19: the DVI probe below
# raises six "LaTeX Error: Cannot determine size of graphic" and the `! `-only detector
# scored it 0.  Both shapes are matched, and src/build_arxiv_bundle.py --selftest holds
# a positive control that makes a real error and demands this function see it.
ERROR_SHAPES = (
    re.compile(r"^!\s.*$", re.M),                        # classic:  ! LaTeX Error: ...
    re.compile(r"^[^\s:]+:\d+:\s.*$", re.M),             # -file-line-error: f.tex:12: ...
    re.compile(r"^.*(?:Emergency stop|Fatal error occurred).*$", re.M),
)


def log_stats(log):
    pages = None
    m = re.search(r"Output written on .*?\((\d+) pages?", log)
    if m:
        pages = int(m.group(1))
    errors = []
    for rx in ERROR_SHAPES:
        errors.extend(rx.findall(log))
    errors = sorted(set(errors))
    undef_ref = len(re.findall(r"Reference `[^']*' on page \d+ undefined", log))
    undef_cit = len(re.findall(r"Citation `[^']*' on page \d+ undefined", log))
    missing = len(re.findall(r"LaTeX Warning: File `[^']*' not found", log)) + \
        len(re.findall(r"I could not locate the file", log))
    return dict(pages=pages, errors=errors, n_errors=len(errors),
                undef_ref=undef_ref, undef_cit=undef_cit, missing_files=missing)


def pdf_text_pages(path):
    try:
        import fitz                                    # PyMuPDF
        doc = fitz.open(path)
        return [re.sub(r"\s+", " ", doc[i].get_text("text")).strip()
                for i in range(doc.page_count)]
    except Exception:
        from pypdf import PdfReader
        r = PdfReader(path)
        return [re.sub(r"\s+", " ", (p.extract_text() or "")).strip() for p in r.pages]


def pdf_media_boxes(path):
    try:
        import fitz
        doc = fitz.open(path)
        return sorted({(round(doc[i].rect.width, 1), round(doc[i].rect.height, 1))
                       for i in range(doc.page_count)})
    except Exception:
        from pypdf import PdfReader
        r = PdfReader(path)
        return sorted({(round(float(p.mediabox.width), 1),
                        round(float(p.mediabox.height), 1)) for p in r.pages})


def compare_pdfs(a, b, label_a, label_b):
    pa, pb = pdf_text_pages(a), pdf_text_pages(b)
    sa = hashlib.sha1("\n".join(pa).encode("utf-8")).hexdigest()[:16]
    sb = hashlib.sha1("\n".join(pb).encode("utf-8")).hexdigest()[:16]
    diffs = [i + 1 for i in range(min(len(pa), len(pb))) if pa[i] != pb[i]]
    res = dict(pages_a=len(pa), pages_b=len(pb), sha_a=sa, sha_b=sb,
               chars_a=sum(len(x) for x in pa), chars_b=sum(len(x) for x in pb),
               diff_pages=diffs, equal=(len(pa) == len(pb) and not diffs))
    print("  %-34s %3d pp  %9d chars  text-sha1 %s" % (label_a, res["pages_a"],
                                                       res["chars_a"], sa))
    print("  %-34s %3d pp  %9d chars  text-sha1 %s" % (label_b, res["pages_b"],
                                                       res["chars_b"], sb))
    if res["equal"]:
        print("  -> IDENTICAL text flow on every page")
    else:
        print("  -> DIFFERENT: %d page(s) differ%s" %
              (len(diffs), (": " + str(diffs[:12])) if diffs else ""))
    return res


def fls_cross_check(workdir, closure_files):
    """TeX's own testimony about what it opened."""
    fls = os.path.join(workdir, "main.fls")
    if not os.path.isfile(fls):
        return None
    opened = set()
    base = os.path.normcase(os.path.abspath(workdir))
    repo = os.path.normcase(os.path.abspath(ROOT))
    leaked = []
    n_input = 0
    for line in open(fls, "r", errors="replace"):
        if not line.startswith("INPUT "):
            continue
        n_input += 1
        p = line[6:].strip()
        ap = os.path.normcase(os.path.abspath(os.path.join(workdir, p)))
        if not ap.startswith(base + os.sep):
            # Not in the clean room.  Almost always a texmf file -- but if TEXINPUTS
            # (or a stray absolute path) reached back into the repository, the compile
            # would succeed on a file the PACKAGE DOES NOT CONTAIN, and this whole
            # check would pass for the wrong reason.  Name it instead of skipping it.
            if ap.startswith(repo + os.sep):
                leaked.append(os.path.relpath(ap, repo).replace("\\", "/"))
            continue
        rel = os.path.relpath(ap, base).replace("\\", "/")
        if os.path.splitext(rel)[1].lower() in (".aux", ".out", ".log", ".fls", ".toc",
                                                ".pdf") and rel.startswith("main"):
            continue
        opened.add(rel)
    have = set(os.path.normcase(f) for f in closure_files)
    opened_n = set(os.path.normcase(f) for f in opened)
    return dict(opened=sorted(opened), n_input=n_input, leaked=sorted(set(leaked)),
                tex_opened_not_shipped=sorted(opened_n - have),
                shipped_never_opened=sorted(have - opened_n))


# ==================================================================================
# 7.  driver
# ==================================================================================
def stale_sources(closure, paper=PAPER):
    """[(file, seconds_newer)] for every packaged source newer than paper/main.pdf."""
    pdf = os.path.join(paper, "main.pdf")
    if not os.path.isfile(pdf):
        return []
    t = os.path.getmtime(pdf)
    out = []
    for rel in closure.files:
        f = os.path.join(paper, rel.replace("/", os.sep))
        if os.path.isfile(f) and os.path.getmtime(f) > t + 1:
            out.append((rel, os.path.getmtime(f) - t))
    return sorted(out, key=lambda x: -x[1])


def leftovers(closure, paper=PAPER):
    """Every file in paper/ that the closure does NOT reach.  Printed so that a human
    can see what is being left out and say whether that is intended -- the opposite
    failure mode from the ten lost .dat files, and just as invisible otherwise."""
    out = []
    for r, dirs, fs in os.walk(paper):
        dirs[:] = [d for d in dirs if d not in ("_superseded", "__pycache__")]
        for f in fs:
            rel = os.path.relpath(os.path.join(r, f), paper).replace(os.sep, "/")
            if rel in closure.seen:
                continue
            if os.path.splitext(rel)[1].lower() in (".aux", ".log", ".out", ".fls",
                                                    ".toc", ".synctex", ".fdb_latexmk"):
                continue
            out.append(rel)
    return sorted(out)


def build(out_path, verbose=True):
    cl, search_dirs = build_closure()
    if verbose:
        print("CLOSURE from paper/%s  (search dirs %r)" % (MAIN, search_dirs))
        print("  files reached : %d" % len(cl.files))
        by_ext = {}
        for f in cl.files:
            by_ext[os.path.splitext(f)[1] or "(none)"] = \
                by_ext.get(os.path.splitext(f)[1] or "(none)", 0) + 1
        print("  by extension  : %s" % ", ".join("%s %d" % (k, v)
                                                 for k, v in sorted(by_ext.items())))
    if cl.unresolved:
        print("  UNRESOLVED REFERENCES (the package would be incomplete):")
        for k, t, f in cl.unresolved:
            print("      %-9s %-40s referenced by %s" % (k, t, f))
    missed = literal_sweep(cl)
    if missed:
        print("  LITERAL SWEEP found filename-shaped targets NOT in the closure:")
        for f, t in missed:
            print("      %s  ->  %s" % (f, t))
    blob, members = make_tarball(cl, out_path)
    blob2, _ = make_tarball(cl, out_path + ".rebuild")
    os.remove(out_path + ".rebuild")
    if verbose:
        print("  wrote %s  (%d members, %.2f MB, sha1 %s)" %
              (os.path.relpath(out_path, ROOT), len(members),
               len(blob) / 1048576.0, hashlib.sha1(blob).hexdigest()[:16]))
        print("  deterministic rebuild byte-identical : %s" % (blob == blob2))
    return cl, out_path, (blob == blob2)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.join(ROOT, "build", "arxiv-v3.tar.gz"))
    ap.add_argument("--passes", type=int, default=3)
    ap.add_argument("--work", default=None,
                    help="clean-room directory (default: a fresh temp dir)")
    ap.add_argument("--no-compile", action="store_true")
    ap.add_argument("--keep", action="store_true", help="keep the clean-room directory")
    ap.add_argument("--selftest", action="store_true",
                    help="negative controls: mutate the package and demand FAIL")
    ap.add_argument("--list", action="store_true",
                    help="print the derived file list and exit")
    args = ap.parse_args(argv)

    if args.selftest:
        return selftest(args)

    if args.list:
        cl, sd = build_closure()
        print("# %d files, derived from paper/main.tex; search dirs %r" % (len(cl.files), sd))
        for f in sorted(cl.files):
            print(f)
        return 0

    print("=" * 86)
    print("BUILD THE arXiv v3 SOURCE PACKAGE")
    print("=" * 86)
    cl, out, deterministic = build(os.path.abspath(args.out))
    bad = list(cl.unresolved) + literal_sweep(cl)
    left = leftovers(cl)
    if left:
        print("  LEFT BEHIND in paper/ on purpose (%d files, none reachable from "
              "main.tex):" % len(left))
        for f in left:
            print("      %s" % f)
    print()

    print("-" * 86)
    print("AUDIT: what is in the package (read from the tarball, not from the tree)")
    print("-" * 86)
    fails, info = audit_bundle(out)
    print("  members %d   uncompressed %.2f MB   gzipped %.2f MB"
          % (info["n_files"], info["bytes_uncompressed"] / 1048576.0,
             info["bytes_gz"] / 1048576.0))
    if fails:
        for f in fails:
            print("  FAIL  %s" % f)
    else:
        print("  PASS  nothing retired, fabricated, superseded or auxiliary is in the "
              "package")
        print("        (checked: 15 retired basenames, 22 forbidden extensions, "
              "7 forbidden directories,")
        print("         2 fabricated sqw_edges rows, 19184/824504, stray CR bytes, "
              "absolute paths,")
        print("         and every member byte-compared against its source in paper/)")
    print()

    print("-" * 86)
    print("arXiv REQUIREMENTS")
    print("-" * 86)
    reqs = arxiv_requirements(out)
    req_fail = [r for r in reqs if r[1] == "FAIL"]
    req_open = [r for r in reqs if r[1] in ("WARN", "ACTION")]
    for name, verdict, detail in reqs:
        print("  %-9s %-42s %s" % (verdict, name, detail))
    print()

    compile_ok = True
    if not args.no_compile:
        print("-" * 86)
        print("CLEAN-ROOM COMPILE  (%d passes, outside the repository)" % args.passes)
        print("-" * 86)
        base = args.work or tempfile.mkdtemp(prefix="arxiv_v3_")
        if not os.path.isdir(base):
            os.makedirs(base)
        room = os.path.join(base, "bundle")
        mirror = os.path.join(base, "tree")
        for d in (room, mirror):
            if os.path.isdir(d):
                shutil.rmtree(d)
            os.makedirs(d)
        with tarfile.open(out, "r:gz") as tf:
            tf.extractall(room)
        # the mirror: the whole paper/ tree, same passes, same binary -- so that any
        # difference is the PACKAGE's fault and not the compiler's
        for r, dirs, fs in os.walk(PAPER):
            dirs[:] = [d for d in dirs if d not in ("_superseded", "__pycache__")]
            for f in fs:
                if os.path.splitext(f)[1].lower() in (".aux", ".log", ".out", ".fls",
                                                      ".toc", ".gz"):
                    continue
                s = os.path.join(r, f)
                dst = os.path.join(mirror, os.path.relpath(s, PAPER))
                if not os.path.isdir(os.path.dirname(dst)):
                    os.makedirs(os.path.dirname(dst))
                shutil.copy2(s, dst)
        print("  clean room : %s" % room)
        rc, log = run_latex(room, args.passes)
        st = log_stats(log)
        print("  pdflatex exit %s | pages %s | errors %d | undefined refs %d | "
              "undefined citations %d | missing files %d"
              % (rc, st["pages"], st["n_errors"], st["undef_ref"], st["undef_cit"],
                 st["missing_files"]))
        for e in st["errors"][:10]:
            print("      %s" % e.strip())
        pdf_room = os.path.join(room, "main.pdf")
        compile_ok = (st["n_errors"] == 0 and st["missing_files"] == 0
                      and os.path.isfile(pdf_room))

        # arXiv's Submission System 1.5 asks the SUBMITTER which processor to use.
        # Does the choice actually matter for this package, or is it boilerplate?
        # Measure it: run the other path (LaTeX in DVI mode) and see what happens.
        dvi = os.path.join(base, "dviroute")
        if os.path.isdir(dvi):
            shutil.rmtree(dvi)
        shutil.copytree(room, dvi)
        for junk in ("main.aux", "main.log", "main.out", "main.fls", "main.pdf"):
            jp = os.path.join(dvi, junk)
            if os.path.isfile(jp):
                os.remove(jp)
        pd = subprocess.Popen(["latex", "-interaction=nonstopmode", "-file-line-error",
                               "main.tex"], cwd=dvi, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT)
        pd.communicate()
        dlog = ""
        if os.path.isfile(os.path.join(dvi, "main.log")):
            dlog = open(os.path.join(dvi, "main.log"), "rb").read().decode("utf-8",
                                                                          "replace")
        dst = log_stats(dlog)
        print("  DVI-route probe (the OTHER processor arXiv offers at upload):")
        print("      latex main.tex -> %d error(s); first: %s"
              % (dst["n_errors"],
                 dst["errors"][0].strip()[:96] if dst["errors"] else "none"))
        print("      %s" % ("=> the upload form MUST be set to 'LaTeX with PDFLaTeX'; "
                            "one error per embedded PDF figure" if dst["n_errors"]
                            else "the package survives the DVI route as well"))
        shutil.rmtree(dvi, ignore_errors=True)

        fx = fls_cross_check(room, cl.files)
        if fx is not None:
            print("  .fls cross-check (TeX's own record of what it opened):")
            print("      files TeX opened inside the package : %d" % len(fx["opened"]))
            print("      opened but NOT shipped              : %d %s"
                  % (len(fx["tex_opened_not_shipped"]), fx["tex_opened_not_shipped"][:8]))
            print("      shipped but never opened            : %d %s"
                  % (len(fx["shipped_never_opened"]), fx["shipped_never_opened"][:8]))
            print("      read from the REPOSITORY (isolation leak) : %d %s"
                  % (len(fx["leaked"]), fx["leaked"][:6] or "- clean room is sealed"))
            print("      (%d INPUT lines total; the rest are texmf: the class, the "
                  "packages and the fonts)" % fx["n_input"])
            if (fx["tex_opened_not_shipped"] or fx["shipped_never_opened"]
                    or fx["leaked"]):
                compile_ok = False

        print("  recompiling the full tree with the same %d passes for comparison..."
              % args.passes)
        run_latex(mirror, args.passes, recorder=False)
        pdf_tree = os.path.join(mirror, "main.pdf")
        print()
        print("  PDF EQUIVALENCE")
        r1 = compare_pdfs(pdf_room, pdf_tree, "package, clean room",
                          "full paper/ tree, same passes")
        r2 = compare_pdfs(pdf_room, os.path.join(PAPER, "main.pdf"),
                          "package, clean room", "committed paper/main.pdf")
        boxes = pdf_media_boxes(pdf_room)
        print("  MediaBox of every page : %s  (%s)"
              % (boxes, "US Letter 612x792 pt" if boxes == [(612.0, 792.0)]
                 else "NON-STANDARD -- check the paper size"))
        compile_ok = compile_ok and r1["equal"] and (r2["pages_a"] == r2["pages_b"])
        if not r2["equal"]:
            # Which of the two is out of date?  Say so with a measurement, not a guess:
            # a committed PDF older than its own sources is stale, and then the
            # authority is the tree recompilation just above, not main.pdf.
            newer = stale_sources(cl)
            print("  committed paper/main.pdf differs on %d page(s) %s"
                  % (len(r2["diff_pages"]), r2["diff_pages"][:8]))
            if newer:
                print("  DIAGNOSIS: paper/main.pdf is STALE -- %d of its own sources "
                      "are newer than it:" % len(newer))
                for f, dt in newer[:8]:
                    print("      %-34s newer by %6.1f min" % (f, dt / 60.0))
                print("  The package matches what the CURRENT tree compiles to "
                      "(identical above); rebuild paper/main.pdf, then re-run this "
                      "script, before uploading.")
            else:
                print("  DIAGNOSIS: no source is newer than paper/main.pdf, so this "
                      "difference is the PACKAGE's fault. Do not upload.")
                compile_ok = False
        if not args.keep and args.work is None:
            shutil.rmtree(base, ignore_errors=True)
        else:
            print("  clean room kept at %s" % base)
        print()

    print("=" * 86)
    ok = (not bad) and (not fails) and (not req_fail) and compile_ok and deterministic
    print("RESULT: %s" % ("PASS -- the package is complete, clean and self-contained."
                          if ok else
                          "FAIL -- see the lines marked FAIL above."))
    if req_fail:
        print("        %d arXiv rule(s) BROKEN -- arXiv would reject or mangle this:"
              % len(req_fail))
        for name, _v, detail in req_fail:
            print("          FAIL   %s" % name)
    if req_open:
        print("        %d item(s) still open (a recommendation, or a choice made in "
              "the upload form rather than in the package):" % len(req_open))
        for name, v, detail in req_open:
            print("          %-6s %s" % (v, name))
    print("=" * 86)
    return 0 if ok else 1


# ==================================================================================
# 8.  negative controls
# ==================================================================================
def selftest(args):
    """Five mutations, five checks that must FAIL.  A PASS proves nothing on its own."""
    print("=" * 86)
    print("NEGATIVE CONTROLS  (each mutation must make its check FAIL)")
    print("=" * 86)
    tmp = tempfile.mkdtemp(prefix="arxiv_v3_negctl_")
    results = []
    try:
        cl, _sd = build_closure()

        # --- control 1: the bracket-requiring parser, the exact old bug -------------
        saved = _consume_args.__doc__
        n_full = len(cl.files)
        strict = _closure_with_required_brackets()
        results.append(("1. parser that REQUIRES table[...] loses .dat files",
                        len(strict) < n_full,
                        "%d files instead of %d -- %d lost: %s"
                        % (len(strict), n_full, n_full - len(strict),
                           sorted(set(cl.files) - set(strict))[:12])))
        del saved

        # --- control 2: a missing .dat must kill the clean-room compile -------------
        out = os.path.join(tmp, "b.tar.gz")
        make_tarball(cl, out)
        room = os.path.join(tmp, "room")
        os.makedirs(room)
        with tarfile.open(out, "r:gz") as tf:
            tf.extractall(room)
        victim = next(f for f in cl.files if f.endswith(".dat"))
        os.remove(os.path.join(room, victim.replace("/", os.sep)))
        rc, log = run_latex(room, 1, recorder=False)
        st = log_stats(log)
        # This is also the POSITIVE CONTROL on log_stats itself: the assertion is on
        # n_errors, not on missing_files, because n_errors is the number that was
        # silently blind to the -file-line-error shape until 2026-09-19.
        results.append(("2. deleting %s is SEEN BY log_stats as an error" % victim,
                        st["n_errors"] > 0,
                        "errors %d (%s), missing-file warnings %d"
                        % (st["n_errors"],
                           (st["errors"][0].strip()[:60] if st["errors"] else "-"),
                           st["missing_files"])))

        # --- control 2b: a missing .tex must be seen too ----------------------------
        room1b = os.path.join(tmp, "room1b")
        os.makedirs(room1b)
        with tarfile.open(out, "r:gz") as tf:
            tf.extractall(room1b)
        os.remove(os.path.join(room1b, "sec_10.tex"))
        _rc, log1b = run_latex(room1b, 1, recorder=False)
        st1b = log_stats(log1b)
        results.append(("2b. deleting sec_10.tex is seen as an error",
                        st1b["n_errors"] > 0,
                        "errors %d (%s)" % (st1b["n_errors"],
                                            st1b["errors"][0].strip()[:60]
                                            if st1b["errors"] else "-")))

        # --- control 3: a forbidden token must be caught by the audit ---------------
        bad_dir = os.path.join(tmp, "bad")
        shutil.copytree(PAPER, bad_dir, ignore=shutil.ignore_patterns(
            "_superseded", "__pycache__", "*.aux", "*.log", "*.out", "*.fls", "*.gz"))
        tgt = os.path.join(bad_dir, "figs", "fig_scaling2_native.tex")
        txt = open(tgt, "rb").read().decode("utf-8", "replace")
        open(tgt, "wb").write(txt.replace("19183", "19184", 1).encode("utf-8"))
        out2 = os.path.join(tmp, "bad.tar.gz")
        _fake_tarball(cl, bad_dir, out2)
        f2, _ = audit_bundle(out2, paper=bad_dir)
        results.append(("3. re-introducing 19184 is caught by the audit",
                        any("19184" in x for x in f2),
                        "; ".join(x for x in f2 if "19184" in x)[:120] or "NOT CAUGHT"))

        # --- control 4: an edited .tex must change the PDF text flow ----------------
        room2 = os.path.join(tmp, "room2")
        os.makedirs(room2)
        with tarfile.open(out, "r:gz") as tf:
            tf.extractall(room2)
        s10 = os.path.join(room2, "sec_10.tex")
        t10 = open(s10, "rb").read().decode("utf-8", "replace")
        open(s10, "wb").write((t10 + "\n\nNEGATIVE CONTROL SENTINEL PARAGRAPH.\n")
                              .encode("utf-8"))
        run_latex(room2, 2, recorder=False)
        pages = pdf_text_pages(os.path.join(room2, "main.pdf"))
        hit = any("NEGATIVE CONTROL SENTINEL" in p for p in pages)
        results.append(("4. the PDF text comparison can see an inserted paragraph",
                        hit, "sentinel found in the extracted text: %s" % hit))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    bad = 0
    for name, fired, detail in results:
        print("  %-6s %-58s %s" % ("FIRED" if fired else "SILENT", name, detail))
        if not fired:
            bad += 1
    print("=" * 86)
    print("RESULT: %s" % (("PASS -- all %d controls fired; the checks have teeth."
                           % len(results)) if bad == 0 else
                          ("FAIL -- %d of %d control(s) stayed silent; those checks "
                           "are decorative." % (bad, len(results)))))
    print("=" * 86)
    return 0 if bad == 0 else 1


def _closure_with_required_brackets():
    """Reproduce the OLD, broken derivation: `table` only counts when it carries an
    option bracket.  Used solely as a negative control."""
    rx_in = re.compile(r"\\(?:input|include)\s*\{([^}]*)\}")
    rx_gr = re.compile(r"\\includegraphics\s*(?:\[[^\]]*\])?\s*\{([^}]*)\}")
    rx_tb = re.compile(r"table\s*\[[^\]]*\]\s*\{([^}]*)\}")
    files, seen, queue = [MAIN], set([MAIN]), [MAIN]
    while queue:
        rel = queue.pop(0)
        if not rel.endswith(".tex"):
            continue
        body = strip_tex_comments(open(os.path.join(PAPER, rel.replace("/", os.sep)),
                                       "rb").read().decode("utf-8", "replace"))
        for rx, kind in ((rx_in, "input"), (rx_gr, "graphic"), (rx_tb, "table")):
            for t in rx.findall(body):
                got = _resolve(t.strip(), kind, ["", "figs"], rel)
                if got and got not in seen:
                    seen.add(got)
                    files.append(got)
                    queue.append(got)
    return files


def _fake_tarball(cl, paper_dir, out_path):
    members = sorted(cl.files)
    buf = io.BytesIO()
    tf = tarfile.open(fileobj=buf, mode="w", format=tarfile.USTAR_FORMAT)
    for rel in members:
        src = os.path.join(paper_dir, rel.replace("/", os.sep))
        if not os.path.isfile(src):
            continue
        data = open(src, "rb").read()
        ti = tarfile.TarInfo(rel)
        ti.size, ti.mtime, ti.mode = len(data), 0, 0o644
        ti.uid = ti.gid = 0
        ti.uname = ti.gname = ""
        tf.addfile(ti, io.BytesIO(data))
    tf.close()
    gz = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=gz, mtime=0) as g:
        g.write(buf.getvalue())
    open(out_path, "wb").write(gz.getvalue())


if __name__ == "__main__":
    sys.exit(main())
