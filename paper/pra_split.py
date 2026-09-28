# -*- coding: utf-8 -*-
r"""Build and check the two files Physical Review A asks for: main_pra.pdf and sm_pra.pdf.

PRA takes the Supplemental Material as a separate file.  main.tex compiles the main text,
the references and the SM as ONE document (the arXiv form, which stays as it is).  Two
drivers next to it split that document without copying any of it:

    main_pra.tex   runs main.tex itself and ends the document where main.tex opens the SM
                   (\input{sm_begin.tex}), i.e. right after \input{bibliography.tex}
    sm_pra.tex     runs main.tex itself, skips everything before \input{sm_begin.tex}, runs
                   the SM \input list of main.tex, and ends with a reference list of the works
                   the SM cites, in the order it first cites them, whose entries are read
                   from bibliography.tex at build time

Neither driver holds a sentence, a number or a reference entry of the paper, so the
guardians -- which read every paper/*.tex -- see nothing twice, and neither driver can
drift from main.tex.  The comments at the top of each driver say what they add.

This script builds both in a temporary copy of paper/ (it writes nothing in paper/ except
main_pra.pdf and sm_pra.pdf) and checks, against a build of main.tex made in the same copy
and against the committed paper/main.pdf:

  * the logs: 0 errors, 0 undefined references or citations, 0 multiply defined labels or
    citations, no rerun request, no overfull box that the main.tex build does not also
    have, no warning line that the main.tex log does not also have;
  * the SM reference list: its entries are exactly the keys the SM cites, in the order of
    their first citation both as typeset (sm_pra.aux) and in the source (main.tex's SM
    \input list, expanded);
  * the APS convention for the SM (added 2026-09-28): the main text cites the SM as a
    reference (\bibitem{supplemental} in bibliography.tex, "See Supplemental Material at
    [URL will be inserted by publisher] for ..."), the main reference list is in first-
    citation order, and the entry's "which includes Refs. [a-b]" names exactly the block of
    works cited only in the SM (check_sm_citation);
  * the pages: every page of main.pdf is in exactly one of the two PDFs.  Each page is
    rendered and compared pixel by pixel below the running head (REVTeX prints only the
    page number there); a page that is not pixel-identical is compared as text -- the
    multiset of its characters without the citation numbers, its number of citation
    groups -- and by the height of every line on it.  '??' anywhere is a failure.

Usage, from the repository root or from paper/:
    python paper/pra_split.py --build      # build, check, copy the two PDFs into paper/
    python paper/pra_split.py --compare    # only the page comparison, on the PDFs in paper/

Needs pdflatex (no BibTeX: the bibliography is hand-written), pdftotext, pdftoppm and Pillow.
xr-hyper's [nocite] option, used by both drivers, needs an xr-hyper of 2023 or later
(built here with v7.01o, 2025-07-12).
"""
from __future__ import print_function

import collections
import hashlib
import io
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
COMMENT = re.compile(r"(?<!\\)%.*")
CITE = re.compile(r"\\(?:cite|citep|citet|onlinecite|citealp|citenum|nocite)\*?"
                  r"(?:\[[^\]]*\])*\{([^}]*)\}")


def read(d, name):
    with io.open(os.path.join(d, name), encoding="utf-8") as fh:
        return fh.read()


def expand(d, name):
    """Inline \\input expansion with comments removed, the way the guardians do it."""
    path = name if os.path.isfile(os.path.join(d, name)) else name + ".tex"
    s = COMMENT.sub("", read(d, path))
    out, pos = [], 0
    for m in re.finditer(r"\\input\{([^}]*)\}", s):
        out.append(s[pos:m.start()])
        out.append(expand(d, m.group(1)))
        pos = m.end()
    out.append(s[pos:])
    return "".join(out)


def first_cites(text):
    order = []
    for m in CITE.finditer(text):
        for k in m.group(1).split(","):
            k = k.strip()
            if k and k not in order:
                order.append(k)
    return order


def check_sm_refs(d):
    """The SM reference list against the SM's citations, typeset and in the source."""
    rep, ok = [], True
    aux = io.open(os.path.join(d, "sm_pra.aux"), encoding="latin-1").read()
    cited = []
    for g in re.findall(r"^\\citation\{([^}]*)\}", aux, re.M):
        for k in g.split(","):
            k = k.strip()
            if k and k not in cited and not k.endswith("Control"):
                cited.append(k)
    listed = re.findall(r"^\\bibcite\{([^}]*)\}\{\{?(\d+)", aux, re.M)
    keys = [k for k, _ in listed]
    nums = [int(n) for _, n in listed]
    body = expand(d, "main.tex")
    src = first_cites(body[body.index("\\end{thebibliography}"):])
    defined = re.findall(r"^\\bibitem\{([^}]*)\}", COMMENT.sub("", read(d, "bibliography.tex")),
                         re.M)
    rep.append("SM reference list: %d entries, numbered %s; the SM cites %d works"
               % (len(keys), "1-%d in order" % len(nums) if nums == list(range(1, len(nums) + 1))
                  else "OUT OF ORDER", len(cited)))
    for what, got in (("the order of first citation as typeset (sm_pra.aux)", cited),
                      ("the order of first citation in the source (main.tex, SM part)", src)):
        same = keys == got
        ok &= same
        rep.append("  list = %s: %s" % (what, "yes" if same else "NO"))
    missing = [k for k in keys if k not in defined]
    ok &= not missing and nums == list(range(1, len(nums) + 1))
    rep.append("  every entry is a \\bibitem of bibliography.tex: %s"
               % ("yes" if not missing else "NO: %s" % missing))
    return ok, rep


def check_sm_citation(d):
    r"""The APS convention for the Supplemental Material, checked in the source.

    The main text must cite the SM as a reference (\bibitem{supplemental}, "See Supplemental
    Material at [URL will be inserted by publisher] for ..."), and the works cited only in the
    SM must be in the main reference list, the SM entry naming them ("which includes Refs.
    [a-b]", written with \citenum{first}--\citenum{last}).  Checked here: the entry exists and
    is cited in the main text; the reference list is in first-citation order (main text, then
    the works cited only in the SM); and the two keys of the range are the first and the last
    work cited only in the SM, so the range is exactly that block."""
    rep, ok = [], True
    body = expand(d, "main.tex")
    b0, b1 = body.index("\\begin{thebibliography}"), body.index("\\end{thebibliography}")
    main_order = first_cites(body[:b0])
    sm_order = first_cites(body[b1:])
    sm_only = [k for k in sm_order if k not in main_order]
    bib_src = COMMENT.sub("", read(d, "bibliography.tex"))
    keys = re.findall(r"^\\bibitem\{([^}]*)\}", bib_src, re.M)
    m = re.search(r"^\\bibitem\{supplemental\}(.*)$", bib_src, re.M)
    rng = re.search(r"\\citenum\{([^}]*)\}--\\citenum\{([^}]*)\}", m.group(1)) if m else None
    cited = "supplemental" in main_order
    ok &= bool(m) and cited
    rep.append("SM as a reference: \\bibitem{supplemental} %s, cited in the main text: %s"
               % ("present" if m else "MISSING", "yes" if cited else "NO"))
    order_ok = keys == main_order + sm_only
    ok &= order_ok
    rep.append("  reference list in first-citation order (main text, then SM only): %s"
               % ("yes" if order_ok else "NO"))
    if rng and sm_only:
        a, b = rng.group(1), rng.group(2)
        same = (a, b) == (sm_only[0], sm_only[-1]) and keys[keys.index(a):keys.index(b) + 1] \
            == sm_only if a in keys and b in keys else False
        ok &= same
        rep.append("  'which includes Refs.' names %s--%s = entries %d-%d; the %d works cited only "
                   "in the SM are %s--%s: %s"
                   % (a, b, keys.index(a) + 1 if a in keys else -1,
                      keys.index(b) + 1 if b in keys else -1, len(sm_only), sm_only[0],
                      sm_only[-1], "same block" if same else "DIFFERENT"))
    else:
        ok = False
        rep.append("  'which includes Refs.' range: NOT FOUND in the entry")
    return ok, rep


def run_latex(d, job):
    p = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                        job + ".tex"], cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return p.returncode


def digest(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest() if os.path.exists(path) else ""


def log_report(d, job):
    log = io.open(os.path.join(d, job + ".log"), encoding="latin-1").read()
    errors = re.findall(r"^! .*", log, re.M)
    undef = re.findall(r"^(?:LaTeX|Package \w+) Warning: .*(?:undefined|Undefined).*", log, re.M)
    undef += re.findall(r"^LaTeX Warning: There were undefined .*", log, re.M)
    multi = re.findall(r"^.*Warning: .*multiply[ -]defined.*", log, re.M | re.I)
    undef += re.findall(r"^Package xr Warning: .*", log, re.M)            # "No file ..."
    undef += re.findall(r"^No file .*\.aux\.", log, re.M)   # (No file <job>.bbl is REVTeX's
    #                                                      own probe, and main.tex's too)
    rerun = re.search(r"Rerun to get|may have changed", log) is not None
    over = re.findall(r"^Overfull \\[hv]box \(([\d.]+)pt too \w+\)[^\n]*", log, re.M)
    overl = re.findall(r"^(Overfull \\[hv]box \([\d.]+pt too \w+\)[^\n]*)", log, re.M)
    warn = re.findall(r"^.*Warning.*$", log, re.M)
    return {"errors": errors, "undefined": undef, "multiply": multi, "rerun": rerun,
            "overfull": overl, "warnings": warn}


def pages_text(pdf):
    """One string per page.  Reading order (no -layout): -layout interleaves the glyphs of
    displayed mathematics differently from one PDF to the next even when the pages are
    pixel-identical, so it cannot be compared as a sequence."""
    out = subprocess.run(["pdftotext", "-enc", "UTF-8", pdf, "-"],
                         stdout=subprocess.PIPE).stdout.decode("utf-8", "replace")
    pages = out.split("\f")
    if pages and not pages[-1].strip():
        pages = pages[:-1]
    return pages


CITEGRP = re.compile(r"\[(?:\d+(?:\s*[,–-]\s*)?)+\]")


def page_chars(page):
    """The multiset of non-blank characters of a page, without its running head (the page
    number, first line) and without the numeric citation groups [n], [n, m], [n-m]."""
    lines = [l.strip() for l in page.split("\n") if l.strip()]
    if lines and re.fullmatch(r"\d+", lines[0]):
        lines = lines[1:]
    t = "\n".join(lines)
    ncite = len(CITEGRP.findall(t))
    return collections.Counter(c for c in CITEGRP.sub("", t) if not c.isspace()), ncite


def render(pdf, first, last, dpi, outdir, tag):
    from PIL import Image
    subprocess.run(["pdftoppm", "-r", str(dpi), "-gray", "-f", str(first), "-l", str(last),
                    pdf, os.path.join(outdir, tag)], check=True)
    names = sorted(n for n in os.listdir(outdir) if n.startswith(tag + "-"))
    return [Image.open(os.path.join(outdir, n)).convert("L") for n in names]


def line_bands(img, top):
    """The vertical line layout of a page below the running head: the runs of rows that
    carry ink, as (first row, last row)."""
    from PIL import ImageOps
    w, h = img.size
    inv = ImageOps.invert(img.crop((0, top, w, h))).point(lambda v: 255 if v > 64 else 0)
    px = inv.load()
    rows = [y + top for y in range(h - top) if any(px[x, y] for x in range(w))]
    bands = []
    for y in rows:
        if bands and y == bands[-1][1] + 1:
            bands[-1][1] = y
        else:
            bands.append([y, y])
    return bands


def same_lines(a, b, tol=4):
    """Same number of lines at the same heights, to within `tol` pixels (a bracket or a
    descender that changes line can move the edge of a band by a pixel or two)."""
    return len(a) == len(b) and all(abs(x[0] - y[0]) <= tol and abs(x[1] - y[1]) <= tol
                                    for x, y in zip(a, b))


def pixel_diff(a, b, top):
    from PIL import ImageChops
    w, h = a.size
    if b.size != a.size:
        return (0, 0, w, h)
    d = ImageChops.difference(a.crop((0, top, w, h)), b.crop((0, top, w, h)))
    return d.point(lambda v: 255 if v > 64 else 0).getbbox()


def compare(main_pdf, mp_pdf, sm_pdf, dpi=100):
    """Page-by-page comparison.  Every page is rendered and compared pixel by pixel below the
    running head (the page number is the only thing REVTeX prints there); a page that is not
    pixel-identical is then compared as text: the multiset of its characters without the
    citation numbers, the number of citation groups, and the rows that carry ink."""
    M, A, B = pages_text(main_pdf), pages_text(mp_pdf), pages_text(sm_pdf)
    head = int(0.6 * dpi)                      # the running head lies above 0.6 in
    rep = ["main.pdf %d pp, main_pra.pdf %d pp, sm_pra.pdf %d pp" % (len(M), len(A), len(B))]
    work = tempfile.mkdtemp(prefix="pra_cmp_")
    nA = len(A)
    nS = len(M) - nA
    imM = render(main_pdf, 1, len(M), dpi, work, "m")
    imA = render(mp_pdf, 1, len(A), dpi, work, "a")
    imB = render(sm_pdf, 1, len(B), dpi, work, "b")
    ok = True
    counts = collections.Counter()
    for label, pairs in (("main_pra", [(i, i) for i in range(nA)]),
                         ("sm_pra", [(nA + i, i) for i in range(nS)])):
        other_img, other_txt = (imA, A) if label == "main_pra" else (imB, B)
        for mi, oi in pairs:
            if oi >= len(other_img):
                ok = False
                rep.append("  %s has no p.%d for main.pdf p.%d" % (label, oi + 1, mi + 1))
                continue
            bb = pixel_diff(imM[mi], other_img[oi], head)
            if bb is None:
                counts[label, "identical"] += 1
                continue
            cm, nm = page_chars(M[mi])
            co, no = page_chars(other_txt[oi])
            miss, extra = cm - co, co - cm
            lines_ok = same_lines(line_bands(imM[mi], head), line_bands(other_img[oi], head))
            # Residue that is only digits, brackets, commas and dashes: a citation group that
            # pdftotext interleaved with a neighbouring formula (so the regex could not lift
            # it out), or a hyphen that became a line-end hyphen when a narrower citation
            # number moved a line break inside the paragraph.  Reported, not hidden.
            citeish = set("0123456789[],–-")
            residue_ok = all(c in citeish for c in list(miss) + list(extra))
            if label == "sm_pra" and residue_ok and lines_ok:
                kind = ("citation numbers only" if not (miss or extra) else
                        "citation numbers only; text extraction differs by '%s' / '%s'"
                        % ("".join(sorted(miss.elements())), "".join(sorted(extra.elements()))))
                counts[label, "citation numbers only"] += 1
                rep.append("  sm_pra p.%d vs main.pdf p.%d: %s (%d vs %d citation groups found; "
                           "same lines, same heights)" % (oi + 1, mi + 1, kind, nm, no))
                continue
            ok = False
            rep.append("  %s p.%d vs main.pdf p.%d: DIFFERENT (pixels in %s; %d characters "
                       "only in main.pdf, %d only in %s; citation groups %d vs %d)"
                       % (label, oi + 1, mi + 1, bb, sum(miss.values()), sum(extra.values()),
                          label, nm, no))
            rep.append("      only in main.pdf: %s" % "".join(sorted(miss.elements()))[:200])
            rep.append("      only in %s: %s" % (label, "".join(sorted(extra.elements()))[:200]))
    rep.append("  main_pra pp. 1-%d vs main.pdf pp. 1-%d: %d pixel-identical below the page "
               "number" % (nA, nA, counts["main_pra", "identical"]))
    rep.append("  sm_pra pp. 1-%d vs main.pdf pp. %d-%d: %d pixel-identical below the page "
               "number, %d differ in citation numbers only (same lines at the same heights; "
               "any residue of the text comparison is listed above)"
               % (nS, nA + 1, len(M), counts["sm_pra", "identical"],
                  counts["sm_pra", "citation numbers only"]))
    extra = B[nS:]
    if extra:
        rep.append("  sm_pra pp. %d-%d: the SM reference list (%d pages; not in main.pdf by "
                   "construction)" % (nS + 1, len(B), len(extra)))
    qq = [(n, i + 1) for n, P in (("main_pra", A), ("sm_pra", B)) for i, p in enumerate(P)
          if "??" in p]
    rep.append("  '??' anywhere in the text of the two PDFs: %s" % (qq or "none"))
    shutil.rmtree(work, ignore_errors=True)
    return ok and not qq, rep


def build():
    tmp = tempfile.mkdtemp(prefix="pra_split_")
    d = os.path.join(tmp, "paper")
    shutil.copytree(HERE, d, ignore=shutil.ignore_patterns(
        "*.aux", "*.log", "*.out", "*.bbl", "*.pdf.tmp", "arxiv-submission.tar.gz",
        "main.pdf", "main_pra.pdf", "sm_pra.pdf"))
    print("building in", d)
    # the reference build: main.tex itself, for the overfull comparison and the page check
    for _ in range(3):
        run_latex(d, "main")
    base = log_report(d, "main")
    # the two drivers: alternate until neither asks for a rerun and neither .aux moves
    prev = None
    for rnd in range(1, 8):
        run_latex(d, "main_pra")
        run_latex(d, "sm_pra")
        cur = (digest(os.path.join(d, "main_pra.aux")), digest(os.path.join(d, "sm_pra.aux")))
        r1, r2 = log_report(d, "main_pra"), log_report(d, "sm_pra")
        if cur == prev and not r1["rerun"] and not r2["rerun"]:
            break
        prev = cur
    print("pdflatex rounds (main_pra + sm_pra each):", rnd)
    ok, rep_refs = check_sm_refs(d)
    ok_cit, rep_cit = check_sm_citation(d)
    ok &= ok_cit
    rep_refs += rep_cit
    base_over = sorted(base["overfull"])
    for job, r in (("main", base), ("main_pra", r1), ("sm_pra", r2)):
        print("%-9s errors %d, undefined %d, multiply-defined %d, rerun requested %s, overfull %d"
              % (job, len(r["errors"]), len(r["undefined"]), len(r["multiply"]), r["rerun"],
                 len(r["overfull"])))
        for x in r["errors"] + r["undefined"] + r["multiply"]:
            print("    " + x)
        for x in r["overfull"]:
            print("    " + x + ("" if x in base_over else "   <-- NOT in the main.tex build"))
        if job != "main":
            new = sorted(set(r["warnings"]) - set(base["warnings"]))
            print("    warning lines not in the main.tex log: %d" % len(new))
            for x in new:
                print("      " + x)
            ok &= not new
            ok &= not (r["errors"] or r["undefined"] or r["multiply"] or r["rerun"])
            ok &= all(x in base_over for x in r["overfull"])
    print("\n".join(rep_refs))
    for job in ("main_pra", "sm_pra"):
        if not os.path.exists(os.path.join(d, job + ".pdf")):
            print("NO PDF for %s: nothing copied; see %s" % (job, d))
            return False
    for job in ("main_pra", "sm_pra"):
        shutil.copyfile(os.path.join(d, job + ".pdf"), os.path.join(HERE, job + ".pdf"))
    print("copied main_pra.pdf and sm_pra.pdf into", HERE)
    good, rep = compare(os.path.join(HERE, "main.pdf"), os.path.join(HERE, "main_pra.pdf"),
                        os.path.join(HERE, "sm_pra.pdf"))
    fresh = pages_text(os.path.join(d, "main.pdf")) == pages_text(os.path.join(HERE, "main.pdf"))
    print("\nPAGE COMPARISON against paper/main.pdf (the committed PDF; a fresh build of "
          "main.tex in the same copy has %s text)" % ("the SAME" if fresh else "DIFFERENT"))
    print("\n".join(rep))
    print("\nbuild tree kept for inspection:", d)
    return ok and good and fresh


if __name__ == "__main__":
    if "--compare" in sys.argv[1:]:          # compare the PDFs already in paper/, build nothing
        good, rep = compare(os.path.join(HERE, "main.pdf"), os.path.join(HERE, "main_pra.pdf"),
                            os.path.join(HERE, "sm_pra.pdf"))
        print("\n".join(rep))
        sys.exit(0 if good else 1)
    if "--build" in sys.argv[1:]:
        sys.exit(0 if build() else 1)
    print(__doc__)
