# -*- coding: utf-8 -*-
r"""Check the manuscript against ITSELF: does every section still say the same thing?

WHY THIS FILE EXISTS
--------------------
On 2026-09-19 an adversarial pass found seven contradictions BETWEEN sections while all
three existing guardians were green, and named the reason exactly:

    "verify.py compares a PRINTED NUMBER against a DEPOSITED FILE, file by file; it does
     not compare a CLAIM against a CLAIM.  A verify PASS is therefore not evidence of
     coherence between sections."

That is true, and it is the blind spot that let Sec. V be rewritten under the abstract,
the introduction, the conclusion, the resource table and the data-availability statement
without any of the five being updated.  The contradictions were: a vacuity range of
"2.8 to 4.9" printed in four places including the abstract, while the new Sec. V table
gives 1.66-8.61; "at the four sizes evaluated" against five rows; "Two sizes are what we
measured" against a support measured at L = 8, 10, 12 and 14; a hardware shot budget that
read 5e4 in three places and 3.5e5 in a fourth; "uniformly" in the conclusion against
thirty lines of Sec. V refuting it; "sampled fraction" against Sec. VII's own written
rule; and a data-availability statement quoting guardian counts that had moved.

This file is the guardian for that class.  Every check has TWO halves -- the correct
statement must be PRESENT at least as often as it was on the clean tree, and the
superseded one must be ABSENT -- because a check that only looks for the right string
passes on a manuscript that says both, which is exactly what happened.

WHAT IT DOES NOT DO
-------------------
It does not check physics (`verify.py`), it does not re-run generators
(`check_figures.py`), and it does not check the documentation (`check_provenance.py`).
It checks only that the manuscript does not contradict itself or the deposited grid.

SELF-TEST FIRST
---------------
A PASS is worth nothing without a control that fires.  Every check gets two: one
occurrence of the correct statement is deleted, and the superseded statement is pasted
back in.  Both must make the check fail.  If any control does not fire, this file exits 2
and reports NO result, rather than a green one.
"""
from __future__ import print_function

import io
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.join(ROOT, "paper")

RE_INPUT = re.compile(r"\\(?:input|include)\s*\{([^}]*)\}")


# ---------------------------------------------------------------------------
#  The body, as the reader sees it: \input expanded, comments removed.
#  Scanning the .tex files raw gives both false alarms (a commented-out sentence
#  prints nothing) and false zeros (a sentence assembled through \input is missed).
# ---------------------------------------------------------------------------
def strip_comment(line):
    r"""Everything from an unescaped % to end of line is not source; `\%` is a percent."""
    out, i = [], 0
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


def read_rel(rel):
    for cand in (rel, rel + ".tex"):
        p = os.path.join(PAPER, cand.replace("/", os.sep))
        if os.path.isfile(p):
            return io.open(p, encoding="utf-8", newline="").read()
    return None


def body(start="main.tex", stack=None):
    stack = stack or []
    if start in stack:
        return ""
    raw = read_rel(start)
    if raw is None:
        return ""
    raw = raw.replace("\r\n", "\n").replace("\r", "\n")
    out = []
    for line in raw.split("\n"):
        code = strip_comment(line)
        pos = 0
        while True:
            m = RE_INPUT.search(code, pos)
            if not m:
                out.append(code[pos:])
                break
            out.append(code[pos:m.start()])
            out.append(body(m.group(1), stack + [start]))
            pos = m.end()
        out.append("\n")
    return "".join(out)


def flat(s):
    return re.sub(r"\s+", " ", s)


# ---------------------------------------------------------------------------
#  Ground truth that is NOT typed here: it is re-derived from the deposit.
# ---------------------------------------------------------------------------
def vacuity_range():
    """The vacuity factor of each of the twenty published-fraction cells.

    Sixteen (L = 6, 8, 10, 12) live in PUB_rows.json; the four L = 14 cells live in
    L14_FR008.json and must be taken at the DEEPEST recursion, which is the depth the
    manuscript quotes.  Reading only the first file gives 8.04 instead of 8.61 and a
    green run for the wrong reason, so both are read and the depth is taken as the
    maximum present rather than assumed.
    """
    d = os.path.join(ROOT, "data", "c3_frontier", "published_fraction")
    vals = {}
    for r in json.load(io.open(os.path.join(d, "PUB_rows.json"), encoding="utf-8")):
        vals[(r["L"], r["eta"])] = r["vac"]
    blk = json.load(io.open(os.path.join(d, "L14_FR008.json"), encoding="utf-8"))[0]
    deepest = max(x["n"] for x in blk["rows"])
    for x in blk["rows"]:
        if x["n"] == deepest:
            vals[(blk["L"], x["eta"])] = x["bound_leak"] / x["bound_trivial"]
    return vals


# ---------------------------------------------------------------------------
#  The checks.  Each is (name, (regex, floor), absent, poison, why).
#
#  floor  -- the count MEASURED on the clean tree of 2026-09-19, not a guess.  That is
#            what makes DELETING one of the places a failure instead of a quietly
#            smaller number.  Raise it when a place is added; never lower it to make a
#            run go green.
#  absent -- the superseded form.  It must appear zero times.
#  poison -- a literal string that must make this check fail; it is its own control.
# ---------------------------------------------------------------------------
def build_checks(txt):
    vals = vacuity_range()
    lo_s = "%.2f" % min(vals.values())
    hi_s = "%.2f" % max(vals.values())
    nsizes = len(set(L for L, _ in vals))
    word = {4: "four", 5: "five", 6: "six"}[nsizes]
    rng = re.escape("$%s$" % lo_s) + r".{0,40}" + re.escape("$%s$" % hi_s)

    return [
        ("the vacuity range is the deposited one, everywhere it is stated",
         (rng, 5),
         [r"\$2\.8\$\s*(?:--|to)\s*\$4\.9\$", r"2\.8 to 4\.9"],
         "vacuous by factors of $2.8$ to $4.9$",
         "the abstract, Sec. I, Sec. V, Sec. IX and Sec. X each state the range over "
         "which the certificate is vacuous. It is %s-%s over the twenty cells of the "
         "deposited published-fraction grid, and the superseded 2.8-4.9 must appear "
         "nowhere." % (lo_s, hi_s)),

        ("the number of published sizes agrees with the deposited grid",
         (r"at all %s sizes|at each of the %s sizes evaluated|at the %s sizes measured"
          r"|%s sizes up to a sector|at %s system sizes|all five sizes and all four"
          % (word, word, word, word, word), 8),
         [r"at all four sizes", r"at each of the four sizes evaluated",
          r"at the four sizes measured", r"four sizes up to a sector",
          r"at four system sizes"],
         "at all four sizes",
         "the deposited grid carries %d sizes, so no section may still count four"
         % nsizes),

        ("the Krylov-support claim is the measured one, not an open question",
         (r"L=8\$, \$10\$, \$12\$ and\s*\$14", 4),
         [r"measured at two sizes and we state it for",
          r"Two sizes are what we",
          r"whether it persists at \$L=10\$ and \$L=12\$ is open"],
         "Two sizes are what we measured",
         "Sec. V measures the coordinate support of K(H,phi) at L = 8, 10, 12 and 14. "
         "The introduction and the conclusion must not still offer it as open, and the "
         "cost of settling it must not still be quoted as unspent."),

        ("the hardware shot budget is one number in all four places",
         (r"3\.5\\times10\^\{5\}", 3),
         [r"pooled over the \$t\{=\}0\$ reference and \$K\{=\}7\$ evolution circuits",
          r"\$5\\times10\^\{4\}\$ computational-basis shots pooled over",
          r"reconstructed from \$50\{,\}000\$ computational-basis bitstrings"],
         "$5\\times10^{4}$ computational-basis shots pooled over the reference",
         "the ibm_fez run is seven circuits at 5e4 shots each = 3.5e5. The provenance "
         "table, Sec. IX A, Sec. IX C and the Fig. 14 caption must agree, and the "
         "factor 14.4 that Sec. IX C derives is computed from 3.5e5."),

        ("the Fig. 14 caption carries the same total",
         (r"\$350\{,\}000\$ computational-basis bitstrings", 1),
         [r"from \$50\{,\}000\$ computational-basis bitstrings"],
         "reconstructed from $50{,}000$ computational-basis bitstrings",
         "the hardware hero caption is where a referee checks the shot count first"),

        ("the crossing is a drift, not 'uniformly'",
         (r"a drift(?:,| and) not scatter", 2),
         [r"uniformly over \$L=6\$--\$12\$"],
         "uniformly over $L=6$--$12$",
         "Sec. V spends thirty lines showing that the factor of three is an ordered "
         "drift and not scatter. 'Uniformly' is the opposite claim, and it stood in "
         "the conclusion."),

        ("the sector fraction is named as Sec. VII's own rule requires",
         (r"published fraction", 17),
         [r"sampled fraction"],
         "at every sampled fraction",
         "Sec. VII states in the printed text that describing the T->infinity curve as "
         "SAMPLED without a shot count is to mislabel it, and that it is labelled the "
         "published curve from there on. The abstract and Sec. X broke that rule."),

        ("the declared gap of Fig. 3(a,b) is closed, and said to be closed",
         (r"has since been run at finite shots", 1),
         [r"Declared gap"],
         "Declared gap, see text",
         "a finite-shot run of the literal configuration of Fig. 3(a,b) now exists "
         "(5.2e6 shots per channel), so the provenance table and Sec. II B must not "
         "still declare it missing."),

        ("the data-availability statement quotes the guardian it actually has",
         # 9843 -> 9844 on 2026-09-19: adding src/check_trim.py moved the README
         # script-count check from 58 to 59 and the assertion total with it.  An
         # assertion count is derived from the tree, so it goes stale whenever a
         # file is added -- which is why this check pins it.  It fired correctly.
         (r"\$831\$ checks over \$9\\,844\$ numeric assertions", 1),
         [r"\$807\$ checks", r"\$7\\,454\$ numeric assertions",
          r"four known open defects"],
         "$807$ checks over $7\\,454$ numeric assertions",
         "Sec. IX F tells the reader what `python src/verify.py` prints. Run it and "
         "compare. If verify.py's counts have moved, this check is the reminder that "
         "Sec. IX F, README.md, docs/DATA_AVAILABILITY.md, docs/REPRODUCE.md and "
         "docs/KNOWN_DISCREPANCIES.md all move with it."),

        ("the largest sector quoted is the largest sector evaluated",
         (r"10\\,306\\,296", 10),
         [r"to dimension \$731\\,808\$ ours",
          r"sector dimensions \$300\$ to \$731\\,808\$",
          r"four sizes up to a sector of dimension \$731\\,808\$"],
         "sector dimensions $300$ to $731\\,808$",
         "the abstract, Sec. I, Sec. IX and Sec. X state the largest sector the "
         "certificate was evaluated on. The L = 14 row made it 10,306,296."),
    ], lo_s, hi_s, nsizes


def run(txt, checks, verbose=True):
    bad = []
    for name, (pat, floor), absent, poison, why in checks:
        n = len(re.findall(pat, txt))
        misses = []
        if n < floor:
            misses.append("the statement appears %d time(s); the clean tree had %d, so "
                          "at least one place that said it no longer does" % (n, floor))
        for a in absent:
            k = len(re.findall(a, txt))
            if k:
                misses.append("the superseded form %r is printed %d time(s)" % (a, k))
        if misses:
            bad.append((name, misses, why))
            if verbose:
                print("  FAIL  %s" % name)
                for m in misses:
                    print("          %s" % m)
                print("          why: %s" % why)
        elif verbose:
            print("  ok    %-62s %3d" % (name, n))
    return bad


def selftest(txt, checks):
    """Two controls per check, and BOTH must fire: delete one occurrence of the correct
    statement, and paste the superseded statement back in."""
    fired, total = 0, 0
    for i, (name, (pat, floor), absent, poison, why) in enumerate(checks):
        total += 1
        if not run(re.sub(pat, "@@", txt, count=1), [checks[i]], verbose=False):
            print("  CONTROL FAILED: deleting one occurrence does not make %r fire. "
                  "Its floor is %d and the tree has more, so the floor is stale."
                  % (name, floor))
            return False
        fired += 1
        total += 1
        if not run(txt + " " + poison, [checks[i]], verbose=False):
            print("  CONTROL FAILED: pasting %r back in does not make %r fire."
                  % (poison, name))
            return False
        fired += 1
    print("  %d of %d synthetic controls fired (two per check)." % (fired, total))
    return fired == total


def main():
    print("check_coherence.py -- does the manuscript still agree with itself?")
    print("repository root: %s" % ROOT)
    txt = flat(body())
    if len(txt) < 100000:
        print("RESULT: ABORT -- the expanded body is only %d characters. The expander "
              "reached nothing and a green run would mean nothing." % len(txt))
        return 2
    checks, lo, hi, nsizes = build_checks(txt)
    print("body: %d characters of composed source" % len(txt))
    print("re-derived from data/c3_frontier/published_fraction/: %d sizes, "
          "vacuity %s to %s over the twenty cells" % (nsizes, lo, hi))
    print()
    print("SELF-TEST")
    if not selftest(txt, checks):
        print("RESULT: ABORT -- a control did not fire; this file certifies nothing.")
        return 2
    print()
    print("CHECKS (%d)" % len(checks))
    bad = run(txt, checks)
    print()
    if bad:
        print("RESULT: FAIL -- %d cross-section claim(s) of the manuscript contradict "
              "another section or the deposit." % len(bad))
        return 1
    print("RESULT: PASS -- every cross-section claim checked here is stated the same "
          "way in every section that states it.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
