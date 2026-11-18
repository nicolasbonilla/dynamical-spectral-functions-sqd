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
         # 2026-09-30: floor 5 -> 4 DOWNWARD, reframe (skeptic V1): Sec. IX C no longer
         # restates the range; its home is Sec. V B (Table published), and the abstract,
         # Sec. I and Sec. X still state it.  The range itself is unchanged.
         (rng, 4),
         [r"\$2\.8\$\s*(?:--|to)\s*\$4\.9\$", r"2\.8 to 4\.9"],
         "vacuous by factors of $2.8$ to $4.9$",
         "the abstract, Sec. I, Sec. V and Sec. X each state the range over "
         "which the certificate is vacuous. It is %s-%s over the twenty cells of the "
         "deposited published-fraction grid, and the superseded 2.8-4.9 must appear "
         "nowhere." % (lo_s, hi_s)),

        ("the number of published sizes agrees with the deposited grid",
         (r"at all %s sizes|at each of the %s sizes evaluated|at the %s sizes measured"
          r"|%s sizes up to a sector|at %s system sizes|all five sizes and all four"
          % (word, word, word, word, word), 5),
         # 2026-09-30: floor 6 -> 5 DOWNWARD, reframe (skeptic V1): the Sec. IX C sentence
         # 'at the five sizes measured' went with the restated vacuity range (see above).
         # 2026-09-25: floor 8 -> 6, DOWNWARD: the rewritten introduction and conclusion state
         # the five sizes fewer times; no section counts four (the forbidden forms below).
         [r"at all four sizes", r"at each of the four sizes evaluated",
          r"at the four sizes measured", r"four sizes up to a sector",
          r"at four system sizes"],
         "at all four sizes",
         "the deposited grid carries %d sizes, so no section may still count four"
         % nsizes),

        ("the Krylov-support claim is the measured one, not an open question",
         # 2026-09-25: the support is every determinant the reflection about the seed site
         # allows (the "whole sector" counts were round-off in the forbidden sector),
         # measured at L = 6 to 14 with z_suppK_sym.py.  Eight places state it (the abstract
         # joined them the same day).
         # 2026-09-25: floor 8 -> 7. DOWNWARD, intro rewrite: the new introduction states the symmetry-allowed support once instead of twice.
         (r"every symmetry-allowed determinant|every determinant (?:that )?the reflection about the seed site allows", 7),
         [r"measured at two sizes and we state it for",
          r"Two sizes are what we",
          r"whether it persists at \$L=10\$ and \$L=12\$ is open"],
         "Two sizes are what we measured",
         "Sec. V measures the coordinate support of K(H,phi) at L = 8, 10, 12 and 14. "
         "The introduction and the conclusion must not still offer it as open, and the "
         "cost of settling it must not still be quoted as unspent."),

        ("the hardware shot budget is one number in all four places",
         # 2026-09-25 (PRA cut): the short Sec. IX A states it once more.
         # 2026-09-25: floor 4 -> 7. 2026-09-26: Sec. IX A and Sec. S8 now also print the budget when disclosing the reversed bit order of the L=6 run
         (r"3\.5\\times10\^\{5\}", 7),
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
         # 2026-09-25: floor 2 -> 3. PRA cut: Sec. V summary in the main text states it once, the full text moved to the SM.
         (r"a drift(?:,| and) not scatter", 3),
         [r"uniformly over \$L=6\$--\$12\$"],
         "uniformly over $L=6$--$12$",
         "Sec. V spends thirty lines showing that the factor of three is an ordered "
         "drift and not scatter. 'Uniformly' is the opposite claim, and it stood in "
         "the conclusion."),

        ("the sector fraction is named as Sec. VII's own rule requires",
         # Floor 17 -> 18 on 2026-09-23. The reframing added one occurrence in
         # sec_1.tex (the obstruction paragraph now names the fraction at which
         # the certificate is evaluated), so the floor went stale in the safe
         # direction and the SELF-TEST caught it: deleting one occurrence left
         # 17, which still cleared a floor of 17, so the control could not fire
         # and the file refused to certify anything. Raise the floor when the
         # tree gains an occurrence; never lower it to make a run go green.
         # Floor 18 -> 19 on 2026-09-25: the corrected Sec. V.E separates the
         # measured vacuity "at the published fractions" from the structural
         # statement, and names the fraction twice where it once named it once
         # (counted with body()+flat(), 18 on the morning tree, 19 after).
         # Floor 19 -> 20 on 2026-09-25: Sec. I now says which reconstructions are
         # uncertified -- those "at the published fractions" -- and which one is
         # not (Fig. akwsampled, at 85%). body()+flat(): 19 before, 20 after.
         # Floor 20 -> 19 on 2026-09-25, the one DOWNWARD move in this file, and why:
         # two of the pinned occurrences (sec_3_1 and the Fig. thm1iii caption) called
         # 0.85 -- the fraction of Fig. akwsampled, which IS certified -- "the published
         # fraction", which made "vacuous at every published fraction" literally false.
         # They were misuses, confirmed by an adversarial review, and now read "the
         # fraction 0.85 of Fig. akwsampled".  body()+flat(): 19 after the correction.
         # 19 -> 20 later on 2026-09-25: the conclusion and Sec. V.E now say where the
         # published fractions lie relative to full captured weight (body()+flat(): 20).
         # 2026-09-25, last change of the day: the abstract, Sec. IX and Sec. X now say
         # "operating fraction(s) of the resource scan", which is the precise name (the
         # 85% of Fig. akwsampled is published too, and certified). Both names count;
         # body()+flat(): 19 "published fraction" + 4 "operating fraction ... scan" = 23.
         # 2026-09-25: floor 23 -> 26. PRA cut: the Sec. V summary names the published and operating fractions again; the verbatim text is in the SM.
         # 2026-09-25: floor 26 -> 27. PRA cut: the Sec. VII summary names the published fractions once more; the verbatim text is in the SM.
         # 2026-09-25: floor 27 -> 24. DOWNWARD, final referee pass: three 'published fraction' claims of vacuity also covered Fig. akwsampled (published at 85% and certified); they now say 'operating fraction(s)' or 'those fractions'.
         # 2026-09-25: floor 24 -> 23. intro rewrite: the new introduction names the operating fraction once.
         # 2026-09-25: floor 23 -> 21. conclusion rewrite: operating fraction named once more.
         # 2026-09-30: floor 21 -> 20. DOWNWARD, reframe (skeptic V1): Sec. IX C no longer restates the vacuity at 'every operating fraction of the resource scan'; Sec. V B, the abstract, Sec. I and Sec. X still do.
         # 2026-10-01: floor 20 -> 21. UPWARD, closing review of the prose pass: Sec. I now defines the operating fractions of the resource scan in words before stating the vacuity there, one more occurrence (body()+flat(): 21).
         (r"published fraction|operating fractions? of the resource scan", 21),
         [r"sampled fraction"],
         "at every sampled fraction",
         "Sec. VII states in the printed text that describing the T->infinity curve as "
         "SAMPLED without a shot count is to mislabel it, and that it is labelled the "
         "published curve from there on. The abstract and Sec. X broke that rule."),

        ("the declared gap of Fig. 3(a,b) is closed, and said to be closed",
         # 2026-09-25: "has since been run" read as version history; the paper now says
         # "has also been run", which states the same fact.
         (r"has also been run at finite shots", 1),
         [r"Declared gap"],
         "Declared gap, see text",
         "a finite-shot run of the literal configuration of Fig. 3(a,b) now exists "
         "(5.2e6 shots per channel), so the provenance table and Sec. II B must not "
         "still declare it missing."),

        ("the data-availability statement quotes the guardian it actually has",
         # This literal has moved THREE times in three days, always for the same benign
         # reason -- adding a .py to src/ moves the README script-count check and the
         # assertion total with it: 9843 -> 9844 (check_trim.py), 9844 -> 9845
         # (abs_paste.py), 9845 -> 9846 (the tab:moments coverage gap).  The pin stays a
         # literal on purpose: its job is to REFUSE a silent change, and it has done that
         # every time.  When it fires, do not edit it by hand first -- run
         # scratchpad/propaga_cuentas.py, which re-runs verify.py and rewrites the eight
         # places that quote the total, and then bring this line to match.
         # 2026-09-28: 872 / 9,889 -> 874 / 9,891, two phrase guards added to verify.py
         # (hw.sm_scope_*, CLAIMS_RETRACTED).
         # 2026-09-28: 874 / 9,891 -> 874 / 11,487, the Fig. 3 frequency grid widened from
         # [-9t, 9t] (600 points) to [-9t, 17t] (866 points): verify.py checks every plotted
         # point of the six curves, +266 each (KNOWN_DISCREPANCIES.md section 31).
         # 2026-09-28 (late): 874 / 11,487 -> 934 / 11,674, verify.py section 9.6b recomputes
         # the finite-shot paragraph of the Fig. 3 caption (and its echoes in Secs. II, VII,
         # S1, S8) from data/fig3_finite_shot_T*.json, the outputs of the deposited
         # src/fig3_finite_shot.py (KNOWN_DISCREPANCIES.md section 31).
         # 2026-09-28 (late, closing the review): 934 / 11,674 -> 1212 / 14,495, verify.py
         # section 9.6c (the free-fermion supports, the spin-peak deviations and the Fig. S5 pool re-classified from the deposited certificate rows), the chi_max of the witness, the LF hash of the Fig. 3 generator and the L=14 re-ranking record
         # (KNOWN_DISCREPANCIES.md section 33); its Fig. S5 coverage gap is closed (xfail 5 -> 4).
         # 2026-09-30: 1212 / 14,495 -> 1214 / 14,497, verify.py section 9.3 now checks every
         # printed exact permutation p (Sec. I, Sec. VI, Fig. 6 caption: three occurrences,
         # each against exact_perm_p_one_sided_F1 at 1e-15) instead of the first one at 5e-15.
         # 2026-09-30 (final review): 1214 / 14,497 -> 1280 / 14,563, verify.py sections 9.6b
         # (Fig. 3 caption panel (b), window and branch thresholds), 9.9b (the printed noise
         # values) and 9.9c (every restated Fig. 3 certificate number).
         # 2026-11-18 (deposit 3.5.0, the manuscript prepared for JCP): 1280 / 14,563 -> 1269 / 16,565.  The
         # statement is in sec_10.tex; verify.py section (12) runs the 21 checks of src/verify_jcp.py, and the
         # first-order checks of 9.6c and 9.9c they replace left section 9 (docs/KNOWN_DISCREPANCIES.md section 36).
         # The counts of 3.4.1 are forbidden below.
         (r"\$1269\$ checks over \$16\\,565\$ numeric assertions", 1),
         # The forbidden phrase used to be the bare "four known open defects", which was
         # unambiguous while the live sentence read "five".  On 2026-09-21 the live
         # sentence became "four known open defects, and two declared gaps in its own
         # coverage" -- a DIFFERENT and TRUE statement that happens to contain the
         # banned substring, so the check fired on correct text.  The ban is therefore
         # narrowed to the superseded SENTENCE, not the phrase inside it.  Narrowing a
         # ban is dangerous and is done here only because the wider form now has a
         # legitimate match; the superseded claim is still caught in full.
         [r"\$807\$ checks", r"\$7\\,454\$ numeric assertions",
          r"four known open defects are printed by name",
          r"\$1280\$ checks", r"\$14\\,563\$ numeric assertions"],
         "$807$ checks over $7\\,454$ numeric assertions",
         "Sec. IX F tells the reader what `python src/verify.py` prints. Run it and "
         "compare. If verify.py's counts have moved, this check is the reminder that "
         "Sec. IX F, README.md, docs/DATA_AVAILABILITY.md, docs/REPRODUCE.md and "
         "docs/KNOWN_DISCREPANCIES.md all move with it."),

        ("the largest sector quoted is the largest sector evaluated",
         # 10 -> 8 on 2026-09-25: two of the ten were the unprojected power-iteration counts
         # "10,306,296 determinants" and "10,306,063 of 10,306,296", round-off in the
         # reflection-forbidden sector, now removed.  The abstract, Secs. I, V, VII, IX, X
         # and both tables still quote the largest sector (body()+flat(): 8).
         (r"10\\,306\\,296", 8),
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


# ---------------------------------------------------------------------------
#  2026-09-30 (final numeric review NR-8).  Three pooled floors above were lowered that day
#  (5 -> 4, 6 -> 5, 21 -> 20) when Sec. IX C stopped restating the vacuity.  A pooled floor
#  lets one place go silently while another gains a copy, so each statement is now also
#  required IN EACH PLACE that carries it, per source file (the abstract is cut out of
#  main.tex).  The per-file floors are the counts measured on the tree of 2026-09-30.  Each
#  check has its own control: deleting its occurrences in that one file must make it fire.
# ---------------------------------------------------------------------------
def part(where):
    if where == "abstract":
        m = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", flat(body("main.tex")))
        return m.group(1) if m else ""
    return flat(body(where))


def located_checks(rng, word):
    sizes = (r"at all %s sizes|at each of the %s sizes evaluated|at the %s sizes measured"
             r"|%s sizes up to a sector|at %s system sizes|all five sizes and all four"
             % (word, word, word, word, word))
    opf = r"published fraction|operating fractions? of the resource scan"
    out = []
    for where in ("abstract", "sec_1.tex", "sec_5_body.tex", "sec_10.tex"):
        out.append(("the vacuity range is stated in %s" % where, where, rng, 1))
    for where in ("abstract", "sec_1.tex", "sec_5_body.tex", "sec_9.tex", "sec_10.tex"):
        out.append(("the number of published sizes is stated in %s" % where, where, sizes, 1))
    # 2026-10-01: sec_1.tex 1 -> 2, UPWARD (closing review of the prose pass): Sec. I now
    # defines the operating fractions of the resource scan before stating the vacuity there.
    for where, n in (("abstract", 1), ("sec_1.tex", 2), ("sec_5_body.tex", 5), ("sec_7.tex", 1),
                     ("sec_10.tex", 2), ("sm_frontier.tex", 8), ("sm_prices.tex", 1)):
        out.append(("the operating/published fraction is named in %s" % where, where, opf, n))
    return out


def run_located(checks, parts, verbose=True):
    bad = []
    for name, where, pat, floor in checks:
        n = len(re.findall(pat, parts[where]))
        if n < floor:
            bad.append(name)
            if verbose:
                print("  FAIL  %s: %d time(s), the tree of 2026-09-30 had %d" % (name, n, floor))
        elif verbose:
            print("  ok    %-62s %3d" % (name, n))
    return bad


def selftest_located(checks, parts):
    for c in checks:
        name, where, pat, floor = c
        cut = dict(parts)
        cut[where] = re.sub(pat, "@@", parts[where])
        if not run_located([c], cut, verbose=False):
            print("  CONTROL FAILED: deleting %r from %s does not make it fire." % (pat[:40], where))
            return False
    print("  %d of %d per-place controls fired." % (len(checks), len(checks)))
    return True


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
    word = {4: "four", 5: "five", 6: "six"}[nsizes]
    rng = re.escape("$%s$" % lo) + r".{0,40}" + re.escape("$%s$" % hi)
    loc = located_checks(rng, word)
    parts = dict((w, part(w)) for w in set(c[1] for c in loc))
    print("PER-PLACE CHECKS (%d)" % len(loc))
    if not selftest_located(loc, parts):
        print("RESULT: ABORT -- a per-place control did not fire; this file certifies nothing.")
        return 2
    bad += [(n, [], "") for n in run_located(loc, parts)]
    print()
    if bad:
        print("RESULT: FAIL -- %d cross-section claim(s) of the manuscript contradict "
              "another section or the deposit." % len(bad))
        return 1
    print("RESULT: PASS -- every cross-section claim checked here is stated the same "
          "way in every section that states it.")
    return 0


# =====================================================================================================================
#  DEPOSIT 3.5.0 (2026-11-18): the registrations of the manuscript prepared for The Journal of Chemical Physics.
#
#  The JCP package checked its sources with a wrapper (check_coherence_jcp.py) that imported this file unchanged,
#  pointed it at the JCP sources and re-registered what the JCP text changed, each re-registration first asserting
#  the 3.4.1 value it replaces.  In 3.5.0 paper/ IS the JCP text, so the wrapper is merged here, word for word in its
#  logic: the functions above are the checks of 3.4.1, and the three wrapped below are what main() now calls.  Each
#  change, with its reason (the dated notes are the wrapper's):
#
#  1. 2026-10-28 (device runs, on the author's instruction one line of Sec. S8):
#       - "the hardware shot budget is one number in all four places": floor 7 -> 1; its absent-list is WIDENED with
#         the absent-pattern of the check retired below;
#       - "the Fig. 14 caption carries the same total": retired (the figure is not printed); its superseded form
#         stays forbidden through the check above.
#  2. 2026-11-08 (results-first abstract): the abstract no longer prints the vacuity factors or the phrase
#     "operating fraction(s) of the resource scan"; both are still stated in Secs. I, V and X, which keep their own
#     per-place floors.  So "the vacuity range is stated in abstract" and "the operating/published fraction is named
#     in abstract" are retired; pooled floors move DOWNWARD by the count of what left (vacuity range 4 -> 3;
#     operating / published fraction 21 -> 16 in all, with items 3 and 4), pooled floors that went UP are raised to
#     the measured count (published sizes 5 -> 6, largest sector 8 -> 9); three NEW per-place pins hold the
#     abstract to its own content (the vacuity in words, the qualifier "infinite-shot Born-ranked", the hypothesis
#     "Given the exact ground state"); two NEW pooled checks (each with its two controls) forbid two overclaims
#     ("rule(s) out any bound indexed", "calibrated").
#  3. 2026-11-08, pass 7: Sec. S4 names the operating fraction 5 times, not 8 (the three others were in the
#     shot-frontier table and its paragraph, withdrawn as unsourced, docs/JCP_WITHDRAWN_NUMBERS.md): that
#     per-place floor 8 -> 5.
#  4. 2026-11-14 (closing reading): the conclusion's two paragraphs that opened with the vacuity are one; the
#     operating fraction of the resource scan is named once there, so the sec_10.tex per-place floor 2 -> 1
#     (DOWNWARD, one duplicate fewer).
#  5. 2026-11-16 (Theorem 1 gains the second-order branch; rules D7, D17, docs/KNOWN_DISCREPANCIES.md section 36): a
#     vacuity factor quoted above the trivial bound is the SMALLER of the two leakage branches.  vacuity_range() is
#     wrapped: the 3.4.1 values are asserted first (1.66, 8.61), then each cell takes min(first-order, second-order),
#     the second-order factor (sqrt(1-w)+Lambda^/eta)^2/(1+w) from the same stored fields.  The range becomes
#     1.38--8.61; the superseded 1.66--8.61 joins the forbidden forms of that check.
# =====================================================================================================================
_V4_VACUITY_RANGE, _V4_BUILD_CHECKS, _V4_LOCATED_CHECKS = vacuity_range, build_checks, located_checks

JCP_BUDGET = 'the hardware shot budget is one number in all four places'
JCP_CAPTION = 'the Fig. 14 caption carries the same total'
JCP_VACUITY = 'the vacuity range is the deposited one, everywhere it is stated'
JCP_SIZES = 'the number of published sizes agrees with the deposited grid'
JCP_OPF = "the sector fraction is named as Sec. VII's own rule requires"
JCP_SECTOR = 'the largest sector quoted is the largest sector evaluated'

# name -> (3.4.1 floor asserted, JCP floor, reason)
JCP_FLOORS = {
    JCP_VACUITY: (4, 3, 'abstract no longer prints 1.66--8.61; Secs. I, V, X do'),
    JCP_SIZES: (5, 6, 'Sec. VIII A Limitations adds "at all five sizes"'),
    JCP_OPF: (21, 16, 'abstract no longer names the operating fraction; Secs. I, V, VII, X, S4, S6 do; pass 7: '
                      'the shot-frontier table of S4, withdrawn as unsourced (docs/JCP_WITHDRAWN_NUMBERS.md), named '
                      'it three times; 2026-11-14 (closing reading): the conclusion states the vacuity at the '
                      'operating fractions once instead of twice (17 -> 16)'),
    JCP_SECTOR: (8, 9, 'the opening paragraph of Sec. I quotes the largest sector'),
}
JCP_VAC_V4 = ('1.66', '8.61')
JCP_FIRST_ORDER_RANGE = re.escape('$1.66$') + r'.{0,40}' + re.escape('$8.61$')

JCP_NEW_POOLED = [
    ('the captured weight forces a bound to its trivial value, it does not "rule it out"',
     (r'forces (?:the constant of )?any bound indexed on', 3),
     [r'rules? out any\s+bound indexed'],
     'an analytic counterexample and a counting argument rule out any bound indexed on it',
     'JCP 2026-11-08: the counterexample forces the constant of any bound indexed on the captured '
     'weight alone to its trivial value (Sec. III D); "rules out" is wider than that, and the fact '
     'check of the APS abstract on this work removed it.'),
    ('the crossing is measured, not "calibrated"',
     (r'The true error at which the bound becomes informative|Where the bound first improves on '
      r'the trivial one', 2),
     [r'\bcalibrated\b'],
     'Where the bound is informative, it is calibrated.',
     'JCP 2026-11-08: at the crossing the bound is 3e2--2.4e3 times the true error (Sec. III C; since 2026-11-16 '
     'the certificate, 1.4e2--7.8e2 times); '
     'what is reproducible is the true error at which it first improves on the trivial bound.'),
]

JCP_RETIRED_PLACES = ('the vacuity range is stated in abstract',
                      'the operating/published fraction is named in abstract',
                      'the operating/published fraction is named in sm_frontier.tex')
JCP_NEW_PLACES = [
    ('the abstract states the vacuity in words', 'abstract', r'no better than the trivial one', 1),
    ('the abstract keeps the qualifier "infinite-shot Born-ranked"', 'abstract',
     r'infinite-shot Born-ranked', 1),
    ('the abstract states the hypothesis of the bound', 'abstract', r'Given the exact ground state', 1),
]


def jcp_second_order_factors(root):
    """2026-11-16 (D7, D17): the second-order vacuity factor (sqrt(1-w)+x)^2/(1+w) of each of the twenty
    published-fraction cells, from the same stored fields vacuity_range() reads: PUB_rows.json (w_S, B, triv;
    x = B/(1+sqrt w) - sqrt(1-w), the stored first-order branch inverted) and L14_FR008.json at the deepest
    depth (w_S, LamHat_over_eta, bound_trivial)."""
    import io
    import json
    import math
    d = os.path.join(root, 'data', 'c3_frontier', 'published_fraction')
    out = {}
    for r in json.load(io.open(os.path.join(d, 'PUB_rows.json'), encoding='utf-8')):
        w = r['w_S']
        x = r['B'] / (1.0 + math.sqrt(w)) - math.sqrt(1.0 - w)
        out[(r['L'], r['eta'])] = (math.sqrt(1.0 - w) + x) ** 2 / r['triv']
    blk = json.load(io.open(os.path.join(d, 'L14_FR008.json'), encoding='utf-8'))[0]
    deepest = max(q['n'] for q in blk['rows'])
    for q in blk['rows']:
        if q['n'] == deepest:
            out[(blk['L'], q['eta'])] = (math.sqrt(1.0 - blk['w_S']) + q['LamHat_over_eta']) ** 2 / q['bound_trivial']
    return out


def _jcp_vacuity_range():
    vals = _V4_VACUITY_RANGE()
    got = ('%.2f' % min(vals.values()), '%.2f' % max(vals.values()))
    assert got == JCP_VAC_V4, ('the 3.4.1 vacuity range changed; re-register by hand', got)
    two = jcp_second_order_factors(ROOT)
    assert set(two) == set(vals), 'the second-order factors do not cover the cells of vacuity_range()'
    return dict((k, min(v, two[k])) for k, v in vals.items())


def _jcp_build_checks(txt):
    res = _V4_BUILD_CHECKS(txt)            # (checks, lo, hi, nsizes)
    checks = res[0]
    names = [c[0] for c in checks]
    assert JCP_BUDGET in names and JCP_CAPTION in names, 'the 3.4.1 checks changed; re-register by hand'
    for nm in JCP_FLOORS:
        assert nm in names, ('the 3.4.1 checks changed; re-register by hand', nm)
    cap = checks[names.index(JCP_CAPTION)]
    assert cap[1][1] == 1, cap[1]
    out = []
    for name, (pat, floor), absent, poison, why in checks:
        if name == JCP_BUDGET:
            assert floor == 7, floor
            out.append((name, (pat, 1), list(absent) + list(cap[2]), poison,
                        why + ' [JCP 2026-10-28: the device runs are one line of Sec. S8, which '
                        'prints 3.5e5 once; floor 7 -> 1; the retired caption check\'s '
                        'superseded form is forbidden here]'))
        elif name == JCP_CAPTION:
            continue
        elif name in JCP_FLOORS:
            rel, new, reason = JCP_FLOORS[name]
            assert floor == rel, (name, floor, rel)
            if name == JCP_VACUITY:
                # 2026-11-16 (D7, D17): the range is the smaller branch's (vacuity_range is re-registered below);
                # the superseded first-order range must appear nowhere
                absent = list(absent) + [JCP_FIRST_ORDER_RANGE]
                why = why + (' [JCP 2026-11-16: the smaller of the two leakage branches over the trivial '
                             'one; the first-order range 1.66-8.61 must appear nowhere]')
            out.append((name, (pat, new), absent, poison,
                        why + ' [JCP 2026-11-08: floor %d -> %d, %s]' % (rel, new, reason)))
        else:
            out.append((name, (pat, floor), absent, poison, why))
    out.extend(JCP_NEW_POOLED)
    return (out,) + tuple(res[1:])


def _jcp_located_checks(rng, word):
    loc = _V4_LOCATED_CHECKS(rng, word)
    have = [c[0] for c in loc]
    for r in JCP_RETIRED_PLACES:
        assert r in have, ('the 3.4.1 per-place checks changed; re-register by hand', r)
    # pass 7 (2026-11-08): S4 names the operating fraction 5 times, not 8: the three others were in the
    # shot-frontier table and its paragraph, withdrawn as unsourced; same pattern, floor 8 -> 5
    sf = [c for c in loc if c[0] == 'the operating/published fraction is named in sm_frontier.tex'][0]
    assert sf[3] == 8, sf
    extra = [(sf[0] + ' (floor 8 -> 5, pass 7)', sf[1], sf[2], 5)]
    # 2026-11-14 (closing reading): the two conclusion paragraphs that opened with the vacuity are merged into one
    # that states the crossing first and the vacuity once; sec_10.tex names the operating fraction of the resource
    # scan once (in the sentence check_trim keeps verbatim as a proof), floor 2 -> 1
    s10 = [c for c in loc if c[0] == 'the operating/published fraction is named in sec_10.tex'][0]
    assert s10[3] == 2, s10
    extra.append((s10[0] + ' (floor 2 -> 1, closing 2026-11-14)', s10[1], s10[2], 1))
    drop = JCP_RETIRED_PLACES + (s10[0],)
    return [c for c in loc if c[0] not in drop] + JCP_NEW_PLACES + extra


vacuity_range = _jcp_vacuity_range
build_checks = _jcp_build_checks
located_checks = _jcp_located_checks


if __name__ == "__main__":
    print('check_coherence: deposit 3.5.0, the JCP registrations of 2026-10-28, 2026-11-08, 2026-11-14 and '
          '2026-11-16 are merged (see the block above main\'s call)')
    sys.exit(main())
