# -*- coding: utf-8 -*-
r"""Guardian for a LENGTH-REDUCTION pass: prove that only prose was removed.

WHY THIS FILE EXISTS
--------------------
On 2026-09-19 the manuscript measured 31 668 words of prose against v2's 12 059 -- a
factor of 2.63 -- and the decision was taken to cut.  A cut is the single most dangerous
edit this project can make, for a reason the other four guardians cannot see:

    verify.py compares a printed number against a deposited file.  check_coherence.py
    compares a claim against a claim.  NEITHER OF THEM NOTICES A DELETION.  If a
    concession is removed outright, no printed number contradicts a file and no section
    contradicts another section: the manuscript is merely quieter, and every guardian
    stays green.

A cut that removes a measured number, a concession, a declaration of limit, a priority
attribution or a retraction turns an honest paper into a dishonest one WITHOUT tripping
a single existing check.  That is the class this file guards.

THE RULE IT ENFORCES
--------------------
    Repetition may go.  Facts may not.

Concretely, against a git baseline:

  * NUMBERS.  Every distinctive numeric literal (decimal, >= 3 digits, or an exponent)
    present in the baseline prose must still be present SOMEWHERE in the manuscript.
    Its COUNT may fall -- that is the entire point of the pass, a number stated five
    times and now stated twice is a success -- but it may not reach zero.
  * CONCESSIONS.  Every sentence in the baseline carrying a concession marker must still
    have a match in the trimmed manuscript, compared on content words rather than
    characters, so that terser phrasing passes and deletion does not.  Rewordings are
    reported by name for the author to read; only disappearances are fatal.

WHAT IT DOES NOT DO
-------------------
It does not judge whether the trim was worth making, it does not check physics, and it
cannot tell an honest rewording from a dishonest one -- it can only guarantee that every
rewording is SEEN.  The reading of the reworded list is the author's job and is not
delegable to this file.

SELF-TEST FIRST
---------------
A PASS is worth nothing without controls that fire.  Before reporting anything, this
file deletes a number and a concession from a copy of the trimmed text and requires both
to fail.  If either control does not fire, it exits 2 and reports NO result.
"""
from __future__ import print_function

import io
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
PAPER = os.path.join(REPO, "paper")

BASELINE = sys.argv[1] if len(sys.argv) > 1 else "bc9670c"

# Files whose prose is under this guardian.  bibliography.tex is excluded: it is
# reference data, not prose, and its numbers are years, volumes and pages.
SKIP = ("bibliography.tex",)

# ---------------------------------------------------------------------------
#  Markers.  A sentence carrying any of these is treated as load-bearing honesty
#  and may not vanish.  The list is deliberately over-inclusive: a false positive
#  costs one line of review, a false negative costs the paper's credibility.
# ---------------------------------------------------------------------------
MARCAS = (
    # -- negations and inability
    r"\bwe do not\b", r"\bwe cannot\b", r"\bwe have not\b", r"\bwe did not\b",
    r"\bdoes not\b", r"\bdo not\b", r"\bcannot\b", r"\bno bound\b",
    r"\bnot certif", r"\buncertified\b", r"\bis not\b", r"\bare not\b",
    r"\bnever\b", r"\bnone of\b", r"\bnothing\b", r"\bno evidence\b",
    # -- vacuity, the paper's central concession
    r"\bvacuous\b", r"\bvacuity\b", r"\btrivial\b", r"\bnon-?vacuous\b",
    # -- explicit limitation language
    r"\blimitation\b", r"\bcaveat\b", r"\bremains open\b", r"\bopen question\b",
    r"\bwe stress\b", r"\bwe emphasi", r"\bhonest\b", r"\bhonestly\b",
    r"\bat present\b", r"\bfor now\b", r"\bso far\b", r"\byet to\b",
    # -- failure, loss, absence of advantage
    r"\bfails?\b", r"\bfailure\b", r"\bloses\b", r"\blost to\b", r"\bweaker\b",
    r"\bno advantage\b", r"\bnot demonstrat", r"\bnot establish",
    r"\bworse than\b", r"\boutperform", r"\bunderperform",
    # -- limit declarations
    r"\bat most\b", r"\bno more than\b", r"\bupper bound\b", r"\bonly if\b",
    r"\bonly when\b", r"\brestricted to\b", r"\bconfined to\b",
    # -- priority attribution
    #    THIS BLOCK WAS WIDENED ON 2026-09-19 AFTER IT FAILED.  The first version
    #    carried only "prior work" and "earlier work", and the paragraph it therefore
    #    ranked as the MOST cuttable in the whole manuscript -- 129 words, no number,
    #    no marker -- was this, in sec_3_body:
    #        "Nor is ``a posteriori'' ours to claim, and we name the prior art rather
    #         than wait to be shown it. [...] That is the framing of this section,
    #         established earlier and by others, and the priority is theirs."
    #    which is a priority attribution naming four prior works.  The manuscript says
    #    "prior art", not "prior work"; "established earlier and by others", not
    #    "earlier work".  A cut driven by the old list would have deleted a priority
    #    concession first.  Phrase-matching honesty is brittle; widen on every miss.
    r"\bprior (?:work|art|literature)\b", r"\bearlier work\b",
    r"\bindependent", r"\bconcurrent", r"\bpriorit", r"\bours to claim\b",
    r"\bnot ours\b", r"\bestablished earlier\b", r"\bby others\b",
    r"\bbefore (?:this|version|our)\b", r"\bword for word\b",
    r"\bwe name the\b", r"\bthe(?:ir)?s\b(?=[^.]{0,20}$)",
    r"\bfirst (?:demonstrat|report|prov|observ)", r"\bfollowing \\cite",
    r"\bdue to \\cite", r"\bas shown in \\cite", r"\bcredit\b",
    r"\battribut", r"\bwho (?:first|already)\b", r"\bpredates?\b",
    r"\banticipated by\b", r"\bset out\b", r"\bsix weeks before\b",
    # -- retraction / correction of the record
    r"\bsupersed", r"\bcorrects?\b", r"\bearlier version\b", r"\bretract",
    r"\bpreviously (?:claimed|stated|reported)\b", r"\ban error\b",
    r"\bincorrect", r"\berroneous",
)
RE_MARCAS = re.compile("|".join(MARCAS), re.I)

# Words too common to identify a sentence by.
VACIAS = set("""a an the and or but of to in on for with by is are was were be been being
that this these those it its as at from we our us they their which who whom whose not no
than then so such if when where while have has had do does did can could may might must
shall should will would there here one two both each any all some more most less least
same other another into over under about between within without""".split())


def sin_comentarios(t):
    """Drop LaTeX comments without being fooled by an escaped percent sign."""
    t = re.sub(r"(?m)^\s*%.*$", "", t)
    t = re.sub(r"(?<!\\)%.*", "", t)
    return t


def texto_de(t):
    """Reduce a .tex body to something sentence-splittable, keeping numbers."""
    t = sin_comentarios(t)
    # Arguments of cross-reference commands are identifiers, not prose or data.
    t = re.sub(r"\\(?:ref|eqref|label|cite[a-zA-Z]*|input|include|bibliography)"
               r"\s*\{[^{}]*\}", " ", t)
    return t


# ---------------------------------------------------------------------------
#  Numbers
# ---------------------------------------------------------------------------
RE_NUM = re.compile(r"(?<![A-Za-z0-9._-])(\d[\d,]*(?:\.\d+)?)(?![0-9])")


def numeros(t):
    """Multiset of numeric literals, thousands separators normalised away."""
    cuenta = {}
    for m in RE_NUM.finditer(texto_de(t)):
        n = m.group(1).replace(",", "")
        n = n.rstrip(".")
        if not n:
            continue
        # Normalise 2.40 and 2.4 to the SAME key?  No -- trailing zeros are
        # significant figures in this manuscript and losing one is a real loss.
        cuenta[n] = cuenta.get(n, 0) + 1
    return cuenta


def distintivo(n):
    """True if losing this number entirely would be a loss of measured content.

    Bare integers below 100 ("two sizes", "Sec. 5", "Eq. 12") are structural and
    their disappearance is normal in a trim; decimals, long integers and anything
    with significant figures are data.
    """
    if "." in n:
        return True
    return len(n.lstrip("0")) >= 3


# ---------------------------------------------------------------------------
#  Concessions
# ---------------------------------------------------------------------------
def frases(t):
    """Split into sentences, with math and commands reduced to placeholders."""
    t = texto_de(t)
    # Displayed math is not prose; inline math is kept as a token so that a
    # sentence built around a formula still has a stable signature.
    t = re.sub(r"\\begin\{(equation|align|gather|eqnarray|widetext|figure|table|"
               r"tabular|tikzpicture|axis)\*?\}.*?\\end\{\1\*?\}", " MATHBLOCK ", t,
               flags=re.S)
    t = re.sub(r"\$\$.*?\$\$", " MATH ", t, flags=re.S)
    t = re.sub(r"\$[^$]*\$", " MATH ", t)
    t = re.sub(r"\\[a-zA-Z]+\*?", " ", t)
    t = re.sub(r"[{}&~^_\\]", " ", t)
    t = re.sub(r"\s+", " ", t)
    # Protect the abbreviations this manuscript actually uses before splitting.
    for a in ("Sec.", "Eq.", "Eqs.", "Fig.", "Figs.", "Ref.", "Refs.", "App.",
              "Tab.", "cf.", "e.g.", "i.e.", "vs.", "Phys.", "Rev.", "Lett.",
              "approx.", "no.", "Nos.", "et al."):
        t = t.replace(a, a.replace(".", "\x00"))
    piezas = [p.replace("\x00", ".").strip()
              for p in re.split(r"(?<=[.!?])\s+", t) if p.strip()]
    # A very short sentence cannot carry a signature of its own, and dropping it
    # loses real content: "It is not." following a claim IS the refutation of that
    # claim, and sec_7 carries exactly that construction.  Attach short pieces to
    # their predecessor rather than discarding them.
    fundidas = []
    for p in piezas:
        if fundidas and len(p) < 25:
            fundidas[-1] = fundidas[-1] + " " + p
        else:
            fundidas.append(p)
    return fundidas


def firma(f):
    """Content-word signature: what the sentence SAYS, not how it is phrased."""
    pal = re.findall(r"[A-Za-z][A-Za-z'-]+|\d[\d.]*", f.lower())
    return frozenset(w for w in pal if w not in VACIAS and len(w) > 1)


def concesiones(t):
    out = []
    for f in frases(t):
        if len(f) < 25:
            continue
        if RE_MARCAS.search(f):
            s = firma(f)
            if len(s) >= 5:          # too short to identify reliably
                out.append((f, s))
    return out


def jaccard(a, b):
    u = len(a | b)
    return (len(a & b) / float(u)) if u else 0.0


# ---------------------------------------------------------------------------
#  Corpus assembly
# ---------------------------------------------------------------------------
def ficheros():
    return sorted(f for f in os.listdir(PAPER)
                  if f.endswith(".tex") and f not in SKIP)


def corpus_base():
    """The baseline text, read from git so no stale copy on disk can be used."""
    out = {}
    listado = subprocess.check_output(
        ["git", "-C", REPO, "ls-tree", "--name-only", BASELINE + ":paper"],
        universal_newlines=True).split()
    for nombre in sorted(listado):
        if not nombre.endswith(".tex") or nombre in SKIP:
            continue
        out[nombre] = subprocess.check_output(
            ["git", "-C", REPO, "show", "%s:paper/%s" % (BASELINE, nombre)],
            universal_newlines=True)
    return out


def corpus_hoy():
    out = {}
    for nombre in ficheros():
        out[nombre] = io.open(os.path.join(PAPER, nombre),
                              encoding="utf-8", errors="replace").read()
    return out


def fundir(corpus):
    return "\n".join(corpus[k] for k in sorted(corpus))


# ---------------------------------------------------------------------------
#  The two checks
# ---------------------------------------------------------------------------
def revisa_numeros(antes, ahora):
    na, nh = numeros(antes), numeros(ahora)
    muertos_d, muertos_c, caidos = [], [], []
    for n, c in sorted(na.items(), key=lambda kv: -kv[1]):
        c2 = nh.get(n, 0)
        if c2 == 0:
            (muertos_d if distintivo(n) else muertos_c).append((n, c))
        elif c2 < c:
            caidos.append((n, c, c2))
    return muertos_d, muertos_c, caidos


def revisa_concesiones(antes, ahora, umbral=0.60):
    ca = concesiones(antes)
    ch = concesiones(ahora)
    firmas_h = [s for _, s in ch]
    perdidas, reescritas = [], []
    for f, s in ca:
        mejor, cual = 0.0, None
        for i, s2 in enumerate(firmas_h):
            j = jaccard(s, s2)
            if j > mejor:
                mejor, cual = j, i
        if mejor >= 0.999:
            continue                       # survives verbatim
        if mejor >= umbral:
            reescritas.append((f, ch[cual][0], mejor))
        else:
            perdidas.append((f, mejor))
    return perdidas, reescritas, len(ca), len(ch)


# ---------------------------------------------------------------------------
#  Controls.  Nothing is reported until both of these fire.
# ---------------------------------------------------------------------------
def controles(antes, ahora):
    fallos = []

    # (1) delete one distinctive number that survives today -> numbers must fail
    na, nh = numeros(antes), numeros(ahora)
    cebo = None
    for n in sorted(nh, key=lambda x: -nh[x]):
        if distintivo(n) and na.get(n, 0) > 0:
            cebo = n
            break
    if cebo is None:
        fallos.append("no distinctive number available to build the control on")
    else:
        mutilado = re.sub(r"(?<![A-Za-z0-9._-])" + re.escape(cebo) + r"(?![0-9])",
                          " QQQ ", ahora)
        md, _, _ = revisa_numeros(antes, mutilado)
        if not any(n == cebo for n, _ in md):
            fallos.append("control 1 did not fire: deleting %s was not detected" % cebo)

    # (2) delete one concession sentence -> concessions must fail
    ch = concesiones(ahora)
    if not ch:
        fallos.append("no concession available to build the control on")
    else:
        # pick one that also exists in the baseline, else the control is vacuous
        ca = concesiones(antes)
        objetivo = None
        for f, s in ch:
            if any(jaccard(s, s2) >= 0.999 for _, s2 in ca):
                objetivo = f
                break
        if objetivo is None:
            fallos.append("no shared concession available to build the control on")
        else:
            # remove the literal sentence from the raw text by its rarest words
            raras = sorted(firma(objetivo), key=len, reverse=True)[:4]
            mutilado = ahora
            for f_raw in re.split(r"(?<=[.!?])\s+", mutilado):
                if all(r in f_raw.lower() for r in raras):
                    mutilado = mutilado.replace(f_raw, " ")
                    break
            p, _, _, _ = revisa_concesiones(antes, mutilado)
            if not p:
                fallos.append("control 2 did not fire: deleting a concession "
                              "was not detected")
    return fallos


# ---------------------------------------------------------------------------
def main():
    print("=" * 78)
    print("  check_trim.py -- only prose may leave.  Baseline: %s" % BASELINE)
    print("=" * 78)

    try:
        base = corpus_base()
    except subprocess.CalledProcessError:
        print("  cannot read baseline %s from git." % BASELINE)
        return 2
    hoy = corpus_hoy()
    antes, ahora = fundir(base), fundir(hoy)

    print("  baseline : %2d files, %7d chars" % (len(base), len(antes)))
    print("  worktree : %2d files, %7d chars" % (len(hoy), len(ahora)))
    print()

    print("  -- controls --")
    fallos = controles(antes, ahora)
    if fallos:
        for f in fallos:
            print("     X %s" % f)
        print("\n  CONTROLS DID NOT FIRE.  No result is reported.")
        return 2
    print("     both controls fire.")
    print()

    md, mc, caidos = revisa_numeros(antes, ahora)
    perdidas, reescritas, nca, nch = revisa_concesiones(antes, ahora)

    print("  -- numbers --")
    print("     distinctive literals in baseline : %d"
          % len([n for n in numeros(antes) if distintivo(n)]))
    print("     repetitions removed (count fell) : %d" % len(caidos))
    if caidos[:8]:
        for n, c, c2 in sorted(caidos, key=lambda x: x[1] - x[2], reverse=True)[:8]:
            print("        %-12s %2d -> %2d" % (n, c, c2))
    print("     DISTINCTIVE LITERALS LOST        : %d" % len(md))
    for n, c in md[:40]:
        print("        X %-14s was stated %d time(s), now absent" % (n, c))
    print("     common integers lost (report)    : %d" % len(mc))
    for n, c in mc[:10]:
        print("        . %-14s (%d)" % (n, c))
    print()

    print("  -- concessions --")
    print("     marked sentences: %d baseline -> %d now" % (nca, nch))
    print("     reworded (READ THESE): %d" % len(reescritas))
    for f, g, j in reescritas[:40]:
        print("        ~ %.2f" % j)
        print("          was: %s" % f[:150])
        print("          now: %s" % g[:150])
    print("     CONCESSIONS LOST: %d" % len(perdidas))
    for f, j in perdidas[:40]:
        print("        X (best match %.2f) %s" % (j, f[:180]))
    print()

    mal = len(md) + len(perdidas)
    print("=" * 78)
    if mal:
        print("  FAIL: %d distinctive number(s) and %d concession(s) left the paper."
              % (len(md), len(perdidas)))
        print("  Repetition may go.  Facts may not.")
        return 1
    print("  PASS: every distinctive number and every concession still present.")
    print("  %d repetitions removed.  %d rewordings to read by hand."
          % (len(caidos), len(reescritas)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
