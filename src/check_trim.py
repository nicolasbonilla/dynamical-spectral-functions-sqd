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
# ---------------------------------------------------------------------------
#  AUDITADO a mano el 2026-09-24, frase por frase, contra la base bc9670c.
#
#  Ese dia el resumen y la conclusion NO se recortaron: se REESCRIBIERON, y el
#  comparador de este fichero lee una reescritura como una perdida.  Cada entrada
#  de aqui abajo se comprobo a mano y solo se perdona SI SU PRUEBA SIGUE EN EL
#  PAPER HOY.  Si alguien borra el hecho, la prueba desaparece y el guardian
#  vuelve a saltar: esto no silencia la comprobacion, la convierte en
#  "la frase se movio, y aqui esta la prueba de que el hecho sigue".
#
#  NADA de lo auditado era un numero medido, una concesion, un limite ni una
#  atribucion.  Los tres literales:
#    843    -> el recuento de aserciones cambio al anadir comprobaciones (9846)
#    1.010  -> repeticion redondeada de (H2) en sec_5_body, sustituida por un
#    1.037     puntero; (H2) sigue entera, con los valores completos, junto al
#              teorema (sec_3_body) y en el apendice (sec_5_app)
# ---------------------------------------------------------------------------
# 2026-09-25 (tarde): la prueba se RESTAURA con su generador depositado (src/gflow_dequant.py)
# y los valores por semilla; los numeros v2 (39.24, 26.98, ...) se sustituyen por los
# regenerados y el acotamiento del t pareado por el t calculado. La prueba de que el
# hecho sigue en el paper es la frase que cita los estadisticos pareados.
DQ = "holds the paired statistics quoted here"
AUDITADO_NUM = {
    "68101214": "at five sizes up to a",   # la lista L=6,...,14 de la introduccion vieja
    # 2026-09-25: los numeros de la prueba de descuantizacion, retirada entera
    **{n: DQ for n in ("4.04", "3.30", "0.66", "39.24", "26.98", "22.94", "2.64", "3.96",
                       "2.281", "3.422", "2.776", "0.0625", "0.031", "12.26", "0.88", "1.25")},
    # 2026-09-25: recuentos caducados sustituidos por los del deposito, el 7.6e7 sin
    # script (y atado a la variable equivocada) sustituido por la ley exacta en eps^-2,
    # y la referencia a versiones anteriores trasladada al campo Comments de arXiv.
    "239":   "subspaces of the pooled sweep",
    "195":   "$378$ deposited subspaces",
    "7.6":   "the energy of the discarded component, does not enter",
    "2608.16436": "@comments:Theorem 1(iii) of v1-v2 does not hold",
    # 2026-09-25: la serie de dos niveles (App. B) integrada sobre toda la recta; los valores
    # anteriores cortaban las colas lorentzianas (src/moments_table.py).
    "1.9503": "1.9504", "1.9832": "1.9835", "1.9943": "1.9950", "1.9974": "1.9983",
    # 2026-09-25: supp(phi) contado por encima de 1e-12 (src/support_witness.py, W3/W8): en
    # L=8 la simetria anula 213 de los 2450 determinantes combinatorios, asi que el suelo es
    # 0.571 y no 0.625, y en L=12, 14 solo se acota; el 36.9% contaba esos ceros como soporte.
    "36.9":    "$42.0\\%$ at",
    "0.625":   "a floor of $0.571$ at $L=8$ ($2237$ of $3920$)",
    "0.583":   "above $0.54$ at $L=12$ and $14$",
    "8101214": "above $0.54$ at $L=12$ and $14$",
    # El recuento del guardian cambia cada vez que verify.py gana una comprobacion; la prueba
    # es que la frase que lo cita siga en el paper, no un numero concreto (2026-09-25).
    "843":   "numeric assertions, exiting non-zero",
    # 2026-09-25: unprojected power-iteration counts, contaminated by round-off in the
    # reflection-forbidden sector; replaced by the projected counts (z_suppK_sym.py).
    "3919": "that count is not a support", "919": "that count is not a support",
    "0.172": "0.174",   # interval of the overhead exponent, refitted over L=6-12 (2026-09-25)
    "807":  "that count is not a support",
    "831":   "numeric assertions, exiting non-zero",
    "1.010": "1.0100",
    "1.037": "1.0373",
    # 2026-09-28 (R3-m13): '126 of 126' was an earlier state of the same pooled sweep (version
    # history, src/leakage_certificate.py:40-45), deleted from Sec. III D; the fact is the
    # 468 / 363-of-363 sentence that stays.
    "126":   "fails on $469$, and on all $363$ with $w_S\ge0.99$",
    # 2026-09-28 (late): Fig. S5 and the counts quoted from it were REGENERATED from the deposited
    # data/cert_{stress,akw,teqsci}.json (src/recovered/build_n3.py read an earlier run of the same
    # scans until then): 468 -> 469 violations, mildest factor 2.10 -> 2.00, survival below
    # w_S = 0.99 on 138 -> 137 of 243 (verify.py section 9.6c re-classifies the pool and checks
    # every printed count).  Regenerated values, not deleted facts; the proof is the sentence of
    # Sec. III D / Sec. S9 that prints the new value.
    "468":   "fails on $469$, and on all $363$ with $w_S\ge0.99$",
    "2.10":  "the mildest by a factor $2.00$ (Fig.~\\ref{fig:thm1iii-violation}",
    "138":   "it survives on $137$ of $243$",
    # 2026-09-28 (late): the finite-shot penalty of the Fig. 3 configuration was 2.80 on the
    # former window [-9t, 9t] (6.12e-3 / 2.19e-3).  Its generator was recovered and deposited
    # (src/fig3_finite_shot.py, reproducing the 2.80 exactly on that window) and re-run on the
    # window of Fig. 3(a,b): 6.16e-3 / 2.36e-3 = 2.61 (data/fig3_finite_shot_T2.6e6.json;
    # KNOWN_DISCREPANCIES.md section 31).  A regenerated value, not a deleted fact: the proof is
    # the Sec. S1 sentence that prints the new factor (the caption, figs/fig5_caption.tex, is
    # not in this corpus; verify.py section 9.6b checks it against the data).
    "2.80":  "a measured $2.61$ and not the $2.05$",
}

AUDITADA = [
    # (fragmento que dejo de estar,            prueba de que el hecho sigue hoy)
    # 2026-09-25 (recorte para PRA): la introduccion de la Sec. VI se reescribio corta; el
    # hecho de cada frase sigue en la version nueva, cuya frase es la prueba.
    ("None of them can play it, for either term", "None can play that role, for either term of"),
    ("The second is about the leakage",          "the coordinate subspace the sampler returned"),
    ("Only the second is an obstruction on the theorem", "not for the leakage $\\Lambda_S(\\eta)/\\eta$ of the"),
    ("Proposition loses its subject",            "not for the leakage $\\Lambda_S(\\eta)/\\eta$ of the"),
    ("What replaces the predictor is not a cheaper functional", "What replaces the predictor is the measured map of"),
    ("Why it fails: the retained mass is displaced", "not discarded, and the two-level example is the extreme case"),
    ("(Relative here and throughout this subsection", "as the repository's own relative $L_1$ does"),
    ("Three claims of the previous version were larger", "the region in which this method"),
    ("The resource statement of Sec. has a consequence", "we make no advantage claim at any size reached"),
    ("That proof is three lines of Rayleigh--Ritz, and it is not the proof", "and it is not the proof this statement is usually given"),
    # 2026-09-25 (revision final de arbitro): 'published fraction' -> 'operating fraction of the
    # resource scan' donde la afirmacion de vacuidad excluye la Fig. akwsampled (publicada y certificada).
    # 2026-09-30 (reframe skeptics): the temporal hedge 'at present' is removed from Sec. V B, Sec. S3
    # and Sec. VII; the concession itself ('uncertified') is unchanged.  Proof updated to match.
    ("The reconstructions reported in this paper are therefore", "are therefore uncertified"),
    ("Evaluated, the new certificate is vacuous at every published fraction",
     "the certificate is vacuous at every operating fraction of the resource scan"),
    ("The text already said that the molecular set does not separate", "does not separate the two candidate resources"),
    # 2026-09-25 (introduccion reescrita, ~1.2 pp): cada hecho sigue en el texto; la prueba es su frase actual.
    # 2026-09-28 (R3-m13): the priority is dated by the cited preprints (both July 2026), not by
    # this work's first posting, so the proof no longer reads "days before"; the attribution
    # sentence of Sec. I ("the priority for it is theirs") is intact.
    # 2026-09-30 (citation audit, C3-2 / LB-1): Tarabunga et al. prove BOTH of their families
    # Gaussian invariant (arXiv:2607.02242v2, App. B, P.2), so the "invariant/basis-dependent
    # division" was never theirs to be credited with; Sec. I now credits them with the computable
    # invariant quantifiers and the covariance-only / full-state separation, and Leone and Bittel
    # with the strong monotone on the invariant side.  The priority sentence stays, reworded.
    ("Leone and Bittel gave a Gaussian monotone", "both in July 2026; the priority for those measures is theirs"),
    ("The priority for the dichotomy is theirs", "both in July 2026; the priority for those measures is theirs"),
    ("It is false---by an analytic two-level counterexample", "forces the constant of any bound indexed on it to its trivial value"),
    ("It fails wherever it was invoked", "fails on $469$, and on all $363$ with $w_S\\ge0.99$"),
    # 2026-09-28 (late), review of the Fig. 3 regeneration round: four sentences reworded for
    # accuracy, each fact still printed; the proof is its current wording.
    #  - the coefficient-insensitivity statement: the raw support IS lower-semicontinuous (it is
    #    upper semicontinuity that fails), and the CP tensor rank fails BOTH, so "like the CP
    #    tensor rank, not lower-semicontinuous" and "the opposite failure" were both wrong;
    #  - Secs. I and X state only what Theorem thm:witness proves (the O(K) endpoint of the
    #    orbit); "runs from O(K) to exponential" is not established and is no longer claimed;
    #  - the raw-support obstruction is on UPPER bounds ("cannot upper-bound the count more
    #    tightly than the sector dimension"), as the proof shows, not on lower bounds.
    ("like the CP tensor rank, this count is not lower-semicontinuous",
     "unlike the count, is not even lower-semicontinuous"),
    ("Along a geminal orbit every order of it is fixed while the determinant support runs from",
     "the orbit contains a point of classical cost $O(K)$"),
    ("can lower-bound the count at all, rank being coefficient-insensitive",
     "can upper-bound the count more tightly than the sector dimension"),
    # 2026-09-28 (late): the coverage gap verify.py declared for Fig. S5 is closed by its
    # section 9.6c (the plotted ratios are re-derived from the certificate rows); Sec. IX D
    # says so in place of declaring the gap.
    ("one declared gap in its own coverage", "re-derived by it from the deposited certificate rows"),
    # 2026-09-28 (late): the K=18 ranking WAS re-run at L=14 (src/frontier/rerank_L14.py ->
    # data/c3_frontier/published_fraction/L14_rerank.json, set overlap 1.000000 at the published
    # fraction and three others), so the caveat "not re-run at that size" is repaired, not
    # deleted; the proof is the sentence that reports the re-run.  The depth-indicator caveat
    # of the same row stays.
    ("re-run at that size", "reproduces the stored selection as a set to $1.000000$"),
    ("Neither caveat is repaired", "It is not repaired, and the row is reported because it makes the"),
    # 2026-09-28 (R1-01, R3-M2, R3-M3, R3-m10, R1-04): Theorem 1 gains (H4), the exact probe;
    # the claim that only the leakage term is a posteriori is narrowed (both terms need Psi_0);
    # the (H1)-(H3) list printed twice (S2 and S11) is kept once, in S11; the momentum block
    # names the site basis.  Each fact below is still printed; the proof is its current wording.
    ("the weight term is the half this work cannot compute without", "needs the full norm of the probe as well as its retained part"),
    ("needs the exact decomposition of the probe norm", "estimated against it is biased in the optimistic direction"),
    ("Only one of its two terms is a posteriori", "which a user does not have"),
    ("The two conventions are never interchanged silently in this paper", "the two conventions are never interchanged silently"),
    ("If the reference is itself a Lanczos reconstruction of finite depth", "measures the reference rather than the reconstruction"),
    ("the discrepancy acquires an additive", "an additive term verified to"),
    ("the two sizes we measured, the Krylov support fills", "the certificate becomes non-vacuous only close to the full block"),
    ("We retract it explicitly", "with $0$ violations across the $606$ subspaces"),
    ("At every published fraction, at each of the five sizes evaluated, the leakage branch",
     "is worse than the trivial bound at every operating fraction of the resource scan"),
    # 2026-09-30 (reader-value reframe, Sec. I): the defensive hedge "The reconstructions
    # themselves are not in question" is gone; the fact (an independent implementation recovers
    # the acceptance fractions of the scan to under 2%) is printed in Secs. I and V.
    ("The reconstruction itself is untouched", "recovers the acceptance fractions of the scan to under"),
    ("An invariant of the occupation spectrum cannot see a boundary", "An invariant of the occupation spectrum cannot see a boundary"),
    ("The statement is one of non-determination and non-vacuity", "The statement is one of non-determination and non-vacuity, not of non-existence"),
    ("That is not an absence of relation but a deterministic monotone", "a deterministic monotone relation whose sign"),
    ("proves the bound and retracts its predecessor", "proves the bound and shows that no bound on the captured weight alone can replace it"),
    # 2026-09-25 (conclusion reescrita y pasada de tono): el hecho sigue; la prueba es su frase actual.
    ("The two terms of the upper bound are two budgets, and only the first", "the captured weight does not control it"),
    ("The error of a spectral function reconstructed on a subspace of sampled configurations splits",
     "bounds the $L_1$ error from below by the Born weight the subspace failed to capture"),
    ("The bound itself stands where a weight-only bound cannot", "an analytic counterexample and a counting argument rule out any bound indexed on it"),
    ("And the method is priced on three axes rather than one", "The method is priced in support, resolution and shots"),
    ("The three momentum-resolved responses reported here have been computed classically",
     # 2026-09-30 (citation audit, C2-2): Benthien and Jeckelmann resolve N(q,w) in momentum at
     # L=60 (cond-mat/0606748, Fig. 5), reaching L=120 only at the zone boundary; 90--120 was wrong.
     "dynamical DMRG computed the same three channels at $60$ to $120$ sites in 2007"),
    ("Moving from the published fraction to the non-vacuous band", "moving from the operating fraction to the non-vacuous band multiplies the shots per snapshot"),
    ("Calling Eq. an error estimate would be refutable", "is therefore an exclusion, not an error estimate"),
    ("The result is negative and we state it first", "Evaluated on the published subspaces themselves"),
    ("The lattice results are one-dimensional Hubbard rings", "which is the exact-diagonalization convention for this observable"),
    # 2026-09-25: la prueba de descuantizacion sale entera (datos sin generador en src/)
    ("The dequantization test, and a resolution floor", DQ),
    ("Whether that further",                   DQ),
    ("The deposit holds five summary rows",    DQ),
    ("The interval brackets the critical value", DQ),
    ("A second limit is not repaired by repairing the deposit", DQ),
    ("With five paired seeds the exact two-sided sign test", DQ),
    ("The generative comparison goes because", DQ),
    ("The amortized companion result goes",    DQ),
    ("We prove a two-sided bound",             "error lies between $1-w$ and"),
    ("No bound indexed on the captured weight","forces any bound indexed on it to its trivial value"),
    ("ours is worse than the trivial bound",   "the bound is vacuous at every operating fraction of the resource scan"),
    ("What survives is a calibration",         "crossing of the trivial bound"),
    # 2026-09-30 (reframe skeptics): Sec. I now says 'a different lattice' instead of the defensive
    # 'the two lattices are never pooled'; the fact (rings and open chain kept apart) is the same.
    ("These are the periodic rings",           "on the open chain of Sec.~\\ref{sec:moments}, a different lattice"),
    ("The obstruction is structural: the",     "coordinate support on every symmetry-allowed determinant"),
    ("And the certificate is calibrated",      "crossing of the trivial bound"),
    ("The error of a spectral function",       "no inequality indexed on the captured weight can"),
    ("The bound itself stands where",          "no inequality indexed on the captured weight can"),
    ("obstruction is structural rather than",  "coordinate support on every symmetry-allowed determinant"),
    ("Whether a certified regime exists",      "well-posed measurement"),
    ("No such object is built anywhere",       "No such object is built anywhere in this work"),
    ("orders measure the distance between",    "not an error in the published reach"),
    ("The conclusion is right, the route",     "the route is not"),
    ("of the sector buys nothing by itself",   "containment, not spanning"),
    ("Moment exactness is a list of",          "finite linear functional"),
    ("depth is certified per point",           "what the depth cap does to"),
    # 2026-09-28 (R1-01): (H4), the exact probe, is added to the hypotheses of Theorem 1, so
    # the sentence that counted three now counts four; the fact (the hypotheses are stated and
    # checked) is the same sentence with the new count.
    ("Three hypotheses carry the bound",       "Four hypotheses carry the bound"),
    ("the displacement is the looseness",      "the displacement is the looseness"),
    ("The reconstructions themselves are untouched", "with nothing shared"),
    ("The two axes usually priced",            "The third price"),
    ("it counts poles",                        "counts poles"),
    ("It covers the deposited data files",     "not a proof that every sentence"),
    ("It is vacuous at every published",       "vacuous at every published fraction"),

    # 2026-09-25. Estas seis NO eran concesiones: eran afirmaciones FALSAS, refutadas
    # por una revision adversarial independiente y comprobadas a mano, sustituidas
    # por la version estrecha y verdadera. La prueba exige que el sustituto siga ahi.
    #  - el "si y solo si" de Lambda=0 falla con w<1 (H=diag(1,2), phi=e0+e1, S={0});
    #    el criterio correcto es K(H, phi_P);
    #  - "lo que los disparos no compran": C3_grid.json muestra Lambda/eta cayendo de
    #    9.87 a 0.109 (L=8, eta=0.18) a lo largo del orden de Born;
    #  - "cannot be lowered" / "in this paper or another": no demostrado; solo se midio
    #    el margen del escalon de Cauchy-Schwarz.
    ("only the first is bought by drawing more bitstrings", "the captured weight controls only the first"),
    ("is what they do not buy",                "but the captured weight does not control it"),
    ("The separation is the economics",        "What the weight does not control"),
    ("Why the threshold cannot be lowered",    "Where the threshold comes from"),
    ("so the coordinate support of that Krylov space is a floor", "the Krylov space of the retained part of the probe"),
    ("in this paper or another, can be non-vacuous", "every symmetry-allowed determinant"),

    # 2026-09-25. El paper se escribe como un trabajo nuevo: la historia de versiones
    # sale del cuerpo y vive en el campo Comments de arXiv. Cada frase de abajo era
    # historia de versiones; lo que tenia de fisica sigue en el texto (prueba) o la
    # nota de version sigue en Comments (@comments).
    ("the retracted one was never evaluated",  "@comments:Theorem 1(iii) of v1-v2 does not hold"),
    ("An inequality of this form was stated as Theorem 1(iii)", "@comments:Theorem 1(iii) of v1-v2 does not hold"),
    ("strictly the failure is unbounded: in a three-level", "the weight term alone fails without bound: in a three-level family"),
    ("violations across the same",             "subspaces of the pooled sweep"),
    ("none of the three is checked anywhere in the deposit of versions", "We check all four here"),
    ("Where the previous version of this work asserted a regime", "This section reports the measured regime of validity"),
    ("Two things were wrong with that, and we retract both", "for a reason of resolution"),
    ("stood behind the sentence asserting that the two-spinon", "none is used for any channel of this work at this size"),
    ("It is the lower end of the computation window, censored", "none is used for any channel of this work at this size"),
    ("the group that established it belongs in the paragraph", "Vaquero-Sabater, Carreras, Broers"),
    ("The sentence of v1--v2 concluding from this suite", "on twelve of the nineteen there is nothing for efficiency to mean"),

    # 2026-09-25. Dos frases que decian MAS de lo que los datos dan, estrechadas:
    #  - "the reconstructions reported in this paper are uncertified": falso para la
    #    Fig. akwsampled (85 %, 16/16 canales no vacuos en cert_akw.json); la concesion
    #    sigue para las fracciones publicadas;
    #  - "Lambda decreases with depth": no es un teorema y sube hasta un 2 % entre las
    #    profundidades 30 y 100 a L=12; lo que la frase necesitaba si se cumple en las
    #    16 series depositadas (a partir de 100, nunca por debajo del valor mas profundo).
    # 2026-09-30: 'at present' removed (temporal hedge; the concession is unchanged).  The sentence now
    # reads 'The reconstructions at the fractions this paper publishes---...---are therefore uncertified'
    # in Sec. S3 and 'at the operating fractions of the resource scan---...---are therefore uncertified'
    # in Sec. V B; the proof is the clause both share.
    ("The reconstructions reported in this paper are therefore, at present, uncertified", "are therefore uncertified"),
    ("decreases with depth, the bound at any shallower depth", "it is not monotone in depth in general"),
    ("by direct integration at",               "by adaptive integration over the whole line"),
    # 2026-09-25, from the final adversarial read: a false sentence (w_S and K_S ARE computed,
    # cert_*.json carries both, and the same paragraph says so), the shots economics that
    # C3_grid.json contradicts, and a version-history clause.
    ("is computed in any file of the reproduction repository", "of which $357$ have $K_S=-1$"),
    ("separates into two terms with different economics", "the first is fixed by the captured weight and vanishes once"),
    ("The molecular set does not discriminate, and now we can show it", "The molecular set does not discriminate."),
    # 2026-09-25: the structural statement corrected for the reflection symmetry.
    ("certify the reconstructions this paper reports", "does not certify the reconstructions at the operating fractions"),
    ("Second, the threshold cannot be lowered by choosing a better subspace", "not the source of the vacuity measured below"),
    ("The two estimators are independent and they agree", "that count is not a support"),
    ("on the leakage is therefore strictly positive for every proper coordinate subspace", "every symmetry-allowed determinant"),
    # 2026-09-25, overclaims narrowed (final adversarial read, confirmed by its refuter):
    ("What is guaranteed in advance is narrower and survives intact", "the second does not imply the first"),
    ("Its three largest fractions lie above the support floor", "for the local seed does not apply to it"),
    ("whose column MATH is the one whose absence from Table", "-electron dimension; the"),
    # 2026-09-28, audited by hand against the committed records:
    #  - Sec. S8 claimed "recovery to the correct particle-number sector ... which is what a
    #    proof of principle is" for the L=6 ibm_fez run, in the paragraph that discloses its
    #    reversed bit order (data/hw_bitorder_check.json: 0 of 350000 noiseless shots kept).
    #    The sentence was an OVERCLAIM, not a concession; it now claims executability only,
    #    says the post-selection was wrong, and keeps "establishes nothing about accuracy".
    #    verify.py refuses both phrases in sm_scope.tex (hw.sm_scope_*).
    #  - the known-gaps sentence of the data statement listed two gaps that the deposit of
    #    2026-09-26 closed (the six standalone figure sources, the builder of the Fig. S5
    #    tables: docs/KNOWN_DISCREPANCIES.md Sec. 5, 29) and "three" scripts without output
    #    where Sec. 13 now leaves two.  It now names what is still open, which is more, not
    #    less: the three rasters, the hand-completed Fig. 5 and two hand-edited captions.
    ("which is what a proof of principle is",  "It establishes nothing about accuracy, neither at the complete coverage it reached"),
    # 2026-09-28 (R1-03): the coefficient-insensitivity proposition refuted an UPPER bound
    # (Psi_eps: mu -> mu(D0) while the raw support stays K+1) but was stated as a lower-bound
    # failure with an inverted CP-rank analogy.  The raw support is lower- but not
    # upper-semicontinuous; the lower-bound failure holds only for F_1 at half filling (the
    # cat), and F_1 gives the weak bound |S| >= N/(N-F_1/4).  Both statements are now proved as such.
    ("cannot lower-bound the raw determinant support", "the opposite failure to that of the CP tensor rank, which is not lower-semicontinuous"),
    ("can lower-bound the count at all", "not upper-semicontinuous (the opposite failure to that of the CP tensor rank"),
    ("six of the seventeen figures enter as PDF with no standalone", "are deposited as rendered and cannot be regenerated from the deposit"),
    # 2026-09-28 (R2-02): the ibm_fez run is no longer placed in the class of device spectral
    # records: its retained subspace is the whole (N+1) sector, all device errors
    # (data/hw_bitorder_check.json), so its A(w) is the ED result (data/heron_spectral.json).
    # The class comparison and 'not the smallest such record' give way to the stronger
    # statement that it is not in that class; 'not offered as a scale result' stays.
    ("so our device record is smaller than every entry in that class", "the run is not an entry in that class"),
    ("It is not the smallest such record either", "four-qubit trapped-ion Hubbard-dimer reconstruction"),
    # 2026-09-30 (citation audit, 39 findings confirmed by 2-3 of 3 skeptics against the sources):
    # four concession sentences reworded to say what the cited work actually shows.  Each
    # concession stays; the proof is its current wording.
    #  - Sec. S5: sierant2026 defines F_k as a MEASURE of fermionic non-Gaussianity (faithful,
    #    Gaussian invariant, (sub)additive) and leaves monotonicity open; collura2026 never
    #    defines F_k (it computes qubit stabilizer Renyi entropies of Gaussian states);
    #  - Sec. S5: Tarabunga et al.'s natural-orbital participation entropies are FGU-invariant
    #    (App. B), so "basis dependent" was false; the sentence now states their covariance-only
    #    / full-state separation and where the basis dependence of |S| actually enters;
    #  - Sec. S2: chen2021slq's Gauss-quadrature W/KS brackets hold for each seeded (weighted)
    #    measure as well as for the density of states (Cor. 2), so "a global density of states,
    #    not the seeded local measure" misstated the source; the novelty is now fenced on the
    #    full-space quadrature versus the leaking determinant-subspace projection.
    ("It is not a refutation of anyone's claim", "no one proposed it as a predictor of the determinant support"),
    ("It is our own no-go, about our own cost", "a corollary consistent with---not a competitor to---the two families of"),
    # 2026-09-30 (final citation review CC3): the sentence was reworded to say what Tarabunga
    # et al. prove about the occupation-entropy family; the fact is still printed and the proof
    # is its current wording.
    # 2026-09-30 (sources verification, SRC-1): CC3 had taken the arXiv listing abstract, which
    # still carries v1 wording ('one member').  arXiv:2607.02242v2 (28 Jul 2026) proves orders 1
    # AND 2 strong pure-state monotones that lower-bound the SWAP (non-Gaussian gate) count
    # (Sec. II, Sec. IV A, Eq. 48, App. A 2), and its order-2 member is the k=1 antiflatness,
    # i.e. F_1.  Leone-Bittel (arXiv:2607.29670) prove the same quantity M_f = F_1 a strong
    # monotone and their Note added credits v2 with an independent proof of their Theorem 1.
    # Wording corrected to 'two members'; the concession (their priority) is still printed.
    ("set out the same division on general grounds", "occupation entropies depend only on the covariance matrix, two members of which (orders $1$ and $2$) are"),
    ("Jia and Lv's title is this section's claim word for word", "Jia and Lv's title is this section's claim word for word"),
    # 2026-09-30 (reader-value reframe, Secs. VII-IX): seven sentences that stated a result as a
    # concession or as commentary on the writing ('It is easy to price ... It is not.', 'Four
    # statements are supported ... and a fifth is not', 'the strongest physical check available to
    # this method', 'We present it here as what it is') now state the result; the fact of each is
    # still printed, and the proof is its current wording.  Sec. VIII no longer says the gaps come
    # 'from the same primitive': S(q,w) and S^zz are full-sector references (Table I), not outputs
    # of the sampler, and the sentence now says what they check.
    ("It is easy to price the first two", "is a property of the question, not of the model"),
    ("The support exponent, and what does and does not survive", "The support exponent and its resolution dependence"),
    ("Four statements are supported by the measured rows", "do not establish that the shot budget grows with a larger exponent than the"),
    ("Supported: the budget is two orders of magnitude", "two orders of magnitude larger than the support it buys"),
    ("The two collective channels of the half-filled Hubbard chain are the strongest", "against which the sampled reconstructions are scored"),
    ("On the same lattice, in the same calculation, from the same primitive", "reproduce the opposite finite-size scaling of the two collective scales"),
    ("We present it here as what it is", "against which the sampled reconstructions are scored"),
    # 2026-09-30 (reader-value reframe, Sec. I): the clause "we make no priority claim against it"
    # was meta-commentary; the priority fact ("independently and concurrently") stays in the
    # same sentence, which is the proof.
    ("Patel, Rangi and Tam construct Green's functions from sampled Krylov subspaces", "independently and concurrently, and cite this work"),
    # 2026-09-30 (reader-value reframe, Secs. II-IV): seven baseline sentences (headings merged
    # with their first sentence, or meta-commentary such as 'and we never claim that it is',
    # 'we claim no novelty', 'We state that claim as what it is') now state the result; each
    # fact is still printed and the proof is its current wording.
    ("This section states the measurement primitive", "full-sector references computed by continued fraction, not reconstructions"),
    ("what the shots buy and what they do not", "the first is fixed by the captured weight and vanishes once"),
    ("The complete bound is therefore not evaluable", "evaluable by a user who does not already have the ground state"),
    # 2026-09-30 (final math review M3): moment exactness needs the probe (built from the exact
    # ground state), so 'decidable from the subspace alone' now reads 'decidable from the subspace
    # and the probe'; the proof string follows the printed wording (Secs. IV and S10).
    ("What the structure does certify", "decidable from the subspace and the probe"),
    ("It rules out nothing about the construction itself", "decidable from the subspace and the probe"),
    ("We state it with the proof the code realizes", "It is not an error bound, and it does not imply convergence in"),
    ("Three rows is all this validation system has", "These three rows are the whole validation"),
    # 2026-09-30 (reader-value reframe, Secs. V-VI): the merged heading+sentence of Sec. V D
    # ('One positive result survives the sweep') is now the heading 'The crossing of the trivial
    # bound is calibrated'; the merged heading+sentence of Sec. VI B ('The obstruction is not a
    # continuity artifact') keeps its non-transfer sentence intact in the moved paragraph.
    ("One positive result survives the sweep", "The crossing of the trivial bound is calibrated"),
    # 2026-09-30 (reframe skeptics V1, C4).  Two sentences went as meta-commentary or repetition,
    # each fact kept where it is proven:
    #  - Sec. IX C restated the vacuity range, 'at the five sizes measured' and 'uncertified';
    #    the home of all three is Sec. V B (Table published and the sentence that ends
    #    'are therefore uncertified'), and the abstract, Sec. I and Sec. X state the range;
    #  - Sec. S2 said 'and we say so before a referee does'; the fact ('Nothing in the proof
    #    of Appendix A is new as analysis') stays, without the aside.
    ("The certificate of Theorem is vacuous at every fraction this work publishes",
     "are therefore uncertified"),
    ("and we say so before a referee does", "is new as analysis. The split,"),
    ("The obstruction is not a continuity artifact", "does \\emph{not} transfer to the operative"),
    # 2026-09-30 (final reviews M2/CC1 and CC4/CC8).  Three sentences reworded, each fact kept:
    #  - Sec. S5 said F_k does not determine the cost 'being constant on an orbit along which the
    #    cost varies exponentially'; the exponential cost at rotated points is not established
    #    in the deposit (the main text says so), so the sentence now says F_k determines the
    #    cost only if those points are in fact easy.  The hardness half of the concession
    #    ('does not lower-bound the cost') is unchanged in the same sentence
    #    [2026-09-30, final math review MA-1: the cost at the pairing point is O(K) while
    #    F_k = 4K, so F_k/4 does lower-bound a linear cost; the half now reads 'does not
    #    lower-bound the bond dimension, nor the cost beyond linear order in $K$', which is what
    #    the witness refutes.  The proof string below is unaffected];
    #  - Sec. S7 'bound it honestly' lost the adverb (meta-commentary); the two fits stay;
    #  - Sec. S7 pointed at 'the list of what this work does not claim', a list Sec. IX C no
    #    longer has; it now names the restriction Sec. IX C records (pole positions at four
    #    sizes).
    ("Two conclusions follow, each in the only direction the construction supports",
     "determines the cost only if the rotated points are in"),
    ("Two fits make that precise and bound it honestly", "bound it. Fitting $"),
    # 2026-09-30 (final math review MA-2): Proposition 2 proves vacuity on pairing-model fibers,
    # i.e. that no orbital-rotation invariant is a useful predictor on EVERY state, not that none
    # is useful anywhere.  Sec. IX's 'four things' sentence now scopes its obstruction item that
    # way ('an obstruction ... useful on every state'); the other three items are unchanged and
    # the sentence falls to 0.56 similarity only through that clause.  The claim list is still
    # printed; the proof is its current wording.
    ("What the work does claim is four things",
     "obstruction showing that no orbital-rotation invariant supplies an a priori predictor of that bound"),
    ("Section states the same concession in the list of what this work does not claim",
     "restriction: the spin--charge statement rests on pole positions at four sizes."),
]


def _prueba_ok(prueba, hoy_crudo):
    """Una prueba es un texto que tiene que seguir en el paper de hoy, o -- si empieza
    por '@comments:' -- en build/arxiv_comments.txt: la historia de versiones vive en
    el campo Comments de arXiv, no en el cuerpo (decision del autor, 2026-09-25)."""
    if prueba.startswith("@comments:"):
        c = os.path.join(REPO, "build", "arxiv_comments.txt")
        return os.path.exists(c) and prueba[len("@comments:"):] in io.open(c, encoding="utf-8").read()
    if prueba.startswith("@doc:"):
        # '@doc:<ruta>:<texto>' -- una retirada DELIBERADA de un resultado entero, con su
        # razon escrita en la documentacion del deposito (2026-09-25: la prueba de
        # descuantizacion, cuyo fichero de datos no tiene generador en src/).
        ruta, texto = prueba[len("@doc:"):].split(":", 1)
        c = os.path.join(REPO, ruta)
        return os.path.exists(c) and texto in " ".join(io.open(c, encoding="utf-8").read().split())
    return prueba in hoy_crudo


def _auditada(frase, hoy_crudo):
    """Perdona una frase SOLO si su prueba sigue en el paper de hoy."""
    for clave, prueba in AUDITADA:
        if clave in frase:
            return _prueba_ok(prueba, hoy_crudo)
    return False


def revisa_numeros(antes, ahora):
    na, nh = numeros(antes), numeros(ahora)
    muertos_d, muertos_c, caidos = [], [], []
    for n, c in sorted(na.items(), key=lambda kv: -kv[1]):
        c2 = nh.get(n, 0)
        if c2 == 0:
            if n in AUDITADO_NUM and _prueba_ok(AUDITADO_NUM[n], " ".join(ahora.split())):
                continue                   # auditado: ver AUDITADO_NUM
            (muertos_d if distintivo(n) else muertos_c).append((n, c))
        elif c2 < c:
            caidos.append((n, c, c2))
    return muertos_d, muertos_c, caidos


def revisa_concesiones(antes, ahora, umbral=0.60):
    # El corpus conserva los saltos de linea del fuente, y una prueba puede cruzar
    # uno; se busca siempre sobre el texto con los espacios ya normalizados.
    plano = " ".join(ahora.split())
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
        elif _auditada(f, plano):
            continue                       # auditada a mano, y su prueba sigue
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
