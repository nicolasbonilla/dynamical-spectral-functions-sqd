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
    ("The reconstructions reported in this paper are therefore", "are therefore, at present, uncertified"),
    ("Evaluated, the new certificate is vacuous at every published fraction",
     "the certificate is vacuous at every operating fraction of the resource scan"),
    ("The text already said that the molecular set does not separate", "does not separate the two candidate resources"),
    # 2026-09-25 (introduccion reescrita, ~1.2 pp): cada hecho sigue en el texto; la prueba es su frase actual.
    ("Leone and Bittel gave a Gaussian monotone", "days before; the priority for it is theirs"),
    ("The priority for the dichotomy is theirs", "days before; the priority for it is theirs"),
    ("It is false---by an analytic two-level counterexample", "forces the constant of any bound indexed on it to its trivial value"),
    ("It fails wherever it was invoked", "first sweep it failed on $126$ of $126$ such subspaces"),
    ("We retract it explicitly", "with $0$ violations across the $606$ subspaces"),
    ("At every published fraction, at each of the five sizes evaluated, the leakage branch",
     "is worse than the trivial bound at every operating fraction of the resource scan"),
    ("The reconstruction itself is untouched", "The reconstructions themselves are not in question"),
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
     "dynamical DMRG computed the same three channels at $90$ to $120$ sites in 2007"),
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
    ("These are the periodic rings",           "the two lattices are never pooled"),
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
    ("Three hypotheses carry the bound",       "Three hypotheses carry the bound"),
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
    ("none of the three is checked anywhere in the deposit of versions", "We check all three here"),
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
    ("The reconstructions reported in this paper are therefore, at present, uncertified", "The reconstructions at the published fractions are therefore, at present, uncertified"),
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
