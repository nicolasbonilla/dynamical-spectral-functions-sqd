# -*- coding: utf-8 -*-
"""Regenerate every data-driven figure fragment into a scratch tree and diff it against
the committed one.

Why this exists.  `make figures` overwrites files in paper/figs/ that are committed and
that the manuscript compiles.  Until 2026-09-18 the only guarantee that those generators
still reproduced the committed fragments was a sentence in the Makefile -- and a sentence
is not a test.  An adversarial pass demonstrated the consequence: reintroducing the
`%`-format bug in make_akw_sampled_honest_fig.py, breaking make_table.py's output directory, and
corrupting a number inside a generated fragment all left `python src/verify.py` reporting
PASS.  Those three repairs were unguarded; this file guards them.

It never writes into the repository.  Everything happens in a temporary directory that is
removed on exit, so running it can neither fix nor break paper/figs/.

    python src/check_figures.py          # exit 0 iff every fragment regenerates identically
    python src/check_figures.py -v       # also list the files that matched

Deliberately NOT run here: src/make_scaling_fig.py.  It emits a different, two-panel
figure and would destroy panel (c) of the committed fig_scaling2_native.tex.  See
docs/KNOWN_DISCREPANCIES.md.

Added 2026-09-26: the two generators recovered into src/recovered/ (see its README.md) --
make_fig_gapscaling.py (Fig. 6) and build_n3.py (the ten n3_*.dat tables of Fig. S5).  Their
eleven artefacts are NOT seeded into the scratch tree before the run, so a generator that
fails to write one is reported MISSING instead of passing on the committed copy.  One
comment line is normalised before comparing, and the run prints it when it does: line 5 of
fig_gapscaling_native.tex records the generated_utc of the gap_scaling.json the fragment was
built from, and the deposited JSON is a later run with identical numbers (NORMALISE below).

Added 2026-09-28 (late): src/recovered/n3_fig5_whiskers.py, run after build_n3.py.  It rewrites
the two whisker blocks of paper/figs/fig_thm1iii_violation_native.tex from the n3_fig5_*.dat
tables and writes data/thm1iii_violation/n3_fig5_summary.json; both are compared, the fragment
byte-for-byte (it is seeded, because the script rewrites it in place and build_n3.py reads its
legend) and the JSON after being removed from the scratch copy of data/.  Nine generators,
twenty-one artefacts.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

GENERATORS = [
    "make_decoupling_native.py",
    "make_method_fig_max.py",
    "make_akw_sampled_honest_fig.py",
    "make_noise_recovery_native.py",
    "make_hardware_hero.py",
    "make_table.py",
    # recovered 2026-09-26 (src/recovered/README.md)
    "recovered/make_fig_gapscaling.py",
    "recovered/build_n3.py",
    # added 2026-09-28 (late): the min/median/max whiskers of Fig. S5, rewritten into the
    # fragment from the two n3_fig5_*.dat build_n3.py has just written (it replaces the
    # one-shot patcher fix_n3_fig5_series.py, which is kept as the record).  Runs after
    # build_n3.py, in this order.
    "recovered/n3_fig5_whiskers.py",
]

# What those generators are expected to (re)write, relative to the repository root (the
# first six generators write the first eight artefacts; the two recovered ones the rest).
# make_table.py writes into rebuild/, a scratch directory, and is therefore compared
# against paper/table_molecules.tex only for its DATA ROWS -- the committed table carries
# five lines of hand-written caption the generator does not produce (see
# docs/KNOWN_DISCREPANCIES.md); that comparison is not attempted here.
ARTEFACTS = [
    "paper/figs/fig_decoupling_native.tex",
    # Added 2026-09-19.  make_decoupling_native.py now emits BOTH fragments of the resource
    # section (it is the deposited port of the generator that rebuilt them for v3), so the
    # second one is checked too.  This widens the guarantee; it does not soften it.
    "paper/figs/fig_resource_master_native.tex",
    "paper/figs/fig_method_native_frag.tex",
    "paper/figs/aw_method.dat",
    "paper/figs/fig_akw_sampled_honest.tex",
    "paper/figs/fig_noise_recovery_native.tex",
    "paper/figs/fig_hardware_hero_frag.tex",
    "paper/figs/heron_hot.dat",
    # Added 2026-09-26: written by the two recovered generators.
    "paper/figs/fig_gapscaling_native.tex",
    "paper/figs/n3_pub_frac_hi.dat",
    "paper/figs/n3_pub_frac_lo.dat",
    "paper/figs/n3_pub_eta_hi.dat",
    "paper/figs/n3_pub_eta_lo.dat",
    "paper/figs/n3_thmB_frac.dat",
    "paper/figs/n3_thmB_eta.dat",
    "paper/figs/n3_inf_frac.dat",
    "paper/figs/n3_inf_eta.dat",
    "paper/figs/n3_fig5_frac.dat",
    "paper/figs/n3_fig5_eta.dat",
    # Added 2026-09-28 (late): written by n3_fig5_whiskers.py.  The fragment is seeded (the
    # script rewrites its two whisker blocks in place and build_n3.py reads its legend), so
    # it is compared byte-for-byte after the rewrite; the summary JSON is removed from the
    # scratch copy of data/ before the run, so the script must write it.
    "paper/figs/fig_thm1iii_violation_native.tex",
    "data/thm1iii_violation/n3_fig5_summary.json",
]

# Artefacts that are compared but NOT copied into the scratch tree first: the generator must
# write them, or they are MISSING.  (The eight older ones are seeded, as they always were.)
NOT_SEEDED = {a for a in ARTEFACTS
              if a.endswith("/fig_gapscaling_native.tex") or a.startswith("paper/figs/n3_")
              or a.endswith("/n3_fig5_summary.json")}

# Read by a generator, never written or compared.  Until 2026-09-28 (late) this held the
# Fig. S5 fragment, whose legend build_n3.py checks; n3_fig5_whiskers.py now rewrites that
# fragment, so it is an ARTEFACT above (seeded, then compared).
INPUTS = []

# Lines normalised on BOTH sides before comparing, per artefact.  Keep this list short and
# every entry a comment line: it is an exemption from a byte-for-byte guard.
NORMALISE = {
    "paper/figs/fig_gapscaling_native.tex": [
        (re.compile(br"^(%   gap_scaling\.json   \(gap_scaling\.py, )"
                    br"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ(\))$", re.M),
         br"\g<1><generated_utc>\g<2>"),
    ],
}


def _normalise(rel, data):
    data = data.replace(b"\r\n", b"\n")
    n = 0
    for rx, rep in NORMALISE.get(rel, ()):
        data, k = rx.subn(rep, data)
        n += k
    return data, n


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    verbose = "-v" in argv or "--verbose" in argv

    tmp = tempfile.mkdtemp(prefix="checkfigs_")
    try:
        shutil.copytree(os.path.join(ROOT, "src"), os.path.join(tmp, "src"))
        shutil.copytree(os.path.join(ROOT, "data"), os.path.join(tmp, "data"))
        os.makedirs(os.path.join(tmp, "paper", "figs"))
        for rel in [a for a in ARTEFACTS if a not in NOT_SEEDED] + INPUTS:
            src = os.path.join(ROOT, rel.replace("/", os.sep))
            if os.path.exists(src):
                shutil.copy2(src, os.path.join(tmp, rel.replace("/", os.sep)))
        for rel in NOT_SEEDED:              # data/ is copied whole: unseed what lives there
            p = os.path.join(tmp, rel.replace("/", os.sep))
            if os.path.exists(p):
                os.remove(p)

        print("check-figures: regenerating %d generators in a scratch tree" % len(GENERATORS))
        failed = []
        for g in GENERATORS:
            p = subprocess.run([sys.executable, os.path.join("src", g)], cwd=tmp,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            if p.returncode != 0:
                failed.append((g, p.returncode, p.stdout.decode("utf-8", "replace")[-1500:]))
        if failed:
            for g, rc, out in failed:
                print("\n  [CRASH] src/%s exited %d\n%s" % (g, rc, out))
            print("\nRESULT: FAIL -- %d figure generator(s) do not run. `make figures` is "
                  "broken." % len(failed))
            return 1

        differ, missing, normalised = [], [], []
        for rel in ARTEFACTS:
            a = os.path.join(ROOT, rel.replace("/", os.sep))
            b = os.path.join(tmp, rel.replace("/", os.sep))
            if not os.path.exists(b):
                missing.append(rel)
                continue
            fa = open(a, "rb").read() if os.path.exists(a) else None
            fb = open(b, "rb").read()
            if fa is None:
                missing.append(rel)
                continue
            na_, ka = _normalise(rel, fa)
            nb_, kb = _normalise(rel, fb)
            if na_ != nb_:
                differ.append(rel)
            elif ka or kb:
                normalised.append(rel)
                print("  [same*] %s  (identical after normalising %d provenance-timestamp "
                      "comment line(s); see NORMALISE)" % (rel, max(ka, kb)))
            elif verbose:
                print("  [same ] %s" % rel)

        if missing:
            print("\n  [MISSING] the generators did not produce: %s" % missing)
        if differ:
            for rel in differ:
                print("\n  [DIFFER] %s" % rel)
                a = open(os.path.join(ROOT, rel.replace("/", os.sep)),
                         encoding="utf-8", errors="replace").read().split("\n")
                b = open(os.path.join(tmp, rel.replace("/", os.sep)),
                         encoding="utf-8", errors="replace").read().split("\n")
                shown = 0
                for i in range(max(len(a), len(b))):
                    la = a[i] if i < len(a) else "<absent>"
                    lb = b[i] if i < len(b) else "<absent>"
                    if la != lb:
                        print("      line %d\n        committed:   %s\n        regenerated: %s"
                              % (i + 1, la[:160], lb[:160]))
                        shown += 1
                        if shown >= 3:
                            print("      ... (%d differing lines in total)"
                                  % sum(1 for j in range(max(len(a), len(b)))
                                        if (a[j] if j < len(a) else None)
                                        != (b[j] if j < len(b) else None)))
                            break
        if differ or missing:
            print("\nRESULT: FAIL -- the committed figure sources are NOT what the deposited "
                  "generators produce from the deposited data. Either a generator drifted or "
                  "a fragment was hand-edited; in both cases the paper and the deposit "
                  "disagree.")
            return 1

        print("RESULT: PASS -- all %d regenerated artefacts are byte-identical to the "
              "committed ones%s." % (len(ARTEFACTS),
                                     (" (%d of them after normalising one provenance-timestamp "
                                      "comment line)" % len(normalised)) if normalised else ""))
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
