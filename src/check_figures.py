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
"""
import os
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
]

# What those six generators are expected to (re)write, relative to the repository root.
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
]


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    verbose = "-v" in argv or "--verbose" in argv

    tmp = tempfile.mkdtemp(prefix="checkfigs_")
    try:
        shutil.copytree(os.path.join(ROOT, "src"), os.path.join(tmp, "src"))
        shutil.copytree(os.path.join(ROOT, "data"), os.path.join(tmp, "data"))
        os.makedirs(os.path.join(tmp, "paper", "figs"))
        for rel in ARTEFACTS:
            src = os.path.join(ROOT, rel.replace("/", os.sep))
            if os.path.exists(src):
                shutil.copy2(src, os.path.join(tmp, rel.replace("/", os.sep)))

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

        differ, missing = [], []
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
            elif fa.replace(b"\r\n", b"\n") != fb.replace(b"\r\n", b"\n"):
                differ.append(rel)
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
              "committed ones." % len(ARTEFACTS))
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
