# -*- coding: utf-8 -*-
"""Compare paper/arxiv-submission.tar.gz with the current paper/ tree, and say so.

The bundle is a FROZEN snapshot of what was posted as arXiv v2.  Every edit to paper/
moves the tree ahead of it, so any sentence of the form "N of its 47 files are identical"
goes stale the moment anybody touches a figure -- which is exactly what happened on
2026-09-18: the pass whose job was to remove stale derived numbers wrote "eight differ"
after measuring the bundle BEFORE its own two edits, and left a false count standing in
three documents.

So the count is not written down anywhere by hand any more.  This script measures it, and
src/verify.py checks the sentences in README.md, docs/ARXIV_SUBMISSION.md and
docs/KNOWN_DISCREPANCIES.md against what this script measures.

    python src/check_tarball.py

Exit status is 0 whatever it finds: a bundle that is behind the tree is expected, not a
defect.  The defect would be a document that states the wrong number, and that is
verify.py's job.
"""
import os
import sys
import tarfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BUNDLE = os.path.join(ROOT, "paper", "arxiv-submission.tar.gz")


def compare(root=ROOT):
    """Returns (identical, differing, only_in_bundle) as sorted lists of paths
    relative to paper/."""
    same, diff, only = [], [], []
    with tarfile.open(os.path.join(root, "paper", "arxiv-submission.tar.gz")) as tf:
        for m in tf.getmembers():
            if not m.isfile():
                continue
            rel = m.name[2:] if m.name.startswith("./") else m.name
            local = os.path.join(root, "paper", rel.replace("/", os.sep))
            data = tf.extractfile(m).read()
            if not os.path.exists(local):
                only.append(rel)
                continue
            here = open(local, "rb").read()
            if here.replace(b"\r\n", b"\n") == data.replace(b"\r\n", b"\n"):
                same.append(rel)
            else:
                diff.append(rel)
    return sorted(same), sorted(diff), sorted(only)


def main():
    same, diff, only = compare()
    total = len(same) + len(diff) + len(only)
    print("paper/arxiv-submission.tar.gz vs the current paper/ tree")
    print("  files in the bundle : %d" % total)
    print("  byte-identical      : %d" % len(same))
    print("  differ              : %d" % len(diff))
    for f in diff:
        print("      %s" % f)
    print("  present only in the bundle (deleted from the tree since v2) : %d" % len(only))
    for f in only:
        print("      %s" % f)
    print("\nThe bundle is the frozen arXiv v2 source.  It must be rebuilt before any v3")
    print("submission; until then it would re-upload the superseded files listed above.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
