# -*- coding: utf-8 -*-
r"""GATE CONDITION 2, checked against the repository's own artefact.

The repo ships work/ckpt/L10_order.npy, produced by release/src/scaling_lanczos_mf.py
stage_rank with K=18, dt=0.5, sub=16 (19 snapshots).  This script reproduces that
ranking with fcore.ranking() and compares, so that the claim "the sweep uses the repo
protocol" is a measurement, not an assertion.  READ-ONLY on the repo.
"""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
# --- deposit paths (repaired 2026-09-19; see docs/C3_FRONTIER.md) ----------
SRC = os.path.dirname(HERE)                       # release/src
ROOT = os.path.dirname(SRC)                       # repository root
OUT = os.environ.get('C3_OUT') or os.path.join(ROOT, 'data', 'c3_frontier')
TMP = os.environ.get('C3_TMP') or os.path.join(ROOT, 'build')
# --------------------------------------------------------------------------
sys.dont_write_bytecode = True
sys.path.insert(0, HERE)
import fcore as F
from fcore import log

CK = os.environ.get('P2_CKPT') or os.path.join(ROOT, 'ckpt')

L = 10
# DECLARED GAP.  work/ckpt/L10_order.npy is a 0.2 MB by-product of
# src/scaling_lanczos_mf.py and is NOT in the deposit (ckpt/ is .gitignore'd).
# Point P2_CKPT at a tree that has it, or regenerate it with
# `python src/scaling_lanczos_mf.py 10` first.  Exit 3 = gap, not failure.
_ref_path = os.path.join(CK, f"L{L}_order.npy")
if not os.path.isfile(_ref_path):
    print("SKIP: %s is not in the deposit; set P2_CKPT to a tree that has it "
          "(see docs/C3_FRONTIER.md, 'what is not deposited')." % _ref_path)
    sys.exit(3)
ref = np.load(_ref_path)
sysd = F.build_system(L)
D = sysd['D']
log(f"repo order file: {ref.shape[0]} entries; my sector D={D}")

for K in (18, 16):
    order, wc, t = F.ranking(sysd, K=K)
    agree = {}
    for FR in (0.18, 0.36, 0.50, 0.70, 0.85, 0.90):
        k = max(1, int(round(FR * D)))
        a = set(order[:k].tolist()); b = set(ref[:k].tolist())
        agree[FR] = len(a & b) / k
    log(f"K={K} [{t:.1f}s]: rank-1 match={order[0]==ref[0]}  "
        + "  ".join("FR=%.2f set-overlap=%.6f" % (f, v) for f, v in agree.items()))
