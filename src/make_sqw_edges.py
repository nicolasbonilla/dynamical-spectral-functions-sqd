# -*- coding: utf-8 -*-
r"""Generator of paper/figs/sqw_edges.dat -- band edges and peak position of the charge
structure factor S(q,omega) of the 1D Hubbard ring (L=12, U/t=8, eta=0.18).

WHY THIS FILE EXISTS
--------------------
Until now paper/figs/sqw_edges.dat had NO generator anywhere in the repository (docs/REPRODUCE.md
line 38 claims sqw_lanczos.py writes it; it does not -- only spinqw_edges.dat is written, by
spin_lanczos.py:42).  Two of its thirteen rows, q/pi = 0.0000 and q/pi = 2.0000, could not be
produced by any computation in this repository and were not measurements of S(q,omega):

    0.0000 4.969 5.000 5.000
    2.0000 4.969 5.000 5.000

  * q = 0 is INACCESSIBLE by construction: the Lanczos driver skips it (sqw_lanczos.py:30,
    "skip q=0") because n_{q=0} is the conserved particle number, so S(q=0,omega>0) = 0
    identically (sqw_field.py:3; sqw_field.py:11-12 pads q=0 with zeros for the raster).
  * 4.969 is the CHARGE GAP Delta = mu^+ - mu^- of the ONE-PARTICLE channel, snapped to the
    nearest S(q,omega) grid point (see src/charge_gap_ed.py: Delta = 4.968759; grid point
    4.969115).  It is not a measurement of S(q,omega).
  * 5.000 is not a point of the omega grid at all (neighbours 4.96911 and 4.99248).

This script reproduces the eleven legitimate rows BYTE FOR BYTE from data/sqw_L12.json using the
documented recipe, and cannot produce the two fabricated ones.

PAPER CLAIM THIS BACKS
----------------------
The lower-edge / upper-edge / peak markers overlaid on the S(q,omega) figure (Fig. 6 of v2).

It also removes the support for the claim of results.tex:99-100 that "as q -> 0 its lower edge
approaches the charge gap Delta".  What the eleven computable rows actually show, at the 3 per
cent threshold this recipe uses (Delta = 4.968759, see src/charge_gap_ed.py):

    q/pi    1/6     1/3     1/2     2/3     5/6      1
    onset   5.600   5.647   6.418   7.003   6.862   7.119
    onset/Delta  1.127   1.137   1.292   1.409   1.381   1.433

TWO HONEST CAVEATS, because an earlier version of this note got both wrong.

(1) The onset is NOT monotone in q.  Reading the table from q/pi = 1 downwards it goes
    7.119 -> 6.862 -> 7.003 -> 6.418 -> 5.647 -> 5.600: it RISES between q/pi = 5/6 and 2/3.
    The defensible statement is that the onset is lowest at the two smallest accessible
    momenta (5.600, 5.647) and flat between them -- not that it "falls as q decreases",
    which the table three lines above contradicts.

(2) The excess over Delta is a function of the threshold and must never be quoted without it.
    Measured on this same dataset by varying only THRESH_FRAC:

        threshold   0.5%    1.0%    2.0%    2.5%    3.0%    5.0%
        min onset   5.063   5.366   5.530   5.577   5.600   5.717
        min excess  +1.9%   +8.0%  +11.3%  +12.2%  +12.7%  +15.1%
        max excess  +8.9%  +21.6%  +34.8%  +39.0%  +43.3%  +57.4%

    At 0.5 per cent the lowest onset is 1.019 Delta and the data essentially DO reach the gap;
    at 3 per cent they do not come within 12.7 per cent of it.  So "+12.7% to +43.3%" is a
    statement about a convention, exactly as the "10-39%" it supersedes was, and this file will
    not sell either as a fact.  What survives every threshold is the point that matters: the
    figure's approach to Delta was carried entirely by the hand-inserted q = 0 row, whose lower
    edge was SET equal to Delta -- the illustration assumed what it illustrated -- and no
    computable row at any threshold sits AT Delta.

RECIPE (identical to spin_lanczos.py:35-41, the only documented edge recipe in the repository)
---------------------------------------------------------------------------------------------
    thr  = 0.03 * max over ALL q of S(q,omega)          (3% of the global maximum)
    lo   = min{ omega : S(q,omega) >  thr }             lower edge  ("onset")
    hi   = max{ omega : S(q,omega) >  thr }             upper edge  ("upper")
    peak = argmax_omega S(q,omega)                      ("peak")
    all three printed with "%.3f"; q printed with "%.4f"

Deterministic: no random numbers are used anywhere (random_seeds_used = null in the provenance).

USAGE
-----
    python src/make_sqw_edges.py            # writes paper/figs/sqw_edges.dat (+ provenance JSON)
    python src/make_sqw_edges.py --check    # verify only, write nothing

Self-verification (aborts with a non-zero exit status on any failure):
  (V1) the input carries no q = 0 and no q = 2 column, and carries EXACTLY the L-1 = 11
       accessible lattice momenta q/pi = 2n/L -- so a DELETED momentum fails too, which
       an earlier version of this check did not notice;
  (V2) every q has signal above the threshold, and onset <= peak <= upper.  (The earlier
       (V2), "every emitted value is a point of the omega grid", was unfalsifiable: those
       values are selected FROM the grid, so it could only ever fire through its own
       empty-signal branch.  It was the 4*K == 4.0*K pattern again.)
  (V3) no row has upper == peak.  This IS the fingerprint of the two fabricated rows, but
       only at the threshold used here: sweeping THRESH_FRAC upwards, this same genuine
       dataset first produces upper == peak at 0.34 (q/pi = 1/2).  The earlier claim that
       "no genuine spectrum here produces" it was false;
  (V4) the eleven rows reproduce the known regression values;
  (V5) --check COMPARES AGAINST THE DEPOSITED FILE.  The earlier --check regenerated the
       rows and compared them with the list in this file, never opening
       paper/figs/sqw_edges.dat: it printed ALL PASS with the two fabricated rows sitting
       in the artefact it claimed to be checking, and with a legitimate row deleted.
"""
import json, os, sys, hashlib, datetime
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC_JSON = os.path.join(REPO, 'data', 'sqw_L12.json')
OUT_DAT = os.path.join(REPO, 'paper', 'figs', 'sqw_edges.dat')
OUT_PROV = os.path.join(REPO, 'data', 'sqw_edges_provenance.json')

HEADER = "qpi onset upper peak"
THRESH_FRAC = 0.03

# (V4) regression: the eleven rows the repository has always had, and the only ones this
# recipe can produce.  Checked verbatim before anything is written.
REGRESSION = [
    "0.1667 5.600 7.704 5.740", "0.3333 5.647 10.415 6.816",
    "0.5000 6.418 10.836 7.750", "0.6667 7.003 10.906 8.405",
    "0.8333 6.862 10.649 8.732", "1.0000 7.119 10.298 8.779",
    "1.1667 6.862 10.649 8.732", "1.3333 7.003 10.906 8.405",
    "1.5000 6.418 10.836 7.750", "1.6667 5.647 10.415 6.816",
    "1.8333 5.600 7.704 5.740",
]


def fail(msg):
    sys.stderr.write("FAIL: " + msg + "\n")
    sys.exit(1)


def build():
    d = json.load(open(SRC_JSON))
    wg = np.array(d['wg'], dtype=float)
    qs = sorted(float(k) for k in d['S'].keys())

    # (V1) q = 0 and its mirror q = 2 must be absent: inaccessible by construction.
    for forbidden in (0.0, 2.0):
        if any(abs(q - forbidden) < 1e-9 for q in qs):
            fail("input contains q/pi = %.1f, which S(q,omega) cannot resolve "
                 "(density conservation; sqw_lanczos.py:30)" % forbidden)
    # (V1b) ... and the momenta must be EXACTLY the L-1 accessible lattice momenta.
    # Insertion was watched before; deletion was not.
    legal = [2.0 * n / d['L'] for n in range(1, d['L'])]
    if len(qs) != len(legal) or any(min(abs(q - g) for g in legal) > 1e-4 for q in qs):
        fail("input momenta are not the %d accessible lattice momenta q/pi = 2n/L of an "
             "L=%d ring: got %s" % (len(legal), d['L'], [round(q, 5) for q in qs]))

    M = np.array([d['S']["%.5f" % q] for q in qs], dtype=float)
    thr = THRESH_FRAC * M.max()

    rows, diag = [], []
    for q in qs:
        A = np.array(d['S']["%.5f" % q], dtype=float)
        sig = wg[A > thr]
        lo = float(sig.min()) if len(sig) else 0.0
        hi = float(sig.max()) if len(sig) else 0.0
        pk = float(wg[int(np.argmax(A))])
        # (V2) the signal must be non-empty and the three markers must be ordered.
        # onset <= peak <= upper is a real constraint on the emitted row -- it is what the
        # fabricated rows (4.969, 5.000, 5.000) violated by putting the peak ON the upper
        # edge -- whereas "every value is a grid point" was true by construction.
        if not len(sig):
            fail("q/pi=%.4f: no S(q,omega) above the %.1f%% threshold; the row would be "
                 "(0,0,0)" % (q, 100 * THRESH_FRAC))
        if not (lo <= pk <= hi):
            fail("q/pi=%.4f: markers out of order, onset=%.3f peak=%.3f upper=%.3f"
                 % (q, lo, pk, hi))
        # (V3) upper == peak: the fingerprint of the two fabricated rows AT THIS THRESHOLD
        # (this dataset does produce it above THRESH_FRAC = 0.34).
        if abs(hi - pk) < 1e-9:
            fail("q/pi=%.4f: upper == peak at the %.1f%% threshold -- the fingerprint of the "
                 "two rows removed on 2026-09-18" % (q, 100 * THRESH_FRAC))
        rows.append("%.4f %.3f %.3f %.3f" % (q, lo, hi, pk))
        diag.append(dict(qpi=q, onset=lo, upper=hi, peak=pk))

    # (V4) regression against the eleven long-standing rows
    if rows != REGRESSION:
        for a, b in zip(rows, REGRESSION):
            if a != b:
                sys.stderr.write("  regenerated: %r\n  on record  : %r\n" % (a, b))
        fail("regenerated rows differ from the eleven legitimate rows on record")

    meta = dict(L=d['L'], U=d['U'], eta=d['eta'], E0=d['E0'],
                nw=len(wg), omega_min=float(wg[0]), omega_max=float(wg[-1]),
                global_max=float(M.max()), threshold=float(thr),
                threshold_fraction=THRESH_FRAC, n_q=len(qs))
    return rows, diag, meta


def main():
    check_only = "--check" in sys.argv
    rows, diag, meta = build()
    text = "\n".join([HEADER] + rows) + "\n"

    print("recipe     : thr = %.2f%% of global max = %.10f" % (100 * THRESH_FRAC, meta['threshold']))
    print("input      : data/sqw_L12.json  (L=%d, U=%.1f, eta=%.2f, %d q-points, %d omega points)"
          % (meta['L'], meta['U'], meta['eta'], meta['n_q'], meta['nw']))
    sys.stdout.write(text)
    print("checks     : V1 momenta are exactly the 11 accessible ones | V2 signal non-empty "
          "and onset<=peak<=upper | V3 no upper==peak | V4 regression -> PASS")

    if check_only:
        # (V5) THE POINT OF --check: compare against the file on disk.  Regenerating the
        # rows and comparing them with the list in this file checks nothing about the
        # deposited artefact.
        if not os.path.exists(OUT_DAT):
            fail("--check: paper/figs/sqw_edges.dat does not exist")
        on_disk = open(OUT_DAT, "rb").read().decode("utf-8").replace(chr(13) + chr(10), chr(10))
        if on_disk != text:
            a = on_disk.split(chr(10))
            b = text.split(chr(10))
            for i in range(max(len(a), len(b))):
                la = a[i] if i < len(a) else "<absent>"
                lb = b[i] if i < len(b) else "<absent>"
                if la != lb:
                    sys.stderr.write("  line %d" % (i + 1) + chr(10)
                                     + "    deposited  : %r" % la + chr(10)
                                     + "    regenerated: %r" % lb + chr(10))
            fail("--check: paper/figs/sqw_edges.dat is NOT what data/sqw_L12.json produces. "
                 "A row has been added, deleted or altered.")
        print("V5 deposited paper/figs/sqw_edges.dat matches byte for byte -> ALL PASS")
        print("--check: nothing written.")
        return

    old = open(OUT_DAT, 'rb').read() if os.path.exists(OUT_DAT) else b''
    f = open(OUT_DAT, 'w', newline='\n')
    f.write(text)
    f.close()
    new = open(OUT_DAT, 'rb').read()
    prov = dict(
        generated_by="src/make_sqw_edges.py",
        generated_utc=datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        input_file="data/sqw_L12.json",
        input_sha256=hashlib.sha256(open(SRC_JSON, 'rb').read()).hexdigest(),
        output_file="paper/figs/sqw_edges.dat",
        output_sha256=hashlib.sha256(new).hexdigest(),
        previous_output_sha256=hashlib.sha256(old).hexdigest() if old else None,
        recipe="thr=0.03*max_{q,w} S(q,w); onset=min{w: S>thr}; upper=max{w: S>thr}; peak=argmax_w S",
        recipe_source="src/spin_lanczos.py:35-41",
        deterministic=True, random_seeds_used=None,
        rows=diag, meta=meta,
        note=("q/pi = 0 and q/pi = 2 are inaccessible by construction (density conservation; "
              "sqw_lanczos.py:30, sqw_field.py:3) and therefore cannot appear. The two rows "
              "'0.0000 4.969 5.000 5.000' and '2.0000 4.969 5.000 5.000' present in earlier "
              "versions of sqw_edges.dat had no computational origin: 4.969 is the one-particle "
              "charge gap Delta snapped to the S(q,omega) grid, and 5.000 is not on the grid."),
    )
    json.dump(prov, open(OUT_PROV, 'w'), indent=1)
    print("WROTE paper/figs/sqw_edges.dat (%d bytes, sha256 %s...)"
          % (len(new), prov['output_sha256'][:16]))
    print("WROTE data/sqw_edges_provenance.json")


if __name__ == '__main__':
    main()
