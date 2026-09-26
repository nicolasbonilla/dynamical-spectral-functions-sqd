# `_superseded/v1_readme_img/` — the README thumbnails of v1

README thumbnails cropped from the v1 PDF on 2026-08-17; superseded 2026-09-26: scaling.png showed withdrawn v1 text, sqw.png the two fabricated sqw_edges.dat rows, spinqw.png the word GAPLESS, hardware.png the 50,000-shot caption.

In more detail, measured against the v3 manuscript before they were moved:

| file | what it shows that v3 does not say |
|---|---|
| `scaling.png` | not a figure: a block of v1 body text asserting the withdrawn selector-benchmark claim (relative-L1 0.042 ± 0.004 against 0.729, "a tenfold reduction") |
| `sqw.png` | the two fabricated `sqw_edges.dat` rows (q = 0 and q/π = 2 at 4.969), "charge (holon–antiholon) continuum" and "(same as A(k,ω))"; v3 says "band" and "(from A(k,ω), not measured here)" |
| `spinqw.png` | the words "GAPLESS" and "continuum", which v3 Sec. VIII does not use at this size |
| `hardware.png` | the v1 "Figure 15" caption with "50,000 computational-basis bitstrings"; v3 prints 3.5×10⁵ = 5×10⁴ on each of seven circuits |
| `resource.png` | v1 caption text: "the orbital invariant F1 is blind to it" |
| `bench.png`, `gallery.png` | thumbnails of `fig:bench` and `fig:gallery`, two figures v3 withdrew |
| `akw.png`, `hero.png`, `witness.png` | crops of the v1 renderings of figures that v3 still prints, replaced so that the gallery is one build |

The current thumbnails in `docs/img/` are rendered from the v3 figures by `docs/img/make_thumbnails.py`
(figure body only, no caption). Nothing reads the files in this folder; they are kept so that the v1
front page stays legible.
