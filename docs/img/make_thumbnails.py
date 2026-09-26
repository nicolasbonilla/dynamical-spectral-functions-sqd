# -*- coding: utf-8 -*-
"""Render the README gallery thumbnails from the committed figures -- figure body only, no caption.

    python docs/img/make_thumbnails.py        # needs PyMuPDF (pip install pymupdf)

Four thumbnails are the whole page of a committed figure PDF in paper/figs/.  The other
four are cropped out of paper/main.pdf: the page is found by searching for the caption tag
("FIG. 5.", ...), and the crop is the union of every drawing and text block between the
running head and the top of that caption, so no caption text enters the thumbnail.
Written 2026-09-26, when the ten thumbnails cropped from the v1 PDF were moved to
_superseded/v1_readme_img/.
"""
import os
import sys

import fitz  # PyMuPDF

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DPI = 150
ZOOM = fitz.Matrix(DPI / 72.0, DPI / 72.0)

# thumbnail -> committed figure PDF (the whole page is the figure body)
WHOLE = {
    "akw.png": "paper/figs/fig_akw_native.pdf",      # Fig. 4, fig:lattice
    "sqw.png": "paper/figs/fig_sqw.pdf",             # Fig. 7, fig:sqw
    "spinqw.png": "paper/figs/fig_spinqw.pdf",       # Fig. 8, fig:spin
    "hero.png": "paper/figs/fig_hero.pdf",           # Fig. S6, fig:hero
}
# thumbnail -> caption tag of a figure in paper/main.pdf
CROPPED = {
    "scaling.png": "FIG. 5.",       # fig:scaling
    "hardware.png": "FIG. 9.",      # fig:heron
    "witness.png": "FIG. S1.",      # fig:witness
    "resource.png": "FIG. S3.",     # fig:master
}
HEAD = 45.0     # pt: below this is the running head (page number), never part of a figure
PAD = 4.0       # pt of margin around the figure body


def crop_rect(page, tag):
    blocks = page.get_text("blocks")
    caps = [b for b in blocks if b[4].strip().startswith(tag)]
    if len(caps) != 1:
        return None
    cap_y0 = caps[0][1]
    rects = [fitz.Rect(b[:4]) for b in blocks if b[1] > HEAD and b[3] <= cap_y0 + 0.5]
    rects += [d["rect"] for d in page.get_drawings()
              if d["rect"].y0 > HEAD and d["rect"].y1 <= cap_y0 + 0.5]
    rects += [fitz.Rect(i["bbox"]) for i in page.get_image_info()
              if i["bbox"][1] > HEAD and i["bbox"][3] <= cap_y0 + 0.5]
    if not rects:
        return None
    r = rects[0]
    for x in rects[1:]:
        r |= x
    r = fitz.Rect(r.x0 - PAD, r.y0 - PAD, r.x1 + PAD, min(r.y1 + PAD, cap_y0 - 1.0))
    return r & page.rect


def main():
    out = []
    for name, rel in WHOLE.items():
        doc = fitz.open(os.path.join(REPO, rel))
        doc[0].get_pixmap(matrix=ZOOM, alpha=False).save(os.path.join(HERE, name))
        out.append((name, rel, "whole page"))
    doc = fitz.open(os.path.join(REPO, "paper", "main.pdf"))
    for name, tag in CROPPED.items():
        hit = None
        for page in doc:
            r = crop_rect(page, tag)
            if r is not None:
                hit = (page, r)
                break
        if hit is None:
            print("caption %r not found in paper/main.pdf" % tag)
            return 1
        page, r = hit
        page.get_pixmap(matrix=ZOOM, clip=r, alpha=False).save(os.path.join(HERE, name))
        out.append((name, "paper/main.pdf p. %d" % (page.number + 1),
                    "clip %.1f %.1f %.1f %.1f pt" % tuple(r)))
    for row in out:
        print("%-13s <- %-34s %s" % row)
    return 0


if __name__ == "__main__":
    sys.exit(main())
