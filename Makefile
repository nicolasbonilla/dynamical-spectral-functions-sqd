# Reproducibility Makefile.  Requires: a TeX distribution (pdflatex) + `pip install -r requirements.txt`.
# Usage:  make verify | make figures | make paper | make all | make clean

PYTHON ?= python
LATEX  ?= pdflatex
PAPER  := paper/main

.PHONY: all ci verify check-figures check-coherence figures data paper clean help

help:
	@echo "make verify        - adversarial guardian: recompute the physics, then check every"
	@echo "                     deposited artefact against it (~15 s, numpy+scipy only)"
	@echo "make check-figures - regenerate every generated fragment in a scratch tree and diff"
	@echo "                     it against the committed one (writes nothing)"
	@echo "make check-coherence - does the manuscript still agree with ITSELF?  the range,"
	@echo "                     the size count, the shot budget and the guardian counts,"
	@echo "                     checked across sections and against the deposited grid"
	@echo "make figures       - regenerate the native pgfplots fragments IN PLACE"
	@echo "make data          - regenerate the deposited tables that have a generator (~4 min)"
	@echo "make ci            - what .github/workflows/ci.yml runs: verify + check-figures"
	@echo "make paper         - compile paper/main.tex -> paper/main.pdf"
	@echo "make all           - figures + paper"
	@echo "make reproduce     - execute notebooks/00_Reproduce_Everything.ipynb (needs pyscf;"
	@echo "                     NOT verified in the 2026-09-18 pass -- see docs/KNOWN_DISCREPANCIES.md)"
	@echo "make clean         - remove LaTeX aux files (keeps main.pdf)"
	@echo ""
	@echo "Where output goes: every compute script writes into THIS clone's data/ by default"
	@echo "(the path is derived from the script's own location), and --out PATH overrides it."
	@echo "Until 2026-09-18 eighteen of them wrote to the container path /w/ and the documented"
	@echo "pipeline did not repopulate data/.  See docs/KNOWN_DISCREPANCIES.md section 2."

# ---- the guardian (numpy+scipy only, ~15 s).  It RECOMPUTES the half-filled Hubbard ground
#      states and the geminal witness from first principles -- open chain F1=9.6047 (the resource
#      datasets and Figs. fig:decoupling, fig:master, fig:ladder) and periodic ring F1=9.3851
#      (fig:scaling); different systems, not variants -- and then checks ~4900 numeric assertions
#      over the deposited .json/.dat/.tex, the abstract included.  It exits non-zero and names the
#      defect.  Known open defects and open claims are printed, do not set the exit code, and are
#      the retraction agenda for v3.  `--fast` exists for editing loops and certifies NOTHING.
verify:
	$(PYTHON) src/verify.py

# ---- do the deposited generators still produce the committed figures?  (writes nothing)
check-figures:
	$(PYTHON) src/check_figures.py

# ---- does the manuscript still agree with ITSELF?  verify.py compares a printed
#      number against a deposited file; it does not compare a claim against a claim,
#      and on 2026-09-19 seven contradictions between sections survived a green run.
#      This is the guardian for that class.  Two synthetic controls per check.
check-coherence:
	$(PYTHON) src/check_coherence.py

# ---- regenerate the deposited tables that DO have a committed generator ----
data:
	$(PYTHON) src/make_sqw_edges.py
	$(PYTHON) src/charge_gap_ed.py
	@echo "OK: paper/figs/sqw_edges.dat and data/charge_gap.json regenerated from the deposited spectra."

# ---- regenerate the data-driven native \input fragments from the committed .json/.dat ----
# Verified 2026-09-18: all six commands exit 0 and every regenerated fragment is BYTE-IDENTICAL to the
# committed one. (Until then this target aborted on make_akw_sampled_fig.py, a %-format bug, and again
# on make_table.py, which wrote into a directory that does not exist. Both are fixed.)
# make_table.py writes rebuild/table_molecules.tex -- a scratch build directory it creates on demand.
# The table SHIPPED with the paper is paper/table_molecules.tex: same data rows, hand-added caption text.
# NOT in this target, deliberately: src/make_scaling_fig.py. It emits a DIFFERENT two-panel figure and
# running it would destroy panel (c) of the committed paper/figs/fig_scaling2_native.tex.
# See docs/KNOWN_DISCREPANCIES.md for that and for every other place the code and the deposit disagree.
# NOTE: the committed paper/figs/*.tex + *.pdf are AUTHORITATIVE and match main.pdf. Two fragments
# received manual annotation placement AFTER generation (resource-master annotation; the scaling
# panel (c)) — for those the committed .tex is the truth. The hardware panel (Fig. 15b) is fully
# regenerated below from data/hw_lucj_n2_result.json (8-seed mean + per-step std error bars).
# The .pdf figures (hero, bench, noise, gflow, gallery, magic-suite, circuit, S(q,w), S^zz, A(k,w))
# are standalone-compiled; see docs/REPRODUCE.md for the two-step recipe.
figures:
	$(PYTHON) src/make_decoupling_native.py
	$(PYTHON) src/make_method_fig_max.py
	$(PYTHON) src/make_akw_sampled_honest_fig.py
	$(PYTHON) src/make_noise_recovery_native.py
	$(PYTHON) src/make_hardware_hero.py
	$(PYTHON) src/make_table.py
	@echo "OK: data-driven fragments regenerated. Committed paper/figs/ remains authoritative."

# ---- reproduce EVERYTHING in one coherent notebook (narrated end-to-end pipeline) ----
# HONEST STATUS: this target was NOT verified in the 2026-09-18 pass.  It executes a notebook
# that imports pyscf, which is not installed in the environment the other targets were checked
# in.  Every other target in this file was run to completion and its result recorded in
# docs/KNOWN_DISCREPANCIES.md section 10.  This one has no measurement behind it.
reproduce:
	jupyter nbconvert --to notebook --execute --inplace notebooks/00_Reproduce_Everything.ipynb
	@echo "OK: full pipeline executed in notebooks/00_Reproduce_Everything.ipynb"

# ---- compile the preprint ----
paper:
	cd paper && $(LATEX) -interaction=nonstopmode main.tex >/dev/null && \
	            $(LATEX) -interaction=nonstopmode main.tex >/dev/null
	@echo "OK: paper/main.pdf built."

all: figures paper

# `make ci` is what .github/workflows/ci.yml runs, in the order it runs it.
ci: verify check-figures check-coherence

clean:
	rm -f paper/*.aux paper/*.log paper/*.out paper/*.toc paper/*.fls paper/*.fdb_latexmk \
	      paper/figs/*.aux paper/figs/*.log
	@echo "cleaned LaTeX aux files."
