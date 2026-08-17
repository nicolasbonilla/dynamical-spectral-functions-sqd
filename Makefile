# Reproducibility Makefile.  Requires: a TeX distribution (pdflatex) + `pip install -r requirements.txt`.
# Usage:  make verify | make figures | make paper | make all | make clean

PYTHON ?= python
LATEX  ?= pdflatex
PAPER  := paper/main

.PHONY: all verify figures paper clean help

help:
	@echo "make verify   - fast exact-diagonalization smoke test (the geminal witness, seconds)"
	@echo "make figures  - regenerate every native pgfplots fragment from committed data"
	@echo "make paper    - compile paper/main.tex -> paper/main.pdf"
	@echo "make all      - figures + paper"
	@echo "make clean    - remove LaTeX aux files (keeps main.pdf)"

# ---- fast reproducibility check (numpy+scipy only, seconds): F1(L=6)=9.3851 + the geminal witness ----
verify:
	$(PYTHON) src/verify.py

# ---- regenerate the data-driven native \input fragments from the committed .json/.dat ----
# NOTE: the committed paper/figs/*.tex + *.pdf are AUTHORITATIVE and match main.pdf. Two fragments
# received manual annotation placement AFTER generation (resource-master annotation; the scaling
# panel (c)) — for those the committed .tex is the truth. The hardware panel (Fig. 15b) is fully
# regenerated below from data/hw_lucj_n2_result.json (8-seed mean + per-step std error bars).
# The .pdf figures (hero, bench, noise, gflow, gallery, magic-suite, circuit, S(q,w), S^zz, A(k,w))
# are standalone-compiled; see docs/REPRODUCE.md for the two-step recipe.
figures:
	$(PYTHON) src/make_decoupling_native.py
	$(PYTHON) src/make_method_fig_max.py
	$(PYTHON) src/make_akw_sampled_fig.py
	$(PYTHON) src/make_noise_recovery_native.py
	$(PYTHON) src/make_hardware_hero.py
	$(PYTHON) src/make_table.py
	@echo "OK: data-driven fragments regenerated. Committed paper/figs/ remains authoritative."

# ---- reproduce EVERYTHING in one coherent notebook (narrated end-to-end pipeline) ----
reproduce:
	jupyter nbconvert --to notebook --execute --inplace notebooks/00_Reproduce_Everything.ipynb
	@echo "OK: full pipeline executed in notebooks/00_Reproduce_Everything.ipynb"

# ---- compile the preprint ----
paper:
	cd paper && $(LATEX) -interaction=nonstopmode main.tex >/dev/null && \
	            $(LATEX) -interaction=nonstopmode main.tex >/dev/null
	@echo "OK: paper/main.pdf built."

all: figures paper

clean:
	rm -f paper/*.aux paper/*.log paper/*.out paper/*.toc paper/*.fls paper/*.fdb_latexmk \
	      paper/figs/*.aux paper/figs/*.log
	@echo "cleaned LaTeX aux files."
