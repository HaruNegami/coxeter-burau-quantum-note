PYTHON ?= python3
TEX = coxeter_burau_quantum

.PHONY: verify verify-coxeter verify-arrangements pdf

verify: verify-coxeter verify-arrangements

verify-coxeter:
	$(PYTHON) verify_coxeter_gates.py

verify-arrangements:
	$(PYTHON) legacy_arrangement/arrangement_burau.py
	$(PYTHON) legacy_arrangement/general_convolution.py
	$(PYTHON) legacy_arrangement/fox_burau.py
	$(PYTHON) legacy_arrangement/monodromy_matrices.py

pdf:
	pdflatex -interaction=nonstopmode -halt-on-error $(TEX).tex
	pdflatex -interaction=nonstopmode -halt-on-error $(TEX).tex
