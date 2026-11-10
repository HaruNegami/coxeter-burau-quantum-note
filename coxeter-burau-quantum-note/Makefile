PYTHON ?= python3
TEX = coxeter_burau_quantum

.PHONY: verify verify-coxeter verify-arrangements verify-mc pdf

verify: verify-coxeter verify-arrangements verify-mc

verify-mc:
	$(PYTHON) mc_arrangement_gates.py --verify --verify-symbolic
	$(PYTHON) mc_arrangement_gates.py --angles 1/10,1/10,1/10 --ell 1/10 --output rational_infinite.json
	$(PYTHON) mc_arrangement_gates.py --angles 1/4,1/4,1/4 --ell 1/4 --output resonant_qubit.json

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
