# Generalized Burau representations and quantum gates

**A note by Haru Negami.**

[Read the note (PDF)](coxeter_burau_quantum.pdf) ·
[LaTeX source](coxeter_burau_quantum.tex) ·
[日本語の説明](README_ja.md)

This repository collects explicit generalized Burau matrices, invariant
Hermitian forms, and their interpretation as quantum gates. It includes
calculations for finite Coxeter types E6, E7, and E8, and separate
middle-convolution calculations for affine hyperplane arrangements.

The note is an exposition with explicit calculations. Squier's representations
and invariant forms are existing constructions. The exceptional-type closure
conclusions follow from the two-eigenvalue classification of
Freedman–Larsen–Wang; Section 6 checks the hypotheses and gives the reduction.
An alternative rank-two proof is included.

## What is included

| File | Purpose |
| --- | --- |
| `coxeter_burau_quantum.pdf` / `.tex` | Note, proofs, conventions, and references |
| `verify_coxeter_gates.py` | Exact E-type relation checks and numerical unitary gate construction |
| `verification_report.json` | Results of the Coxeter checks |
| `e8_gates_theta_pi_over_60.json` | Eight explicit three-qubit gates and state-preparation matrices |
| `legacy_arrangement/` | Arrangement convolution, quotient, Fox-calculus, and based-monodromy scripts |
| `CITATION.cff` | Citation metadata for this note |

## Install and reproduce

Use Python 3.10 or later. From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python verify_coxeter_gates.py
```

On Windows, activate the environment with `.venv\Scripts\activate`.
The main script requires NumPy; the arrangement scripts also require SymPy.
`requirements-lock.txt` records the versions used for the bundled check results.

The default command verifies the invariant-form identities and Artin relations
for E6, E7, and E8 exactly over integer Laurent polynomials. It also checks 12
numerical specializations and exports the E8 gates at `theta/pi = 1/60`.
Generated matrix entries are floating-point approximations. Their defining
formulas and the exact rational phase appear in the note and JSON metadata.

To choose another example:

```bash
python verify_coxeter_gates.py --rank 7 --theta-over-pi 1/40 \
  --output verification_e7.json --export e7_gates.json
```

The parameter must satisfy `abs(theta/pi) < 1/h`, where
`h = 12, 18, 30` for E6, E7, E8 respectively. The script checks all three types
regardless of which rank is selected for the gate export.

Run the arrangement checks separately:

```bash
python legacy_arrangement/arrangement_burau.py
python legacy_arrangement/general_convolution.py
python legacy_arrangement/fox_burau.py
python legacy_arrangement/monodromy_matrices.py
```

These commands regenerate the corresponding reports in `legacy_arrangement/`.
They check middle quotients, descent, Fox derivatives, pure-braid formulas,
and the six-line example. The geometric path-tracking check is numerical.
The symbolic relation checks are exact.

With GNU Make, `make verify` runs all five scripts. No large group enumeration
is needed. The density and finite-image conclusions rely on the written proofs.

## Build the PDF

With a LaTeX installation including `pdflatex` and the packages used in the source:

```bash
pdflatex -interaction=nonstopmode -halt-on-error coxeter_burau_quantum.tex
pdflatex -interaction=nonstopmode -halt-on-error coxeter_burau_quantum.tex
```

Alternatively, run `make pdf`. The compiled PDF is included for readers who do
not have LaTeX installed.

## Scope

The generalized Burau formula is defined for finite-rank Coxeter matrices.
The compact-closure and quantum-universality statements in the note have the
stated finite-Coxeter and parameter hypotheses. The code explicitly constructs
the E6, E7, and E8 families.

An arbitrary arrangement need not have a Dynkin type. Its directional
middle-convolution construction may add hyperplanes, and descent to the original
complement is a separate condition. The note does not identify every generalized
Coxeter–Burau representation with a middle-convolution output.

## References

- C. C. Squier, *Matrix representations of Artin groups* (1988).
- M. H. Freedman, M. J. Larsen, Z. Wang,
  [The two-eigenvalue problem and density of Jones representation of braid groups](https://arxiv.org/abs/math/0103200).
- H. Negami,
  [Quantum gates from the middle convolution of twisted Burau representations](https://arxiv.org/abs/2610.00293).
- S. Adachi and K. Hiroe,
  [Middle convolution for affine hyperplane arrangements and Aomoto realizability](https://arxiv.org/abs/2601.00544).
- K. Hiroe,
  [Middle convolution for Lie algebra representations](https://arxiv.org/abs/2605.09828).

See the note for the complete bibliography and precise attribution of each result.
