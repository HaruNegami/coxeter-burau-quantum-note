# Generalized Burau representations and quantum gates

**A note by Haru Negami — version 2 (10 November 2026).**

[Note (PDF)](coxeter_burau_quantum.pdf) · [TeX](coxeter_burau_quantum.tex) ·
[日本語](README_ja.md) · [Changes](CHANGELOG.md)

Version 2 adds the construction from a general affine hyperplane arrangement
through rank-one middle convolution (MC) to unitary quantum gates, using
Negami's prior invariant-form and signature results. The earlier Coxeter-type
results and arrangement appendices remain in the same note.

| Construction | Location | Scope |
| --- | --- | --- |
| Squier's Coxeter–Burau representations | Sections 1–6 | Artin–Tits representations, finite-type positivity, E6/E7/E8 closures and E8 three-qubit gates |
| Rank-one directional MC | Section 7 | Actual MC quotient on the directional closure of an affine arrangement; form, signature, gates, density and finite-image criteria |
| Additive convolution and based monodromy | Appendices A–H | Arbitrary-rank formulas, middle quotients, Burau pullbacks, six-line example and descent |

The two constructions have different dimensions and parameter conditions.
No general identification between them is assumed. Squier's representations
and forms are known. Section 6 derives the exceptional Coxeter closure
conclusions from Freedman–Larsen–Wang; a direct rank-two proof is also included.
Section 7 explicitly specializes the author's prior theorems and supplies the
additional arguments. No priority claim is made.

## Files

| File | Purpose |
| --- | --- |
| `coxeter_burau_quantum.pdf` / `.tex` | Complete v2 note, retaining its filename |
| `mc_arrangement_gates.py` | MC matrices, invariant form, radical quotient, common unitarization, rational finite-image test |
| `mc_gate_example.json` | Three MC qutrit gates and verification results |
| `resonant_qubit.json` | Rank drops from three to two at product resonance |
| `rational_infinite.json` | Positive-definite rational MC example with infinite image |
| `verify_coxeter_gates.py` | Original exact E-type relations and numerical gate construction |
| `verification_report.json` / `e8_gates_theta_pi_over_60.json` | Original Coxeter results and three-qubit gates |
| `legacy_arrangement/` | Original four convolution, quotient, Fox-calculus and based-monodromy scripts |
| `versions/v1/` | Unchanged initial PDF and TeX |
| `CHANGELOG.md` / `RELEASE_NOTES_v2.md` | Version history and release description |

## Install and reproduce

Use Python 3.11 or later, NumPy and SymPy. From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python mc_arrangement_gates.py --verify --verify-symbolic
```

On Windows, activate with `.venv\Scripts\activate`.
`requirements-lock.txt` records the versions used for the bundled results.

The default MC example uses the three lines `y=0, y=x, y=1-x`, adds collision
lines `x=0, x=1, x=1/2`, and has angles `(1/20,1/15,1/12)` and
`ell=sqrt(2)/20`. Its actual MC quotient is three-dimensional. The exported
fiber generators give qutrit gates; full unitary closure follows from the proof.
Numerical entries are approximations, with exact input expressions in metadata.

```bash
# Rational MC examples
python mc_arrangement_gates.py --angles 1/10,1/10,1/10 --ell 1/10 --output rational_infinite.json
python mc_arrangement_gates.py --angles 1/4,1/4,1/4 --ell 1/4 --output resonant_qubit.json
# Original Coxeter example
python verify_coxeter_gates.py
```

The Coxeter script checks the E6/E7/E8 relations exactly over integer Laurent
polynomials, checks 12 numerical specializations, and exports E8 at
`theta/pi=1/60`. A different export can be requested as follows:

```bash
python verify_coxeter_gates.py --rank 7 --theta-over-pi 1/40 \
  --output verification_e7.json --export e7_gates.json
```

Its condition `abs(theta/pi)<1/h`, with `h=12,18,30`, is a Coxeter condition,
not the MC condition below. With GNU Make:

```bash
make verify              # Both constructions and all four original appendix scripts
make verify-mc           # MC checks and three exports
make verify-coxeter
make verify-arrangements
make pdf                 # Compile the complete v2 note twice
```

Without Make, run the commands above and each of
`legacy_arrangement/arrangement_burau.py`, `general_convolution.py`,
`fox_burau.py`, and `monodromy_matrices.py` with Python. To compile directly,
run `pdflatex -interaction=nonstopmode -halt-on-error coxeter_burau_quantum.tex`
twice. The compiled PDF is included.

## MC parameters and based loops

For `0<=a_i<1` and `0<ell<1`, let `p` count nonzero `a_i`,
`T=ell+sum(a_i)`, `m=floor(T)`, and `delta=1` if T is an integer (zero
otherwise). The MC quotient has dimension `p-delta` and signature
`(p-m,m-delta)`. For a nonzero quotient, the form is positive definite
exactly for `0<T<=1`, and negative definite for `p<=T<p+1`.
Either sign yields unitary gates after a possible overall sign change.

For dimension at least two, definiteness and irrationality of every active
`ell+a_i` imply full unitary closure and non-Clifford fiber generators.
This sufficient condition does not classify every rational or mixed-phase
closure. The rational finite-image test assumes finite-image base character.
Numerical Weyl-normalizer residuals in the exported basis are diagnostics,
not exact Clifford certificates.

The default export contains only actual fiber generators. Other loops can be
evaluated by supplying their based pure root-braid words with `--base-words
path.json`. Input format:

```json
{"words":{"chosen_base_loop":[1,1]},"phases":{"chosen_base_loop":"0"}}
```

This is a format example, not a meridian prescribed for every arrangement.
Words are chronological, the next matrix multiplies on the left, and colors
are swapped at every crossing. Non-pure words are rejected. Clockwise fiber
loops follow Negami's algebraic MC convention; compare this with the
counterclockwise convention of the earlier based-Burau appendix.

An arbitrary arrangement need not have a Dynkin type. MC may add collision
hyperplanes; descent requires trivial monodromy around each added hyperplane.
Reflection arrangements give pure Artin groups; extension to the full Artin
group requires equivariance. A generic E-type MC rank is not automatically
6, 7, or 8.

The text and code were prepared with substantial assistance from OpenAI
ChatGPT. The arguments are provided for review; symbolic and numerical checks
do not establish novelty or certify every proof.

## References

- C. C. Squier, *Matrix representations of Artin groups* (1988).
- Freedman–Larsen–Wang, [The two-eigenvalue problem](https://arxiv.org/abs/math/0103200).
- H. Negami, [KLM and multiplicative MC](https://arxiv.org/abs/2503.14840), Theorems 8–9.
- H. Negami, [KLM signature formula](https://arxiv.org/abs/2609.11793), Theorem 1.1.
- H. Negami, [Finite images and algebraic KZ-type solutions](https://arxiv.org/abs/2609.25805).
- H. Negami, [Quantum gates from MC of twisted Burau](https://arxiv.org/abs/2610.00293).
- Hiroe–Negami, [Long–Moody construction and Katz MC](https://arxiv.org/abs/2303.05770).
- Adachi–Hiroe, [MC for affine hyperplane arrangements](https://arxiv.org/abs/2601.00544).

See the note for the complete bibliography and precise attribution.
