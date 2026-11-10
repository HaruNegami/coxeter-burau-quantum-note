#!/usr/bin/env python3
"""Rank-one Haraoka/Adachi-Hiroe convolution for affine arrangements.

Exact rational geometry; symbolic scalar weights. Coordinate y selects the
projection direction. Output matrices are residues, not based monodromy.
Reference: arXiv:2601.00544v4, (3.2)--(3.6).
Run: python arrangement_burau.py
"""
from dataclasses import dataclass
from itertools import combinations
import json
from pathlib import Path
import sympy as s


def zero(M):
    return all(s.simplify(v) == 0 for v in M)


def affine_row(f, variables):
    p = s.Poly(s.expand(f), *variables)
    if p.total_degree() > 1:
        raise ValueError("Hyperplanes must be affine linear.")
    row = [p.coeff_monomial(v) for v in variables] + [p.coeff_monomial(1)]
    if not all(c.is_Rational for c in row):
        raise ValueError("Use exact rational geometry; weights may be symbolic.")
    return tuple(row)


def normalize(f, variables):
    row = affine_row(f, variables)
    pivot = next((c for c in row[:-1] if c != 0), None)
    if pivot is None:
        if row[-1] == 0:
            raise ValueError("Coincident hyperplanes / identically zero difference.")
        return None  # a nonzero constant has no zero hyperplane
    return tuple(s.cancel(c / pivot) for c in row)


def expression(row, variables):
    return sum(a * v for a, v in zip(row[:-1], variables)) + row[-1]


@dataclass
class Convolution:
    variables: tuple
    y: object
    graphs: list
    weights: list
    parameter: object
    equations: list
    residues: list
    added: list
    original: set


def convolve(equations, weights, variables, y, parameter):
    """Return full convolution on the minimal y-closed arrangement.

    Geometry is supplied by equations f=0. The weight is the scalar residue
    of dlog(f), hence rescaling f makes no change. Coordinate y may be any
    member of variables. For a general direction first change coordinates.
    """
    variables = tuple(variables)
    if y not in variables or len(equations) != len(weights):
        raise ValueError("Supply one weight per equation and y in variables.")
    originals = set()
    transverse, parallel = [], {}
    for f, a in zip(equations, weights):
        key = normalize(f, variables)
        if key is None or key in originals:
            raise ValueError("Supply distinct nonempty hyperplanes.")
        originals.add(key)
        a = s.sympify(a)
        cy = s.diff(f, y)
        if cy == 0:
            parallel[key] = a
        else:
            graph = s.expand(y - f / cy)
            transverse.append((graph, a, key))
    n = len(transverse)
    lam = s.sympify(parameter)
    residues = {}
    for i, (g, a, key) in enumerate(transverse):
        C = s.zeros(n)
        for j, (_, b, _) in enumerate(transverse):
            C[i, j] = b
        C[i, i] += lam
        residues[key] = C
    for key, a in parallel.items():
        residues[key] = a * s.eye(n)
    for i, j in combinations(range(n), 2):
        gi, ai, _ = transverse[i]
        gj, aj, _ = transverse[j]
        key = normalize(gi - gj, variables)
        if key is None:
            continue
        if key not in residues:
            residues[key] = s.zeros(n)
        residues[key][i, i] += aj
        residues[key][i, j] -= aj
        residues[key][j, i] -= ai
        residues[key][j, j] += ai
    keys = list(residues)
    return Convolution(variables, y, [v[0] for v in transverse],
                       [v[1] for v in transverse], lam,
                       [expression(k, variables) for k in keys],
                       [residues[k] for k in keys],
                       [i for i, k in enumerate(keys) if k not in originals],
                       originals)


def middle_quotient(model, substitutions=None):
    """Compute K+L and the induced matrices after parameter specialization.

    Without substitutions, nullspaces are over the symbolic fraction field:
    generic rank only. For exceptional strata supply substitutions here.
    Returns Q, R, quotient residues, and a basis matrix N for K+L, with
    Q*R=I, ker(Q)=im(N), and Q*C=B*Q for every residue.
    """
    sub = substitutions or {}
    n = len(model.weights)
    cs = [C.subs(sub).applyfunc(s.simplify) for C in model.residues]
    lam = s.simplify(model.parameter.subs(sub))
    a = [s.simplify(v.subs(sub)) for v in model.weights]
    if n == 0:
        return s.zeros(0), s.zeros(0), cs, s.zeros(0)
    vectors = [s.eye(n)[:, i] for i, v in enumerate(a) if v == 0]
    fibre_sum = s.ones(n, 1) * s.Matrix([a]) + lam * s.eye(n)
    vectors += fibre_sum.nullspace()
    N0 = s.Matrix.hstack(*vectors) if vectors else s.zeros(n, 0)
    basis = N0.columnspace()
    N = s.Matrix.hstack(*basis) if basis else s.zeros(n, 0)
    T = N
    for i in range(n):
        candidate = T.row_join(s.eye(n)[:, i])
        if candidate.rank() > T.cols:
            T = candidate
    Q = T.inv()[N.cols:, :]
    R = T[:, N.cols:]
    out = [(Q * C * R).applyfunc(s.simplify) for C in cs]
    for C, B in zip(cs, out):
        assert zero(Q * C - B * Q), "Middle subspace not invariant"
    assert zero(Q * R - s.eye(Q.rows))
    return Q, R, out, N


def kohno_check(equations, matrices, variables):
    """Check all codimension-two holonomy relations, including affine cases."""
    rows = [affine_row(f, variables) for f in equations]
    flats = set()
    for i, j in combinations(range(len(rows)), 2):
        pair = s.Matrix([rows[i], rows[j]])
        if pair[:, :-1].rank() != 2:
            continue  # parallel pair, no codimension-two intersection
        incident = tuple(k for k, r in enumerate(rows)
                         if pair.col_join(s.Matrix([r])).rank() == 2)
        flats.add(incident)
    count = 0
    for incident in flats:
        total = s.zeros(matrices[0].rows)
        for i in incident:
            total += matrices[i]
        for i in incident:
            assert zero(matrices[i] * total - total * matrices[i]), incident
            count += 1
    return {"codimension_two_flats": len(flats), "commutators_checked": count}


def reduced_burau(strands, index, t):
    """sigma_index, 1-based, on C^strands / C*(1,...,1).

    Unreduced block [[1-t,t],[1,0]]. Generic quotient model of reduced Burau.
    Products act on columns. This is a based braid matrix, not a residue.
    """
    if not 1 <= index < strands:
        raise ValueError("Require 1 <= index < strands")
    U = s.eye(strands)
    i = index - 1
    U[i:i+2, i:i+2] = s.Matrix([[1-t, t], [1, 0]])
    Q = s.eye(strands - 1).row_join(-s.ones(strands - 1, 1))
    R = s.eye(strands)[:, :strands-1]
    return Q * U * R


def burau_word(strands, word, t):
    """Word [1,-2,1] means sigma_1 sigma_2^-1 sigma_1."""
    M = s.eye(strands - 1)
    for i in word:
        B = reduced_burau(strands, abs(i), t)
        M = M * (B if i > 0 else B.inv())
    return M.applyfunc(s.factor)


def check_equal_weight_identification(n):
    """Verify constant conjugacy to the standard A_n reflection connection."""
    # On C^(n+1)/C*1 use u_i=w_i-w_(n+1). Then T_(i,n+1)
    # has column i equal to (1,...,1)+e_i; AH has its transpose.
    S = s.eye(n) + s.ones(n)
    for i in range(n):
        ei = s.eye(n)[:, i]
        C = ei * (s.ones(1, n) + ei.T)
        T = (s.ones(n, 1) + ei) * ei.T
        assert zero(S * C * S.inv() - T)
    for i, j in combinations(range(n), 2):
        v = s.eye(n)[:, i] - s.eye(n)[:, j]
        T = v * v.T
        assert zero(S * T * S.inv() - T)


def main():
    x, y, z = s.symbols('x y z')
    a, b, c, d, lam, t = s.symbols('a b c d lambda t')
    triangle = convolve([y, y-x, y+x-1], [a,b,c], (x,y), y, lam)
    checks = {"triangle_generic": kohno_check(triangle.equations,
                                              triangle.residues, (x,y))}
    Q, R, reduced, N = middle_quotient(triangle, {lam: -a-b-c})
    assert Q.rows == 2
    checks['triangle_middle'] = kohno_check(triangle.equations, reduced, (x,y))
    # λ=0 and a zero residue test the actual K+L rather than generic rank.
    assert middle_quotient(triangle, {lam: 0})[0].rows == 1
    assert middle_quotient(triangle, {a: 0, lam: -b-c})[0].rows == 1
    empty = convolve([x], [a], (x,y), y, lam)
    assert middle_quotient(empty)[0].rows == 0
    for equations, weights, variables in [
        ([y,y-x,y+x,x], [a,b,c,d], (x,y)),  # triple collision, existing pole
        ([y,y-1,y-x], [a,b,c], (x,y)),      # noncolliding parallel fibres
        ([y,y-x,y-z], [a,b,c], (x,z,y)),    # three ambient coordinates
    ]:
        m = convolve(equations,weights,variables,y,lam)
        checks[str(equations)] = kohno_check(m.equations,m.residues,variables)
    for n in [2,3,4]:
        check_equal_weight_identification(n)
    for n in [3,4,5]:
        bs = [reduced_burau(n,i,t) for i in range(1,n)]
        for i in range(n-2):
            assert zero(bs[i]*bs[i+1]*bs[i]-bs[i+1]*bs[i]*bs[i+1])
        for i,j in combinations(range(n-1),2):
            if abs(i-j)>1:
                assert zero(bs[i]*bs[j]-bs[j]*bs[i])
        full_twist = burau_word(n, list(range(1,n))*n,t)
        assert zero(full_twist-t**n*s.eye(n-1))
    output = {
        'convention': 'dV = sum C_H dlog(f_H) V; residues are not based monodromy',
        'equations': [str(f) for f in triangle.equations],
        'added_equations': [str(triangle.equations[i]) for i in triangle.added],
        'generic_residues': [s.sstr(C) for C in triangle.residues],
        'middle_parameter': '-a-b-c', 'middle_projection':s.sstr(Q),
        'middle_residues': [s.sstr(C) for C in reduced],
        'burau_B3_sigma1':s.sstr(reduced_burau(3,1,t)),
        'burau_B3_sigma2':s.sstr(reduced_burau(3,2,t)),
        'checks':checks,
        'extra_checks':'Exceptional quotient ranks, A_n conjugacy n=2,3,4, braid relations and full twists B3,B4,B5: PASS',
    }
    path = Path(__file__).with_name('results.json')
    path.write_text(json.dumps(output, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(output,ensure_ascii=False,indent=2))


if __name__ == '__main__':
    main()
