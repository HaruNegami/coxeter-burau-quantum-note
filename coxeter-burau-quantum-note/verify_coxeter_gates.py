#!/usr/bin/env python3
"""Small, reproducible Coxeter--Burau checks. Requires Python 3 + NumPy.

The exact checks use integer Laurent polynomials, without SymPy.
Numerical checks are diagnostics and are NOT density/finite-image proofs.
Run: python verify_coxeter_gates.py --output verification_report.json
Export a gate set: python verify_coxeter_gates.py --export e8_gates.json
"""
import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
import numpy as np

DATA = {
    6: (12, [1, 4, 5, 7, 8, 11]),
    7: (18, [1, 5, 7, 9, 11, 13, 17]),
    8: (30, [1, 7, 11, 13, 17, 19, 23, 29]),
}

def edges(n):
    return {(i, i + 1) for i in range(n - 2)} | {(2, n - 1)}

def padd(a, b):
    c = a.copy()
    for k, v in b.items():
        c[k] = c.get(k, 0) + v
    return {k: v for k, v in c.items() if v}

def pmul(a, b):
    c = {}
    for i, x in a.items():
        for j, y in b.items():
            c[i + j] = c.get(i + j, 0) + x * y
    return {k: v for k, v in c.items() if v}

def eye(n):
    return [[{0: 1} if i == j else {} for j in range(n)] for i in range(n)]

def mm(a, b):
    n = len(a)
    c = [[{} for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for k in range(n):
            if a[i][k]:
                for j in range(n):
                    if b[k][j]:
                        c[i][j] = padd(c[i][j], pmul(a[i][k], b[k][j]))
    return c

def star(a):
    return [[{-k: v for k, v in a[j][i].items()} for j in range(len(a))]
            for i in range(len(a))]

def exact_checks(n):
    ed = edges(n)
    h = [[{1: 1, -1: 1} if i == j else
          {0: -1} if tuple(sorted((i, j))) in ed else {}
          for j in range(n)] for i in range(n)]
    generators = []
    for i in range(n):
        r = eye(n)
        r[i] = [{2: -1} if i == j else
                {1: 1} if tuple(sorted((i, j))) in ed else {}
                for j in range(n)]
        assert mm(mm(star(r), h), r) == h, ('invariance', n, i)
        generators.append(r)
    for i in range(n):
        for j in range(i + 1, n):
            a, b = generators[i], generators[j]
            if (i, j) in ed:
                assert mm(mm(a, b), a) == mm(mm(b, a), b), ('braid', n, i, j)
            else:
                assert mm(a, b) == mm(b, a), ('commute', n, i, j)
    return {'type': f'E{n}', 'ring': 'Z[t,t^-1]',
            'invariant_form_identities': n,
            'Artin_relations': n * (n - 1) // 2, 'passed': True}

def construct(n, theta):
    hnum, _ = DATA[n]
    if abs(theta) >= math.pi / hnum:
        raise ValueError('Parameter must be strictly inside the positive interval.')
    a = np.zeros((n, n))
    for i, j in edges(n):
        a[i, j] = a[j, i] = 1
    h = 2 * np.cos(theta) * np.eye(n) - a
    vals, vecs = np.linalg.eigh(h)
    s = (vecs * np.sqrt(vals)) @ vecs.T
    sinv = (vecs * (1 / np.sqrt(vals))) @ vecs.T
    t = np.exp(1j * theta)
    lam = -t * t
    p = s / np.sqrt(2 * np.cos(theta))
    rs, us, vs = [], [], []
    target = np.eye(n)[:, -1]
    for i in range(n):
        r = np.eye(n, dtype=complex)
        r[i, :] -= t * h[i, :]
        u = np.eye(n) + (lam - 1) * np.outer(p[:, i], p[:, i])
        w = target - p[:, i]
        v = np.eye(n) if np.linalg.norm(w) < 1e-14 else (
            np.eye(n) - 2 * np.outer(w, w) / np.dot(w, w))
        rs.append(r); us.append(u); vs.append(v)
    return h, s, sinv, p, rs, us, vs, lam

def norm(x):
    return float(np.linalg.norm(x, ord='fro'))

def numerical_checks(n, theta):
    hnum, exponents = DATA[n]
    h, s, sinv, p, rs, us, vs, lam = construct(n, theta)
    expected = np.sort(2 * np.cos(theta) -
                       2 * np.cos(np.pi * np.array(exponents) / hnum))
    d = np.eye(n, dtype=complex); d[-1, -1] = lam
    report = {
        'type': f'E{n}', 'theta_over_pi': theta / np.pi,
        'smallest_form_eigenvalue': float(np.linalg.eigvalsh(h)[0]),
        'spectrum_error': float(np.max(abs(np.linalg.eigvalsh(h) - expected))),
        'invariance_error': max(norm(r.conj().T @ h @ r - h) for r in rs),
        'unitarity_error': max(norm(u.conj().T @ u - np.eye(n)) for u in us),
        'common_conjugation_error': max(norm(s @ r @ sinv - u)
                                       for r, u in zip(rs, us)),
        'determinant_error': max(float(abs(np.linalg.det(u) - lam)) for u in us),
        'trace_error': max(float(abs(np.trace(u) - (n - 1 + lam))) for u in us),
        'phase_gate_factorization_error': max(norm(v @ d @ v.T - u)
                                             for v, u in zip(vs, us)),
        'encoded_three_qubit_trace_squared': float(abs(7 + lam) ** 2),
        'nonscalar_three_qubit_Clifford_trace_bound': 32,
    }
    for key, value in report.items():
        if key.endswith('_error'):
            assert value < 2e-10, (key, value)
    assert report['encoded_three_qubit_trace_squared'] > 32
    return report

def complex_array(x):
    return {'real': np.asarray(x).real.tolist(), 'imag': np.asarray(x).imag.tolist()}

def export_gates(path, n, theta_fraction):
    theta = float(theta_fraction) * np.pi
    h, s, sinv, p, rs, us, vs, lam = construct(n, theta)
    phase_fraction = (Fraction(1, 2) + theta_fraction) % 1
    order = phase_fraction.denominator
    data = {
        'author': 'Haru Negami', 'type': f'E{n}',
        'theta_over_pi_exact': str(theta_fraction),
        'lambda_phase_over_2pi_exact': str(phase_fraction),
        'generator_order_exact': order,
        'edges_one_based': [[i + 1, j + 1] for i, j in sorted(edges(n))],
        'orthonormal_basis': 'Coordinates after S=positive square root of H',
        'three_qubit_basis_if_E8': ['000','001','010','011','100','101','110','111'],
        'H': h.tolist(), 'S': s.tolist(),
        'unit_vectors_columns': p.tolist(),
        'gates': [complex_array(u) for u in us],
        'state_preparation_Householders': [v.tolist() for v in vs],
        'note': 'Floating-point matrix entries; exact definitions are in the paper. '
                'Order is computed from the exact rational phase, not estimated numerically.',
        'order_error': max(norm(np.linalg.matrix_power(u, order) - np.eye(n)) for u in us),
    }
    Path(path).write_text(json.dumps(data, indent=2) + '\n')

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', default='verification_report.json')
    ap.add_argument('--export', default='e8_gates_theta_pi_over_60.json')
    ap.add_argument('--rank', type=int, choices=DATA, default=8)
    ap.add_argument('--theta-over-pi', default='1/60', help='Exact rational number, e.g. 1/60')
    args = ap.parse_args()
    report = {'exact': [exact_checks(n) for n in DATA], 'numerical': []}
    for n, (hnum, _) in DATA.items():
        for factor in [0, 0.5, -0.5, 0.9]:
            report['numerical'].append(numerical_checks(n, factor * np.pi / hnum))
    report['warning'] = 'Density and finite-image classification use the written proofs; '
    report['warning'] += 'no numerical density inference or finite group enumeration is performed.'
    Path(args.output).write_text(json.dumps(report, indent=2) + '\n')
    export_gates(args.export, args.rank, Fraction(args.theta_over_pi))
    largest = max(value for entry in report['numerical'] for key, value in entry.items()
                  if key.endswith('_error'))
    print(f'Exact E6/E7/E8 invariance and all Artin relations: PASS')
    print(f'{len(report["numerical"])} numerical specializations: PASS; max error {largest:.3e}')
    print(f'Wrote {args.output} and {args.export}')

if __name__ == '__main__':
    main()
