#!/usr/bin/env python3
"""Rank-one middle convolution, Negami's form, quotient and unitary gates.

Angles are turns, not additive residues. Fiber loops follow Negami's
algebraic convention (clockwise); reverse geometric loops by inverting.
Braid words are chronological: the next matrix multiplies on the left.
Only NumPy is required; --verify-symbolic additionally uses SymPy.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import numpy as np


def signature(angles, ell):
    """Exact for Fraction inputs. Returns positive, negative, zero counts."""
    a = [x % 1 for x in angles]
    ell = ell % 1
    if ell == 0:
        raise ValueError('This formula assumes lambda != 1.')
    total = ell + sum(a)
    if isinstance(total, Fraction):
        resonance = int(total.denominator == 1)
        m = total.numerator // total.denominator
    else:
        nearest = round(float(total))
        resonance = int(abs(float(total) - nearest) < 1e-10)
        m = nearest if resonance else math.floor(float(total))
    p = sum(x != 0 for x in a)
    return p-m, m-resonance, len(a)-p+resonance


def ambient(q, lam, root):
    n = len(q)
    H = np.array([[(1/root if i <= j else root) *
                   (1/q[i] - (lam if i == j else 1)) * (q[j]-1)
                   for j in range(n)] for i in range(n)], complex)
    G = []
    for i in range(n):
        g = np.eye(n, dtype=complex)
        g[i] = [lam*(q[j]-1) if j < i else
                lam*q[i] if j == i else q[j]-1 for j in range(n)]
        G.append(g)
    return H, G


def adjacent_transport(q, i):
    r = np.eye(len(q), dtype=complex)
    r[i:i+2, i:i+2] = [[0, q[i]], [1, 1-q[i+1]]]
    return r


def colored_braid(q, word):
    """Return transport and final colors; pure words return original colors."""
    colors = np.array(q, dtype=complex).copy()
    B = np.eye(len(q), dtype=complex)
    for letter in word:
        i = abs(letter)-1
        if not 0 <= i < len(q)-1 or letter == 0:
            raise ValueError('Braid letter outside the available strands.')
        swapped = colors.copy()
        swapped[i], swapped[i+1] = colors[i+1], colors[i]
        R = adjacent_transport(colors, i) if letter > 0 else np.linalg.inv(adjacent_transport(swapped, i))
        B = R @ B
        colors = swapped
    return B, colors


def pure_word(i, j):
    """One-based pair; chronological stem, double crossing, reverse stem."""
    if not 1 <= i < j:
        raise ValueError('Require 1 <= i < j.')
    stem = list(range(j-1, i, -1))
    return stem + [i, i] + [-s for s in reversed(stem)]


def build(angles, ell, base_words=None, base_phases=None):
    a = np.array([float(x % 1) for x in angles])
    ell = float(ell % 1)
    if abs(ell) < 1e-14:
        raise ValueError('lambda = 1 is outside this implementation.')
    q = np.exp(2j*np.pi*a)
    lam, root = np.exp(2j*np.pi*ell), np.exp(1j*np.pi*ell)
    H, fiber = ambient(q, lam, root)
    matrices = {f'fiber_{i+1}': g for i, g in enumerate(fiber)}
    if base_words is None:
        base_words = {f'base_A{i}_{j}': pure_word(i, j)
                      for i in range(1, len(a)) for j in range(i+1, len(a)+1)}
    base_phases = base_phases or {}
    for name, word in base_words.items():
        permutation=list(range(len(q)))
        for letter in word:
            i=abs(letter)-1
            if letter==0 or not 0<=i<len(q)-1:
                raise ValueError('Invalid base braid letter.')
            permutation[i],permutation[i+1]=permutation[i+1],permutation[i]
        if permutation!=list(range(len(q))):
            raise ValueError(f'{name} is not a pure root braid.')
        B, end = colored_braid(q, word)
        if not np.allclose(end, q, atol=1e-10):
            raise ValueError(f'{name} does not preserve the seed colors.')
        matrices[name] = np.exp(2j*np.pi*float(base_phases.get(name, 0))) * B
    vals, vecs = np.linalg.eigh(H)
    tol = 1e-10 * max(1, np.max(abs(vals)))
    Q, N = vecs[:, abs(vals)>tol], vecs[:, abs(vals)<=tol]
    h = Q.conj().T @ H @ Q
    quotient = {name: Q.conj().T @ B @ Q for name, B in matrices.items()}
    d = Q.shape[1]
    sign = 1 if d and np.min(np.linalg.eigvalsh(h)) > tol else -1 if d and np.max(np.linalg.eigvalsh(h)) < -tol else 0
    gates = {}
    if sign:
        # np.linalg.cholesky gives L L*, whereas we need S* S = sign*h.
        S = np.linalg.cholesky(sign*h).conj().T
        invS = np.linalg.inv(S)
        gates = {name: S @ B @ invS for name, B in quotient.items()}
    err = lambda A: float(np.max(abs(A))) if A.size else 0.0
    errors = {
        'hermitian': err(H-H.conj().T),
        'ambient_invariance': max(err(B.conj().T@H@B-H) for B in matrices.values()),
        'kernel_invariance': max(err(Q.conj().T@B@N) for B in matrices.values()),
        'quotient_invariance': max(err(B.conj().T@h@B-h) for B in quotient.values()),
        'unitarity': max((err(U.conj().T@U-np.eye(d)) for U in gates.values()), default=None)
    }
    return {'q':q, 'lambda':lam, 'H':H, 'Q':Q, 'H_quotient':h,
            'ambient':matrices, 'quotient':quotient, 'gates':gates,
            'dimension':d, 'sign_used':sign, 'eigenvalues':vals, 'errors':errors}


def rational_finite_test(angles, ell):
    """Exact rank-one criterion; assumes the parallel/base seed has finite image.

    Because the arrangement group contains its punctured-line fiber group,
    definiteness at every cyclotomic embedding is necessary and sufficient.
    Dimension zero is reported separately. No group enumeration is used.
    """
    a, ell = [Fraction(x) % 1 for x in angles], Fraction(ell) % 1
    if not ell:
        raise ValueError('lambda = 1 is excluded.')
    if sum(signature(a, ell)[:2]) == 0:
        return {'status':'zero_output_no_gate_space',
                'finite_if_base_character_finite':True,
                'reason':'The zero-dimensional representation is trivial.'}
    N = math.lcm(*(x.denominator for x in a), ell.denominator)
    if N > 100000:
        raise ValueError('Modulus too large for the default exhaustive check.')
    rows = []
    for k in range(1, N):
        if math.gcd(k,N) == 1:
            sig = signature([k*x % 1 for x in a], k*ell % 1)
            rows.append({'k':k, 'signature':sig, 'definite':sum(sig[:2])>0 and min(sig[:2])==0})
    return {'modulus':N, 'finite_if_base_character_finite':all(x['definite'] for x in rows),
            'galois_signatures':rows}


def weyl_clifford_residual(U):
    """Numerical normalizer check in one fixed d-level Weyl frame, d <= 8."""
    d = len(U)
    if d > 8 or d < 2:
        return None
    X = np.roll(np.eye(d), 1, axis=0)
    Z = np.diag(np.exp(2j*np.pi*np.arange(d)/d))
    basis = [np.linalg.matrix_power(X,a)@np.linalg.matrix_power(Z,b)
             for a in range(d) for b in range(d)]
    distances=[]
    for P in (X,Z):
        V = U@P@U.conj().T
        best=min(np.linalg.norm(V-(np.trace(W.conj().T@V)/d)/abs(np.trace(W.conj().T@V)/d)*W)
                 if abs(np.trace(W.conj().T@V)/d)>1e-12 else math.inf for W in basis)
        distances.append(float(best))
    return max(distances)


def verify_symbolic():
    import sympy as sp
    q=sp.symbols('q0:3', nonzero=True); s=sp.symbols('s',nonzero=True); lam=s*s
    inv={x:1/x for x in (*q,s)}
    star=lambda A:A.T.applyfunc(lambda x:x.subs(inv, simultaneous=True))
    simp=lambda A:all(sp.cancel(x)==0 for x in A)
    H=sp.Matrix(3,3,lambda i,j:(1/s if i<=j else s)*(1/q[i]-(lam if i==j else 1))*(q[j]-1))
    assert simp(star(H)-H)
    for i in range(3):
        G=sp.eye(3)
        for j in range(3):G[i,j]=lam*(q[j]-1) if j<i else lam*q[i] if j==i else q[j]-1
        assert simp(star(G)*H*G-H)
    for i in range(2):
        R=sp.eye(3);R[i:i+2,i:i+2]=sp.Matrix([[0,q[i]],[1,1-q[i+1]]])
        Hswap=H.subs({q[i]:q[i+1],q[i+1]:q[i]},simultaneous=True)
        assert simp(star(R)*Hswap*R-H)
        def generators(colors):
            out=[]
            for k in range(3):
                G=sp.eye(3)
                for j in range(3):G[k,j]=lam*(colors[j]-1) if j<k else lam*colors[k] if j==k else colors[j]-1
                out.append(G)
            return out
        qs=list(q);qs[i],qs[i+1]=qs[i+1],qs[i]
        G,Gp=generators(q),generators(qs)
        for j in range(3):
            rhs=Gp[i+1] if j==i else Gp[i+1].inv()*Gp[i]*Gp[i+1] if j==i+1 else Gp[j]
            assert simp(R*G[j]-rhs*R)
    return 'PASS: formal Laurent-rational fiber, colored-form and semidirect-action identities (n=3)'


def verify():
    rng=np.random.default_rng(20261110)
    worst=0.; checks=0
    for n in range(1,9):
        for _ in range(20):
            a=rng.uniform(.02,.98,n);ell=rng.uniform(.02,.98)
            result=build(a,ell)
            expected=signature(a,ell)
            ev=result['eigenvalues']; got=(int(sum(ev>1e-8)),int(sum(ev< -1e-8)),int(sum(abs(ev)<=1e-8)))
            assert got==expected,(got,expected)
            worst=max(worst,*[x for x in result['errors'].values() if x is not None])
            q=result['q']
            if n>=3:
                B,c=colored_braid(q,[1,2,1]);C,e=colored_braid(q,[2,1,2])
                assert np.max(abs(B-C))<1e-10 and np.max(abs(c-e))<1e-10
            checks+=1
    for a,ell in [([Fraction(1,4)]*3,Fraction(1,4)),
                  ([Fraction(0),Fraction(1,4),Fraction(1,4)],Fraction(1,2)),
                  ([Fraction(1,6)]*3,Fraction(1,3)),
                  ([Fraction(3,4)]*3,Fraction(3,4))]:
        r=build(a,ell);sig=signature(a,ell)
        assert r['dimension']==sum(sig[:2])
        assert r['sign_used'] != 0
        assert max(x for x in r['errors'].values() if x is not None)<1e-10
        checks+=1
    assert rational_finite_test([Fraction(1,6)]*2,Fraction(1,6))['finite_if_base_character_finite']
    assert rational_finite_test([Fraction(1,4)]*3,Fraction(1,4))['finite_if_base_character_finite']
    assert not rational_finite_test([Fraction(1,10)]*3,Fraction(1,10))['finite_if_base_character_finite']
    assert worst<1e-9,worst
    return {'status':'PASS','tuple_cases':checks,'max_matrix_error':worst,
            'rational_finiteness_examples':3,
            'scope':'Numerical tuple checks supplement the proofs; the finite-image test uses exact Fractions.'}


def encode(x):
    if isinstance(x,np.ndarray):
        return {'real':x.real.tolist(),'imag':x.imag.tolist()} if np.iscomplexobj(x) else x.tolist()
    if isinstance(x,(complex,np.complexfloating)):return [float(x.real),float(x.imag)]
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [encode(v) for v in x]
    return x


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--angles',default='1/20,1/15,1/12')
    parser.add_argument('--ell',default='sqrt2/20',help='rational turns or sqrt2/20')
    parser.add_argument('--output',default='mc_gate_example.json')
    parser.add_argument('--base-words',help='JSON with words and optional rational phases for actual based root braids')
    parser.add_argument('--verify',action='store_true')
    parser.add_argument('--verify-symbolic',action='store_true')
    args=parser.parse_args()
    a=[Fraction(x) for x in args.angles.split(',')]
    ell=math.sqrt(2)/20 if args.ell=='sqrt2/20' else Fraction(args.ell)
    base_data=json.loads(Path(args.base_words).read_text()) if args.base_words else {'words':{}}
    r=build(a,ell,base_words=base_data['words'],base_phases={k:Fraction(v) for k,v in base_data.get('phases',{}).items()})
    r['input']={'angles':[str(x) for x in a], 'ell':args.ell,
                'angle_unit':'turns',
                'base_character':'supplied rational phases' if base_data.get('phases') else 'trivial on supplied base loops',
                'base_words':base_data['words'],
                'base_phases':base_data.get('phases',{}),
                'loop_convention':'Negami clockwise fiber loops; chronological braid words',
                'base_scope':'Only supplied root braids are exported; by default only the actual fiber generators are exported.'}
    r['signature']=signature(a,ell)
    r['weyl_clifford_residuals']={k:weyl_clifford_residual(U) for k,U in r['gates'].items()}
    if isinstance(ell,Fraction):r['finite_image']=rational_finite_test(a,ell)
    if args.verify:r['verification']=verify()
    if args.verify_symbolic:r['symbolic_verification']=verify_symbolic()
    Path(args.output).write_text(json.dumps(encode(r),indent=2)+'\n')
    print('MC quotient dimension:',r['dimension'],'; signature:',r['signature'])
    print('Invariance / quotient / unitarity errors:',r['errors'])
    if args.verify:print(r['verification'])
    if args.verify_symbolic:print(r['symbolic_verification'])
    print('Wrote',args.output)


if __name__=='__main__':main()
