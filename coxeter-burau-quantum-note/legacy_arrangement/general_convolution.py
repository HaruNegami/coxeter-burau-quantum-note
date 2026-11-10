#!/usr/bin/env python3
"""Arbitrary-rank algebraic convolution and maximal descended quotients.

Input: a finite-dimensional holonomy Lie algebra representation, given by
matrices satisfying all codimension-two relations. Output lives on the
directional closure of the arrangement. No differential equation is solved.
Geometry uses exact rational coefficients; matrix entries may be symbolic.
"""
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
import json
import sympy as s
from arrangement_burau import (convolve, normalize, kohno_check, zero,
                               reduced_burau)


def column_basis(M):
    cols = M.columnspace()
    return s.Matrix.hstack(*cols) if cols else s.zeros(M.rows, 0)


def quotient(matrices, subspace):
    """Return Q,R,induced matrices, with ker(Q)=subspace and Q*R=I."""
    d = subspace.rows
    N = column_basis(subspace)
    T = N
    for i in range(d):
        candidate = T.row_join(s.eye(d)[:, i])
        if candidate.rank() > T.cols:
            T = candidate
    Q = T.inv()[N.cols:, :] if d else s.zeros(0)
    R = T[:, N.cols:]
    Bs = [(Q*A*R).applyfunc(s.simplify) for A in matrices]
    for A, B in zip(matrices, Bs):
        assert zero(Q*A-B*Q), 'Subspace is not invariant'
    assert zero(Q*R-s.eye(Q.rows))
    return Q, R, Bs


@dataclass
class MatrixConvolution:
    variables: tuple
    equations: list
    matrices: list
    transverse_inputs: list
    input_rank: int
    parameter: object
    added: list


def matrix_convolve(equations, input_matrices, variables, y, parameter,
                    check_input=True):
    """Haraoka/Adachi-Hiroe block construction, on the directional closure."""
    if not input_matrices or len(equations) != len(input_matrices):
        raise ValueError('Provide one nonempty square matrix per hyperplane')
    As = [s.Matrix(A) for A in input_matrices]
    r = As[0].rows
    if r == 0 or any(A.shape != (r,r) for A in As):
        raise ValueError('All input matrices must be square of the same rank')
    if check_input:
        kohno_check(equations, As, variables)
    # Reuse only the geometry and its stable hyperplane ordering.
    geom = convolve(equations, [0]*len(equations), variables, y, parameter)
    keys = [normalize(f, variables) for f in geom.equations]
    where = {k:i for i,k in enumerate(keys)}
    transverse, parallel = [], []
    for f,A in zip(equations, As):
        key = normalize(f, variables)
        if s.diff(f,y) == 0:
            parallel.append((key,A))
        else:
            transverse.append((s.expand(y-f/s.diff(f,y)),A,key))
    n = len(transverse)
    lam = s.sympify(parameter)
    Cs = [s.zeros(n*r) for _ in keys]
    for i,(_,Ai,key) in enumerate(transverse):
        C = Cs[where[key]]
        for j,(_,Aj,_) in enumerate(transverse):
            C[i*r:(i+1)*r,j*r:(j+1)*r] = Aj+(lam*s.eye(r) if i==j else s.zeros(r))
    for key,A in parallel:
        Cs[where[key]] = s.Matrix(Cs[where[key]] + s.kronecker_product(s.eye(n),A))
    for i,j in combinations(range(n),2):
        gi,Ai,_ = transverse[i]
        gj,Aj,_ = transverse[j]
        key = normalize(gi-gj, variables)
        if key is None:
            continue
        C = Cs[where[key]]
        C[i*r:(i+1)*r,i*r:(i+1)*r] += Aj
        C[i*r:(i+1)*r,j*r:(j+1)*r] -= Aj
        C[j*r:(j+1)*r,j*r:(j+1)*r] += Ai
        C[j*r:(j+1)*r,i*r:(i+1)*r] -= Ai
    return MatrixConvolution(tuple(variables), geom.equations, Cs,
                             [A for _,A,_ in transverse],r,lam,geom.added)


def matrix_middle(model, substitutions=None):
    """Form E^n/(K+L), specializing first to preserve exceptional ranks."""
    sub = substitutions or {}
    r = model.input_rank
    As = [A.subs(sub) for A in model.transverse_inputs]
    Cs = [C.subs(sub) for C in model.matrices]
    n = len(As)
    vectors = []
    for i,A in enumerate(As):
        for v in A.nullspace():
            lift = s.zeros(n*r,1)
            lift[i*r:(i+1)*r,0] = v
            vectors.append(lift)
    total = s.zeros(n*r)
    for C in Cs[:n]:
        total += C
    vectors += total.nullspace()
    N = column_basis(s.Matrix.hstack(*vectors)) if vectors else s.zeros(n*r,0)
    Q,R,Bs = quotient(Cs,N)
    return Q,R,Bs,N


def maximal_descent(matrices, killed_indices, mode='lie'):
    """Largest quotient on which specified Lie/group generators vanish.

    Lie: kill A_D. Group: kill M_D-I. The seed is saturated under ALL
    generator matrices. Finite-dimensional invariance under invertible M
    also implies invariance under M^-1. Caller must supply every generator
    of the represented algebra/group and the correct normally killed set.
    Generic symbolic ranks are used unless matrices are specialized first.
    """
    if mode not in ('lie','group') or not matrices:
        raise ValueError('Use a nonempty generator list and lie/group mode')
    d = matrices[0].rows
    if any(A.shape != (d,d) for A in matrices):
        raise ValueError('Generator dimensions must agree')
    if mode == 'group':
        for M in matrices:
            if s.simplify(M.det()) == 0:
                raise ValueError('Group generators must be invertible')
    seed = [matrices[i]-(s.eye(d) if mode=='group' else s.zeros(d))
            for i in killed_indices]
    N = column_basis(s.Matrix.hstack(*seed)) if seed else s.zeros(d,0)
    dimensions = [N.cols]
    while N.cols < d:
        enlarged = column_basis(s.Matrix.hstack(N,*[A*N for A in matrices]))
        if enlarged.cols == N.cols:
            break
        N = enlarged
        dimensions.append(N.cols)
    Q,R,Bs = quotient(matrices,N)
    for i in killed_indices:
        target = s.eye(Q.rows) if mode=='group' else s.zeros(Q.rows)
        assert zero(Bs[i]-target)
    return {'Q':Q,'R':R,'matrices':Bs,'killed_submodule':N,
            'saturation_dimensions':dimensions,'quotient_dimension':Q.rows}


def main():
    x,y,u,v=s.symbols('x y u v')
    a,b,c,lam,t=s.symbols('a b c lambda t')
    tri=[y,y-x,y+x-1]
    old=convolve(tri,[a,b,c],(x,y),y,lam)
    new=matrix_convolve(tri,[s.Matrix([[z]]) for z in (a,b,c)],(x,y),y,lam)
    assert old.equations == new.equations
    assert all(zero(A-B) for A,B in zip(old.residues,new.matrices))
    checks={'rank_one_comparison':'PASS'}
    # A genuinely noncommuting rank-two input on a central three-line arrangement.
    A1=s.Matrix([[2,1],[0,0]])/7
    A2=s.Matrix([[0,0],[1,2]])/7
    A3=s.Matrix([[1,-1],[-1,1]])/7
    assert not zero(A1*A2-A2*A1)
    closed=matrix_convolve([y,y-x,x],[A1,A2,A3],(x,y),y,s.Rational(1,11))
    Q,R,Bs,N=matrix_middle(closed)
    assert Q.rows==2
    checks['noncommuting_closed_full']=kohno_check(closed.equations,closed.matrices,(x,y))
    checks['noncommuting_closed_middle']=kohno_check(closed.equations,Bs,(x,y))
    assert maximal_descent(Bs,closed.added)['quotient_dimension']==2
    # Same arrangement, a direction transverse to all original lines.
    changed=matrix_convolve([u+v,2*v,u-v],[A3,A1,A2],(u,v),v,s.Rational(1,11))
    Q,R,Bs,N=matrix_middle(changed)
    assert Q.rows==3 and len(changed.added)==1
    checks['noncommuting_changed_full']=kohno_check(changed.equations,changed.matrices,(u,v))
    checks['noncommuting_changed_middle']=kohno_check(changed.equations,Bs,(u,v))
    d_changed=maximal_descent(Bs,changed.added)
    # Equal-weight Burau-type specialization: added residue saturation kills all.
    alpha=s.Rational(1,7)
    equal=matrix_convolve(tri,[s.Matrix([[alpha]])]*3,(x,y),y,alpha)
    d_equal=maximal_descent(equal.matrices,equal.added)
    assert d_equal['saturation_dimensions']==[2,3]
    assert d_equal['quotient_dimension']==0
    # A single added divisor already suffices in this equal-weight case.
    assert maximal_descent(equal.matrices,equal.added[:1])['quotient_dimension']==0
    # At lambda=0 the rank-one input survives and all added residues vanish.
    identity=matrix_convolve(tri,[s.Matrix([[alpha]])]*3,(x,y),y,0)
    Q,R,Bs,N=matrix_middle(identity)
    assert Q.rows==1
    assert maximal_descent(Bs,identity.added)['quotient_dimension']==1
    # Group variant: actual pure-braid matrices, not exponentiated residues.
    s1,s2=(reduced_burau(3,i,t) for i in [1,2])
    pure=[s1**2,s2**2,s2*s1**2*s2.inv()]
    group_generic=maximal_descent(pure,[0],mode='group')
    assert group_generic['quotient_dimension']==0
    group_at_one=maximal_descent([M.subs(t,1) for M in pure],[0],mode='group')
    assert group_at_one['quotient_dimension']==2
    output={
        'checks':checks,
        'changed_direction_middle_rank':3,
        'changed_direction_descended_rank':d_changed['quotient_dimension'],
        'equal_weight_triangle':{
            'rank':3,'saturation_dimensions':d_equal['saturation_dimensions'],
            'descended_rank':d_equal['quotient_dimension']},
        'lambda_zero_middle_descended_rank':1,
        'group_P3_kill_A12_generic_rank':group_generic['quotient_dimension'],
        'group_P3_kill_A12_at_q1_rank':group_at_one['quotient_dimension'],
        'status':'ALL CHECKS PASS',
    }
    Path(__file__).with_name('general_results.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__':
    main()
