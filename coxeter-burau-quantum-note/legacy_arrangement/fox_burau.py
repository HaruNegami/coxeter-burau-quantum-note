#!/usr/bin/env python3
"""Burau from Fox derivatives: a direct group-algebra construction.

Supply the free-group automorphism of a based loop, e.g. obtained from its
braid monodromy. This does NOT determine braid monodromy from arrangement
equations. Words use signed 1-based generator indices.
"""
import sympy as s
from arrangement_burau import reduced_burau, zero


def fox_column(word, rank, q):
    """Evaluate all Fox derivatives under x_1=...=x_rank=q."""
    result=s.zeros(rank,1)
    exponent=0
    for letter in word:
        i=abs(letter)-1
        if not 0<=i<rank:
            raise ValueError('Free-group generator outside the selected rank')
        if letter>0:
            result[i]+=q**exponent
            exponent+=1
        else:
            result[i]-=q**(exponent-1)
            exponent-=1
    return result.applyfunc(s.expand)


def fox_homology(images, q):
    """Action on H_1(F_m,K_q)=ker((q-1,...,q-1)), over K=Q(q).

    images[j] is the word for phi(x_(j+1)). The caller supplies an
    automorphism. Character preservation and invertibility are checked,
    but invertibility of the evaluated Jacobian alone does not prove that
    arbitrary supplied words define a free-group automorphism.
    Basis: e_i-e_m, i=1,...,m-1. q is generic, not specialized to 1.
    """
    m=len(images)
    if m<2 or q==1:
        raise ValueError('Require m>=2 and generic q; q=1 changes twisted H_1 rank')
    for word in images:
        if sum(1 if i>0 else -1 for i in word)!=1:
            raise ValueError('Automorphism must preserve the total-exponent character')
    J=s.Matrix.hstack(*[fox_column(w,m,q) for w in images])
    assert zero(s.ones(1,m)*J-s.ones(1,m))
    if s.simplify(J.det())==0:
        raise ValueError('Evaluated Jacobian is not invertible')
    R=s.eye(m-1).col_join(-s.ones(1,m-1))
    B=(J*R)[:m-1,:].applyfunc(s.factor)
    assert zero(J*R-R*B)
    return B


def artin_generator(rank,index):
    """sigma_i: x_i -> x_i*x_(i+1)*x_i^-1; x_(i+1) -> x_i."""
    if not 1<=index<rank:
        raise ValueError('Require 1<=index<rank')
    images=[[j] for j in range(1,rank+1)]
    images[index-1]=[index,index+1,-index]
    images[index]=[index]
    return images


def main():
    q=s.symbols('q')
    for m in (3,4,5):
        # Generic change from the quotient model in arrangement_burau.py
        # to this twisted-homology model. Denominator sum q^j is nonzero
        # in Q(q); specializations require separate treatment.
        D=s.diag(*[q**j for j in range(m)])
        ell=s.Matrix([[q**j for j in range(m)]])
        P=s.eye(m)-s.ones(m,1)*ell/sum(q**j for j in range(m))
        R=s.eye(m)[:,:m-1]
        T=(D*P*R)[:m-1,:]
        for i in range(1,m):
            B=fox_homology(artin_generator(m,i),q)
            U=reduced_burau(m,i,q)
            assert zero(B*T-T*U)
    print('Fox construction matches reduced Burau over Q(q), B3/B4/B5: PASS')
    print('B3 sigma_1 =',fox_homology(artin_generator(3,1),q))
    print('B3 sigma_2 =',fox_homology(artin_generator(3,2),q))


if __name__=='__main__':
    main()
