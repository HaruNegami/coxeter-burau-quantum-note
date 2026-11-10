#!/usr/bin/env python3
"""Based monodromy matrices for the equal-weight Burau-type construction.

All matrix computations are exact over Q(q). Path tracking is an independent
floating-point check of the explicitly described crossing words, not a
numerical integration of the connection and not an interval certificate.
Group products act on columns: the right-hand path is traversed first.
"""
from itertools import combinations
from pathlib import Path
import cmath
import json
import sympy as s
from arrangement_burau import reduced_burau, burau_word, convolve, zero


def pure_word(i,j):
    if not 1<=i<j:
        raise ValueError('Require 1<=i<j')
    return list(range(j-1,i,-1))+[i,i]+[-k for k in range(i+1,j)]


def pure_factors(strands,i,j,q):
    """beta(A_ij)=I+(q-1)*u*v, quotient basis w_k-w_strands."""
    if not 1<=i<j<=strands:
        raise ValueError('Require 1<=i<j<=strands')
    n=strands-1
    e=s.eye(n)
    if j<=n:
        u=q*e[:,i-1]-e[:,j-1]
        for k in range(i,j-1):
            u+=(q-1)*e[:,k]
        v=e[i-1,:]-e[j-1,:]
    else:
        u=s.zeros(n,1)
        for k in range(n):
            u[k]=1 if k<i-1 else q+1 if k==i-1 else q
        v=e[i-1,:]
    return u,v


def pure_matrix(strands,i,j,q):
    u,v=pure_factors(strands,i,j,q)
    return (s.eye(strands-1)+(q-1)*u*v).applyfunc(s.expand)


def line(a,b,count=97):
    return [a+(b-a)*k/count for k in range(count+1)]


def arc(c,r,start,finish,count=397):
    return [c+r*cmath.exp(1j*(start+(finish-start)*k/count))
            for k in range(count+1)]


def join(*paths):
    answer=[]
    for path in paths:
        answer+=path if not answer else path[1:]
    return answer


def based_loop(stem,local_circle):
    return join(stem,local_circle,list(reversed(stem)))


def crossing_word(samples):
    """Chronological positive CCW half twists, then reverse for group product."""
    chronological=[]
    before=samples[0]
    order=sorted(range(len(before)),key=lambda i:before[i].real)
    for after in samples[1:]:
        next_order=sorted(range(len(after)),key=lambda i:after[i].real)
        changed=[i for i in range(len(order)) if order[i]!=next_order[i]]
        if changed:
            assert len(changed)==2 and changed[1]==changed[0]+1, 'Refine path sampling'
            k=changed[0]
            left,right=order[k],order[k+1]
            imagdiff=((before[left]-before[right])+(after[left]-after[right])).imag/2
            assert abs(imagdiff)>1e-5, 'Projection crossing is too close to collision'
            chronological.append(k+1 if imagdiff<0 else -(k+1))
        before,order=after,next_order
    assert order==list(range(len(order))), 'Loop is not pure at this basepoint'
    return list(reversed(chronological))


def geometric_triangle_words():
    """Base (x,y)=(1/4,2), labelled configuration (0,x,1-x,y)."""
    pi=cmath.pi
    r=.05
    words={}
    # Moving y has stems below the real axis, then positive local circles.
    for label,target in [('y=0',0),('y=x',.25),('y=1-x',.75)]:
        stem=join(line(2,2-.1j),line(2-.1j,target+r-.1j),line(target+r-.1j,target+r))
        ys=based_loop(stem,arc(target,r,0,2*pi))
        words[label]=crossing_word([(0,.25,.75,y) for y in ys])
    # x=0 and x=1/2 have real stems. The stem to x=1 bypasses 1/2 above.
    stem0=line(.25,r)
    stemhalf=line(.25,.5-r)
    stem1=join(line(.25,.4),arc(.5,.1,pi,0),line(.6,1-r))
    for label,stem,center,start in [('x=0',stem0,0,0),
                                   ('x=1/2',stemhalf,.5,pi),
                                   ('x=1',stem1,1,pi)]:
        xs=based_loop(stem,arc(center,r,start,start+2*pi))
        words[label]=crossing_word([(0,x,1-x,2) for x in xs])
    return words


def main():
    q=s.symbols('q')
    checks=0
    for m in (3,4,5,6):
        for i,j in combinations(range(1,m+1),2):
            M=pure_matrix(m,i,j,q)
            assert zero(M-burau_word(m,pure_word(i,j),q))
            assert zero((M-s.eye(m-1))*(M-q*q*s.eye(m-1)))
            assert s.factor(M.det())==q*q
            checks+=1
    labels={'y=0':(1,4),'y=x':(2,4),'y=1-x':(3,4),
            'x=0':(1,2),'x=1':(1,3),'x=1/2':(2,3)}
    path_words=geometric_triangle_words()
    for label,(i,j) in labels.items():
        assert path_words[label]==pure_word(i,j),(label,path_words[label])
    # First-order comparison with every residue in the original convention.
    x,y=s.symbols('x y')
    model=convolve([y,y-x,y+x-1],[1,1,1],(x,y),y,1)
    S=s.eye(3)+s.ones(3)
    ordered_pairs=[(1,4),(2,4),(3,4),(1,2),(1,3),(2,3)]
    for C,(i,j) in zip(model.residues,ordered_pairs):
        infinitesimal=pure_matrix(4,i,j,q).diff(q).subs(q,1)
        assert zero(infinitesimal-S*C*S.inv())
    full=s.eye(3)
    for j in range(2,5):
        for i in range(1,j):
            full=full*pure_matrix(4,i,j,q)
    assert zero(full-q**4*s.eye(3))
    records={}
    for label,(i,j) in labels.items():
        M=pure_matrix(4,i,j,q)
        u,v=pure_factors(4,i,j,q)
        records[label]={'pure_generator':f'A_{i}{j}', 'word':pure_word(i,j),
                        'matrix':[[str(z) for z in M.row(k)] for k in range(3)],
                        'latex':s.latex(M), 'u':[str(z) for z in u],
                        'v':[str(z) for z in v]}
    out={'parameter':'q = exp(2*pi*i*alpha), a_i=lambda=alpha, parallel input weights zero',
         'basepoint':['1/4','2'],'labelled_configuration':['0','x','1-x','y'],
         'convention':'columns; right-hand path first; positive CCW local circles',
         'monodromy_matrices':records,
         'checks':{'pure_generators_verified':checks,
                   'geometric_meridian_words_verified':6,
                   'residue_first_derivatives_verified':6,
                   'full_twist':'q^4 I_3','status':'PASS'},
         'scope':'These are the based matrices on the enlarged six-line complement. Generic MC/Burau identification; Burau matrices themselves are defined for q != 0.'}
    Path(__file__).with_name('monodromy_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
