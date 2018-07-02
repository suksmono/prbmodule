# -*- coding: utf-8 -*-
"""
MODULE: redbool
symbolic manipulation module
"""

import numpy as np
from sympy import *
'''
------------------------------------------------------------
k-body to 2-body transform
------------------------------------------------------------
'''
def H2sub(x1,x2,y,d12):
    return(d12*(3*y+x1*x2-2*x1*y-2*x2*y))

'''
There are two binary representation
    s={-1,+1} and q={0,1}
The default in Ising simulation is s-domain
------------------------------------------------------------
symbolic: q-to-s and s-t-q transform
------------------------------------------------------------
'''
def q2s(x):
    return(1/2 -x/2)

def s2q(x):
    return(1 - 2*x)

'''
------------------------------------------------------------
numerical: q-to-s and s-to-q transforms
------------------------------------------------------------
'''
 # define function vq2s
def vq2s(x):
    return(1-2*x)
# define function vs2q
def vs2q(x):
    return(1/2-x/2)

'''
------------------------------------------------------------
EXTRACT ISING PARAMETERS THE RESULTS FROM SYMBOLIC SOLUTION
copy coefficients into {b, hi, Jij} of Ising parameters
input: Hamiltonian in s-domain H(s0, s1, ...)
output: Ising parameters {b, hi, Jij}
------------------------------------------------------------
'''
def isingCoeffs(Hs, NQ):
    '''
    since we work in s-domain, assume the symbols are
    s0, s1, ...., s_(NQ-1)
    '''
    hi=np.zeros(NQ)
    Jij=np.zeros((NQ,NQ))
    dc=Hs.as_coefficients_dict()
    # list all symbols: s0, s1, ...
    ss=symbols('s0:%d'%NQ)

    # extract b
    b=dc[1]
    # extract hi
    for m in range(NQ):
        hi[m]=dc[ss[m]];
    # extract Jij
    for m in range(NQ):
        for n in range(m+1,NQ):
            print('m=',m,'n=',n)
            Jij[m,n]=dc[ss[m]*ss[n]]
    # return the results
    return(b, hi, Jij)
#
