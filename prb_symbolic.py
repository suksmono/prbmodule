"""
MODULE: redbool
symbolic manipulation module
"""

import numpy as np
from sympy import *

'''
------------------------------------------------------------
maximum value of Hq, used to estimate delta > max(Hq)
------------------------------------------------------------
'''
def maxHq(Hq, NQ):
    qq=symbols('q0:%d'%NQ)
    vE=zeros(2**NQ,1)
    for m in range(0,2**NQ):
        vq=np.binary_repr(m, width=NQ)
        tE=Hq
        for n in range(0,NQ):
            tE=tE.subs({qq[n]:float(vq[n]) })
        # put energy spectra to vE
        vE[m]= float(tE)
    return(max(vE))
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
            #print('m=',m,'n=',n)
            Jij[m,n]=dc[ss[m]*ss[n]]
    # return the results
    return(b, hi, Jij)
#
'''
------------------------------------------------------------
    SYMBOLIC COMPUTATION
------------------------------------------------------------
  >input:  Hq in q-domain
  >output: Hs in s-domain
 ALGORITHM
  1> Expand H
  2> Simplify-1: substitute  qi**2 <- qi
  3> Simplify-2: transform k-body terms-> 2-body+compenste
  4> Transform Hq-> Hs
Process:
    Hq0 -> Hq -> Hq2b -> H2s
 Hq0: original input
 Hq1: expanded Hq0
 Hq2b: Hq with only 2-body interaction
 H2s: H in s-domain with 2-body interaction
------------------------------------------------------------
'''
#def q2s_symbolic(H,NQ,d):
def q2s_symbolic(H,qSub,NQ):
    #def Hq2s(H):
    Hq0=expand(H);
    # display result for checking
    # print('Hq0=\n',Hq0);
    # identify all parameters based on NQ
    qq=symbols('q0:%d'%NQ)
    # substitute qi^2->qi
    Hq=Hq0;
    # substiture qi**2 <- qi recursively
    for m in range(NQ):
        Hq=simplify( \
            Hq.subs({qq[m]**2:qq[m]}) \
        )
    
    #print('Hq=\n',Hq);
    # extract terms of EQ1
    dc=Hq.as_coefficients_dict()
    
    """
    put attention to q0*q1*...
     q0*q1 -> q_{NQ}
     q1*q2 -> q_{NQ+1}
     ...
     q_{NQ-2}*q_{NQ-1} -> q_{floor(NQ/2)}
    ----
     1. detect and get high order term (qi*qj*qk...)
     2. substitute (qi*qj->qm) and compensate
     3. goto 1
    ....
    aq=tq.atoms(Symbol)
    len(aq)
    ====
    check all combination C(NQ,NQ),C(NQ,NQ-1), .., C(NQ,3) 
    """
    d=2*maxHq(Hq,NQ)
       
    # do substitution iteratively 
    H2b=Hq
    for m in range(0,len(qSub)):
        H2b=simplify( \
                  H2b.subs({qq[qSub[m][0]]*qq[qSub[m][1]]:qq[qSub[m][2]]}) \
                + H2sub(qq[qSub[m][0]],qq[qSub[m][1]],qq[qSub[m][2]],d) \
            );

    #print('H2b= \n',H2b);
    
    '''
    ------------------------------------------------------------
    transform q={0,1} to s{-1,1}
    ------------------------------------------------------------
    '''
    #s0, s1, s2, s3 = symbols('s0 s1 s2 s3')
    # define s-domain variables
    ss=symbols('s0:%d'%NQ)
    '''
    ------------------------------------------------------------
    TRANSFORM H in q-domain TO s-DOMAIN FOR SIMULATION
    ------------------------------------------------------------
    '''
    H2s=simplify( H2b.subs( \
                   {qq[0]:q2s(ss[0]), qq[1]:q2s(ss[1]), qq[2]:q2s(ss[2]),\
                    qq[3]:q2s(ss[3]), }
                  ) \
         );
    #print('H2s= \n',H2s);
    return(H2b, H2s)
