"""
 ------------------------------------------------------------
 symbolic compution of factorization problem
 Created by: Suksmono@{STEI-ITB, MDR Inc.}
 problem adopted from:
 https://qiita.com/YuichiroMinato/items/9c9fba2b5bc978ec9e73
 the idea: 
 (1) a user enter his/her problem
 (2) server do symbolic calculation and simulation: 
     (a)P0 inc k-body -> P1 with only 2-bdy in s-domain
     (b) extract coefficients {b, hi, Jij}
           -> fed to simulator
           -> deliver to user for D-Wave input
 (3) server send back simulation results and parameter to user
 (4) user run his problem (b, hi, Jij) in D-wave, compare the
     result from simulator
 ------------------------------------------------------------
 """
from sympy import *
import numpy as np
from prb_symbolic import *
from prb_simulation import*

# define used symbols
q0, q1, q2, q3, d = symbols('q0 q1 q2 q3 d')

'''
we will calculate the Hamiltonian
H=(15 - p*q)^2, where p=(1+2*q0 + 4*q1), q=(1+2*q2)
------------------------------------------------------------
 USER WRITES THE PROBLEM
it is known that the prime factors are odd numbers
------------------------------------------------------------
'''
Nc=21;  # a composite number to be factored
p=1+2*q0+4*q1;
q=1+2*q2;
# Hamiltonian of the problem
H=(Nc-p*q)**2;

'''
------------------------------------------------------------
 SYMBOLIC COMPUTATION
 input: H in q-domain
 output: Hs in s-domain
------------------------------------------------------------
'''
#def Hq2s(H):
EQ0=expand(H);
# display result for checking
print('EQ0=',EQ0);

# substitute qi^2->qi
EQ1=simplify( \
        EQ0.subs({q0**2:q0, q1**2:q1, q2**2:q2}) \
    );

print('EQ1=\n',EQ1);

H2b=simplify( \
          EQ1.subs({q1*q2:q3}) \
        + H2sub(q1,q2,q3,d) \
    );
print('H2b= \n',H2b);

'''
------------------------------------------------------------
transform q={0,1} to s{-1,1}
------------------------------------------------------------
'''
s0, s1, s2, s3 = symbols('s0 s1 s2 s3')

'''
------------------------------------------------------------
considering d>max(H), letting q0=q1=q2=0, 
H=(Nc-1)^2 => d=Nc^2
------------------------------------------------------------
'''
H2b1=simplify( \
          H2b.subs({d:Nc*Nc}) \
    );
print('H2b1= \n',H2b1);

'''
------------------------------------------------------------
TRANSFORM H in q-domain TO s-DOMAIN FOR SIMULATION
------------------------------------------------------------
'''
H2s=simplify( H2b1.subs( \
               {q0:q2s(s0), q1:q2s(s1), q2:q2s(s2), q3:q2s(s3), }
              ) \
     );
print('H2s= \n',H2s);

'''
------------------------------------------------------------
EXTRACT ISING PARAMETERS THE RESULTS FROM SYMBOLIC SOLUTION
------------------------------------------------------------
'''
NQ=4 # number of qubits q0, ..., q3
b, hi, Jij = isingCoeffs(H2s,NQ)

'''
------------------------------------------------------------
SEND THE RESULTS (Jij, hi, b) TO SIMULATOR (S.A.)
------------------------------------------------------------
'''

SUBITER=100
MAXITER=1000
# calculate the solution using simulated annealing prb_sa
vSpin=prb_sa(b, hi, Jij, NQ, SUBITER, MAXITER)

'''
 ------------------------------------------------------------
 DISPLAY THE RESULTS
 ------------------------------------------------------------
 '''
# qubits
qSpin=vs2q(vSpin)
print('\nRESULTS')
print('q0=', int(qSpin[0]),',q1=', int(qSpin[1]), \
      ',q2=', int(qSpin[2]), ',q3=', int(qSpin[3]), )
#

vp=int(p.subs({q0:qSpin[0],q1:qSpin[1]}))
vq=int(q.subs({q2:qSpin[2]}))
print('Prime factors of ', Nc, ' are', vp, 'and', vq)
