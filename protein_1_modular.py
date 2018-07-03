"""
------------------------------------------------------------
 symbolic compution of factorization problem
 Created by: Suksmono@{STEI-ITB, MDR Inc.}
 Applied to protein folding problem
------------------------------------------------------------
 """
from sympy import *
import numpy as np
from prb_symbolic import *
from prb_simulation import *

# define used symbols
NQ=4
qq=symbols('q0:%d'%NQ)

'''
we will calculate the Hamiltonian
H=(Nc- p*q)^2,where p=(1+2*q0+4*q1),q=(1+2*q2);Nc=15, 21 ...
------------------------------------------------------------
    1. USER WRITES THE PROBLEM
------------------------------------------------------------
'''
# Hamiltonian of the problem
Hq = -1 -4*qq[2] + 9*qq[0]*qq[2] + 9*qq[1]*qq[2] \
     - 16*qq[0]*qq[1]*qq[2]
 
'''
------------------------------------------------------------
    2. DO SYMBOLIC COMPUTATION AND TRANSFORM THE PROBLEM 
       FROM Q-DOMAIN TO S-DOMAIN
------------------------------------------------------------
'''

NQ=4   # Number of qubits
# define substition of k-body -> 2-body interaction
qSub=[[1,2,3]] # q1*q2->q3

H2b, H2s = q2s_symbolic(Hq, qSub, NQ)
print('H2b= \n',H2b);
print('H2s= \n',H2s);

'''
------------------------------------------------------------
    3. EXTRACT ISING PARAMETERS FROM SYMBOLIC SOLUTION
------------------------------------------------------------
'''
b, hi, Jij = isingCoeffs(H2s,NQ)

'''
------------------------------------------------------------
    4. OBTAIN A SOLUTION BY SIMULATION (S.A.)
------------------------------------------------------------
'''

SUBITER=100
MAXITER=1000
# calculate the solution using simulated annealing prb_sa
vSpin=prb_sa(b, hi, Jij, NQ, SUBITER, MAXITER, minHq(Hq, NQ))

'''
 ------------------------------------------------------------
     5. DISPLAY THE RESULTS
 ------------------------------------------------------------
 '''
# qubits
qSpin=vs2q(vSpin)
print('\n******************************')
print('            RESULTS')
print('******************************\n')
print('q0=', int(qSpin[0]),',q1=', int(qSpin[1]), \
      ',q2=', int(qSpin[2]), ',q3=', int(qSpin[3]), )

print('\n******************************')
print('  Parameter for D-Wave input')
print('     (to be normalized)')
print('******************************\n')
print('Bias b=',b)
print('hi=',hi)
print('Jij=\n',Jij)

