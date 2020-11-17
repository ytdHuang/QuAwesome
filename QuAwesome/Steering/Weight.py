# This code is part of QuAwesome.
#
#    MIT License
#
#    Copyright (c) 2020 and later, Yi-Te Huang
#
#    Permission is hereby granted, free of charge, to any person obtaining a copy
#    of this software and associated documentation files (the "Software"), to deal
#    in the Software without restriction, including without limitation the rights
#    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#    copies of the Software, and to permit persons to whom the Software is
#    furnished to do so, subject to the following conditions:
#
#    The above copyright notice and this permission notice shall be included in all
#    copies or substantial portions of the Software.
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#    SOFTWARE.
#######################################################################################
import sys
import numpy as np
import qutip as qu 
import cvxopt as cvx
from picos import Problem, sum, trace, value
from picos.expressions.variables import HermitianVariable
from QuAwesome import QuAwesomeError as ERROR
from QuAwesome.Steering.genDeterministic import genDeterministicArray as genD

def steeringWeight(assemb, solver, returnF=False) :
    """
    Calculate Steering Weight \n
    Inputs:
        assemb - a 4-D array containing the assemblage members 
        Each element is an unnormalized N * N density matrix: sigma_a|X where 'a' denote outcome condition on measurement 'X'
            example : if there are totally M+1 measurements, and each has A+1 outcomes
                assemblage = [
                    [ sigma_0|0 , sigma_1|0 , ... , sigma_A|0 ],
                    [ sigma_0|1 , sigma_1|1 , ... , sigma_A|1 ],
                        .
                        .
                        .
                    [ sigma_0|M , sigma_1|M , ... , sigma_A|M ]
                ]
        
        solver  - a string of solver (mosek or cvxopt)
        returnF - choose to return F_a|X or not [Default as False]
    """

    # Get dimension info. of assemblage and check if it is valid
    try:
        (M, A, N, N1) = np.shape(assemb)
    
    except ValueError:
        raise ERROR("The dimension of input assemblage is incorrect")

    if N != N1 :
        raise ERROR("The dimension of unnormalized density matrix is incorrect.\nIt should be N x N, with N >= 2")

    # number of local hidden variable
    Nlamb = A ** M

    # start to solve steerable weight by Semidefinite program
    SP = Problem()

    # add variable (F_a|X)
    F = []
    for x in range(M):
        F.append( [ HermitianVariable('F_{0}|{1}'.format(a, x), (N, N)) for a in range(A) ] )

    # generate Deterministic probability distribution array
    D = genD(M, A)

    # add constraints
    SP.add_list_of_constraints([F[x][a] >> 0  for x in range(M) for a in range(A)])
    for l in range(Nlamb):

        # sum_aX D(a|X, lambda) * F_a|X
        summation = sum( [ D[l][x][a] * F[x][a] for x in range(M) for a in range(A) ] )
        
        # iden - sum_aX D(a|X, lambda) * F_a|X <= 0
        SP.add_constraint(
            np.eye(N, dtype=int) - summation << 0
        )

    # sum_aX F_a|X * Sigma_a|X
    summation = sum( [ F[x][a] * assemb[x][a] for x in range(M) for a in range(A) ] )

    # find the solution
    SP.set_objective(
        'max',
        (1 - np.real(trace(summation)))
    )

    # solve the problem
    SP.solve(solver=solver)

    # return results
    if returnF == True:
        F_ax = []
        for x in range(M):
            F_ax.append( [ value(F[x][a], numpy=True) for a in range(A) ] )

        return SP.value, F_ax

    else:
        return SP.value