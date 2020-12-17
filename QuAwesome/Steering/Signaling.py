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
from QuAwesome import QuAwesomeError as ERROR

def Signaling(assemb, solver):
    """
    Calculate Signaling Effect \n
    Input:
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
    """

    # Get dimension info. of assemblage and check if it is valid
    try:
        (M, A, N, N1) = np.shape(assemb)

    except ValueError:
        raise ERROR("The dimension of input assemblage is incorrect")

    if N != N1 :
        raise ERROR("The dimension of unnormalized density matrix is incorrect.\nIt should be N x N, with N >= 2")

    # sum all the outcomes of assemb for different measurement
    rho_X = []
    for x in range(M):
        rho = np.zeros([N, N])

        for a in range(A):
            rho = rho + assemb[x][a]

        rho_X.append(rho)

    # start to solve by Semidefinite program
    SP = Problem()
    sigma = HermitianVariable('\sigma', (N, N))
    
    SP.set_objective(
        'min',
        trace(sigma) - 1
    )
    
    # add constraints
    SP.add_constraint( sigma >> 0 )
    SP.add_list_of_constraints(
        [(sigma - rho_X[x]) >> 0 for x in range(M)]
    )

    SP.solve(solver=solver)

    return SP.value