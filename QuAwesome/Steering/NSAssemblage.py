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
from QuAwesome import QuAwesomeError as ERROR

from picos import Problem, HermitianVariable, RealVariable, sum, value
from itertools import  combinations
from numpy import eye

def Map_to_NS_Assemblage(assemb, solver='mosek', **extra_options):
    """
    Map the input assemblage to a new one which satisfies the no-signaling condition \n
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
        extra_options - options for solver
    """

    # convert the type of assemb from Qobj to ndarray
    for x, sigma_x in enumerate(assemb):
        for a in range(len(sigma_x)):
            if isinstance(assemb[x][a], Qobj):
                assemb[x][a] = (assemb[x][a]).full()

    # Get dimension info. of assemblage and check if it is valid
    try:
        (M, A, N, N1) = np.shape(assemb)
    
    except ValueError:
        raise ERROR("The dimension of input assemblage is incorrect")

    if N != N1 :
        raise ERROR("The dimension of unnormalized density matrix is incorrect.\nIt should be N x N, with N >= 2")

    II   = eye(N, dtype=int)
    case = list(combinations(list(range(M)), 2))

    P = Problem()

    # variable and constraints
    mu = []
    sigma = []
    for x in range(M):
        mu.append([])
        sigma.append([])
        for a in range(A):
            # variable mu
            mu[x].append(RealVariable("mu_{0}|{1}".format(a,x)))

            # variable sigma (NS assemblage)
            s = HermitianVariable("sigma_{0}|{1}".format(a,x), (N, N))
            P.add_constraint(s >> 0)
            sigma[x].append(s)

    # constraints for infinity norm
    for x in range(M):
        for a in range(A):
            Op = assemb[x][a] - sigma[x][a]
            P.add_constraint(-mu[x][a] * II << Op)
            P.add_constraint( mu[x][a] * II >> Op)

    # constraints for no-signaling condition
    for c in case:
        P.add_constraint(sum(sigma[c[0]]) == sum(sigma[c[1]]))

    # objective func.
    P.set_objective(
        'min',
        sum([mu[x][a] for x in range(M) for a in range(A)])
    )
    
    # solve
    P.solve(solver=solver, **extra_options)
    
    # get the new assemblage
    for x in range(M):
        for a in range(A):
            sigma[x][a] = value(sigma[x][a], numpy=True)
            
    return sigma