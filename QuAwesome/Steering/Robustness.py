# imports
import sys
import numpy as np
import qutip as qu 
import cvxopt as cvx
import picos as pic
from picos.expressions.variables import HermitianVariable
from QuAwesome import QuAwesomeError as ERROR
from QuAwesome.Steering.genDeterministic import genDeterministicArray as genD

def steeringRobustness(assemb, solver) :
    """
    Calculate Steering Robustness \n
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
        
        solver - a string of solver (mosek or cvxopt)
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

    # start to solve steerable robustness by Semidefinite program
    SP = pic.Problem()

    # add variable (sigma_lambda)
    s_lamb = []
    for l in range(Nlamb):
        s_lamb.append(
            HermitianVariable('sigma_lambda{0}'.format(l), (N, N))
        )

    # generate Deterministic probability distribution array
    D = genD(M, A)

    # add constraints
    SP.add_list_of_constraints([s_lamb[i] >> 0  for i in range(Nlamb)])
    for x in range(M):
        for a in range(A):
            
            # summation of D(a|X, lambda) * sigma_lambda
            summation = D[0][x][a] * s_lamb[0]
            for l in range(1, Nlamb, 1):
                summation += D[l][x][a] * s_lamb[l]
            
            SP.add_constraint(
                (summation - assemb[x][a]) >> 0
            )

    # summation of sigma_lambda
    summation = s_lamb[0]
    for l in range(1, Nlamb, 1):
        summation += s_lamb[l]

    # find the solution
    SP.set_objective(
        'min',
        (pic.trace(summation) - 1)
    )

    # solve the problem
    SP.solve(solver=solver)
    return SP.value