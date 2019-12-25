import numpy as np
from QuAwesome import ERROR
from QuAwesome.Steering.SDProgram import Robustness_SDP

def TSR(assemblage, solver) :
    """
    Calculate Temporal Steering Robustness \n
    assemblage - a 3x2 matrix with each element is a 2x2 density matrix
        example : [[X-, X+],
                   [Y-, Y+],
                   [Z-, Z+]
                ]
    solver - a string of solver (mosek or cvxopt)
    """
    # check for the dimension of given density matrix
    assem = np.array(assemblage)
    if(assem.shape != (3, 2, 2, 2)):
        ERROR("The assemblage should be a 3x2 matrix with each element is a 2x2 density matrix")

    # run SDP
    return Robustness_SDP(assem, solver)