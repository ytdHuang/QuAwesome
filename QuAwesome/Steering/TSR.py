import numpy as np
from QuAwesome import QuAwesomeError as ERROR
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

    try:
        # run SDP
        return Robustness_SDP(assemblage, solver)
    except IndexError:
        raise ERROR("The assemblage should be a 3x2 matrix with each element is a 2x2 density matrix")