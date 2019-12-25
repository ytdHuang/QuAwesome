import numpy as np
import qutip as qu
from QuAwesome import ERROR
from QuAwesome.Steering.SDProgram import meas_a_X, Weight_SDP

def SSW(state, solver) :
    """
    Calculate Spatial Steering Weight \n
    state - a 4x4 (bipartite state) density matrix
    solver - a string of solver (mosek or cvxopt)
    """
    # check for the dimension of given density matrix
    try :
        state = qu.Qobj(state, dims=[[2, 2], [2, 2]])
    except :
        ERROR("Uncorrect parameter \" state \", it should be a 4x4 density matrix")

    # construct known sigma_a_X
    sigma_a_X = []
    for i in range(3) :
        sigma_a_X.append([])
        for j in range(2) :
            sigma_a_X[i].append( np.matrix( (meas_a_X(j, i) * state).ptrace(1).full() ) )
        
    # run SDP
    return Weight_SDP(sigma_a_X, solver)