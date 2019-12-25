# imports
import sys
import numpy as np
import qutip as qu 
import cvxopt as cvx
import picos as pic
from QuAwesome import QuAwesomeError as ERROR

def __Pauli(x) :
    """
    x = 0, 1, 2 means x, y, z
    """
    if   x == 0 : return qu.sigmax()
    elif x == 1 : return qu.sigmay()
    elif x == 2 : return qu.sigmaz()
    else : raise ERROR("Wrong parameters for Pauli matrix.")

def meas_a_X(a, X) :
    """
    a = 0 or 1 \n
    X = 0, 1, 2 means x, y, z \n
    ******************************\n
    measurement_a_X equals (Alice's  0.5 * (identity +- Pauli matrix) tensor (Bob's identity)
    """
    if (a in range(2)) and (X in range(3)) :
        return qu.tensor((0.5 * (qu.qeye(2) - ((-1) ** a) * __Pauli(x=X))), qu.qeye(2))
    else : raise ERROR("Wrong parameters for Measurement_a_X .")

def __krondel(l1, l2) :
    if l1 == l2 : return 1
    else        : return 0

def __D(lamb, a, x) :
    """
    lamb = 0~7 \n
    a  = 0 or 1 \n
    x  = 0, 1, 2  means x, y, z
    """
    if   lamb == 0 :
        if   x == 0 : return __krondel(a, 0)
        elif x == 1 : return __krondel(a, 0)
        elif x == 2 : return __krondel(a, 0)
    elif lamb == 1 :
        if   x == 0 : return __krondel(a, 0)
        elif x == 1 : return __krondel(a, 0)
        elif x == 2 : return __krondel(a, 1)
    elif lamb == 2 :
        if   x == 0 : return __krondel(a, 0)
        elif x == 1 : return __krondel(a, 1)
        elif x == 2 : return __krondel(a, 0)
    elif lamb == 3 :
        if   x == 0 : return __krondel(a, 0)
        elif x == 1 : return __krondel(a, 1)
        elif x == 2 : return __krondel(a, 1)
    elif lamb == 4 :
        if   x == 0 : return __krondel(a, 1)
        elif x == 1 : return __krondel(a, 0)
        elif x == 2 : return __krondel(a, 0)
    elif lamb == 5 :
        if   x == 0 : return __krondel(a, 1)
        elif x == 1 : return __krondel(a, 0)
        elif x == 2 : return __krondel(a, 1)
    elif lamb == 6 :
        if   x == 0 : return __krondel(a, 1)
        elif x == 1 : return __krondel(a, 1)
        elif x == 2 : return __krondel(a, 0)
    elif lamb == 7 :
        if   x == 0 : return __krondel(a, 1)
        elif x == 1 : return __krondel(a, 1)
        elif x == 2 : return __krondel(a, 1)
    else : raise ERROR("Wrong parameter for D function.")

def Weight_SDP(sigma_a_X, solver) :

    # start to solve steerable weight by Semidefinite program
    SW = pic.Problem()

    # add variable (sigma_lambda)
    s_lamb = []
    for  i in range(8):
        s_lamb.append( SW.add_variable('s_lamb{0}'.format(i), (2, 2), vtype='hermitian') )

    # given sigma_a_X
    S_a_x = pic.new_param('S_a_x', sigma_a_X)

    # add constraints
    SW.add_list_of_constraints([s_lamb[i] >> 0  for i in range(8)])
    for i in range(3):
        for j in range(2):
            SW.add_constraint( (S_a_x[i][j] - ( __D(0, j, i) * s_lamb[0] + __D(1, j, i) * s_lamb[1] + 
                                                __D(2, j, i) * s_lamb[2] + __D(3, j, i) * s_lamb[3] + 
                                                __D(4, j, i) * s_lamb[4] + __D(5, j, i) * s_lamb[5] + 
                                                __D(6, j, i) * s_lamb[6] + __D(7, j, i) * s_lamb[7]
                                              ) 
                               ) >> 0
                            )

    # find the solution
    SW.set_objective('min', (1 - pic.trace(s_lamb[0] + s_lamb[1] + s_lamb[2] + s_lamb[3] + 
                                           s_lamb[4] + s_lamb[5] + s_lamb[6] + s_lamb[7]
                                          )
                            )
                    )

    # solve the problem
    SW.solve(solver=solver)
    return SW.obj_value()

def Robustness_SDP(sigma_a_X, solver) :

    # start to solve steerable weight by Semidefinite program
    SR = pic.Problem()

    # add variable (sigma_lambda)
    s_lamb = []
    for  i in range(8):
        s_lamb.append( SR.add_variable('s_lamb{0}'.format(i), (2, 2), vtype='hermitian') )

    # given sigma_a_X
    S_a_x = pic.new_param('S_a_x', sigma_a_X)

    # add constraints
    SR.add_list_of_constraints([s_lamb[i] >> 0  for i in range(8)])
    for i in range(3):
        for j in range(2):
            SR.add_constraint( (__D(0, j, i) * s_lamb[0] + __D(1, j, i) * s_lamb[1] + 
                                __D(2, j, i) * s_lamb[2] + __D(3, j, i) * s_lamb[3] + 
                                __D(4, j, i) * s_lamb[4] + __D(5, j, i) * s_lamb[5] + 
                                __D(6, j, i) * s_lamb[6] + __D(7, j, i) * s_lamb[7] - S_a_x[i][j]
                               ) >> 0
                            )

    # find the solution
    SR.set_objective('min', (pic.trace(s_lamb[0] + s_lamb[1] + s_lamb[2] + s_lamb[3] + 
                                       s_lamb[4] + s_lamb[5] + s_lamb[6] + s_lamb[7]
                                      ) - 1
                            )
                    )

    # solve the problem
    SR.solve(solver=solver)
    return SR.obj_value()