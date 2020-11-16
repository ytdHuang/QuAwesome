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
from QuAwesome.exceptions import QuAwesomeError as ERROR
from QuAwesome.NPASteering.entry import entry
from QuAwesome.NPASteering.readAssemblage import readAssemblage
from numpy import eye, zeros
from picos import Problem, Constant, value
from picos.expressions.variables import ComplexVariable

def isPostQuantum(assemb, solver, returnM=False):
    """
    Check if an assemblage is post quantum.

    Inputs:
        assemb - a dictionary containing the assemblage members 
            Each element is an unnormalized dim * dim density matrix: sigma_ab|xy where 'a(b)' denote outcome condition on measurement 'x(y)'
            example : if there are totally X(Y) measurement settings and has A(B) outcomes:
                assemb = {
                    '0,0|1,1': Sigma_00|11,
                        .
                        .
                        .
                    'a,b|x,y': Sigma_ab|xy,
                        .
                        .
                        .
                    'A,B|X,Y': Sigma_AB|XY
                }

                You can use the following code to construct the assemblage:
                assemb['{},{}|{},{}'.format(a, b, x, y)] = Sigma_ab|xy (unnormallized density matrix)
        
        solver  - a string of solver (mosek or cvxopt)
        returnM - choose to return moment matrix or not [Default as False]
    
    Outputs:
        True/False - Whether the assemblage is Post Quantum or not
        MomentMatrix[Optional] - return moment matrix (type: numpy.array) if 'returnM' is 'True'
    """
    A1, A2, M1, M2, dim, Sigma = readAssemblage(assemb)
    N1 = (A1 - 1) * M1 + 1
    N2 = (A2 - 1) * M2 + 1

    # Build moment matrix
    moment = None
    row    = None
    VCount = 0
    VIndex = []
    VList  = []
    for i in range(N1):
        for k in range(N2):
            for j in range(N1):
                for l in range(N2):
                    # (A_j * A_i) \tensor (B_l * B_k)
                    E = (entry(j, A1) * entry(i, A1)) & (entry(l, A2) * entry(k, A2))
                    T = E.Type()

                    # Zero matrix
                    if T == 'Z':
                        element = Constant('0', zeros(dim, dim))

                    # Variable
                    elif T == 'V':
                        tag   = E.tag()
                        tag_d = reverse(tag)

                        # variable exist
                        if tag in VIndex:
                            element = VList[VIndex.index(tag)]

                        # Hermitian conjugate of variable exist
                        elif tag_d in VIndex:
                            element = VList[VIndex.index(tag_d)].H

                        # New Variable
                        else:
                            element = ComplexVariable('x_{0}'.format(VCount), (dim, dim))
                            VIndex.append(tag)
                            VList.append(element)
                            VCount = VCount + 1

                    # Assemblage
                    elif T == 'S':
                        A_i, B_i = E.tag()

                        # Sigma_r
                        if (A_i == 0) and (B_i == 0):
                            summation = []
                            for a in range(A1):
                                for b in range(A2):
                                    summation.append( Sigma['{0},{1}|{2},{3}'.format(a, b, 1, 1)] )
                            element = sum(summation)

                        # Sigma_b|y (Sum_a)
                        elif A_i == 0:
                            b, y = i2ax(B_i, A2)

                            summation = []
                            for a in range(A1):
                                summation.append( Sigma['{0},{1}|{2},{3}'.format(a, b, 1, y)] )
                            element = sum(summation)

                        # Sigma_a|x (Sum_b)
                        elif B_i == 0:
                            a, x = i2ax(A_i, A1)

                            summation = []
                            for b in range(A2):
                                summation.append( Sigma['{0},{1}|{2},{3}'.format(a, b, x, 1)] )
                            element = sum(summation)
                        
                        # Sigma_ab|xy
                        else:
                            a, x = i2ax(A_i, A1)
                            b, y = i2ax(B_i, A2)
                            element = Sigma['{0},{1}|{2},{3}'.format(a, b, x, y)]

                    # align the elements in each row
                    if (j == 0) and (l == 0): row = element
                    else: row = row & element
            
            # align the rows
            if (i == 0) and (k == 0): moment = row
            else: moment = moment // row
            
    ### Semi-definite program ###
    post = Problem()
    post.add_constraint(moment>>0)
    post.set_objective('find')

    # solve the problem
    solution = post.solve(solver=solver, primals=None)
    if solution.status == 'primal feasible':
        result = False
    elif solution.status == 'primal infeasible':
        result = True
    else:
        raise ERROR("The PICOS problem status:", solution.status)
    
    # return
    if returnM:
        return result, value(moment, numpy=True)
    else:
        return result


# a function to reverse tag string of entry,
# also consider as dagger operation of entry tag
def reverse(S):
    s = ''
    [sa, sb] = S.split('&')
    sa = sa.split('*')[::-1]
    sb = sb.split('*')[::-1]
    for a in sa:
        s = s + '*' + a
    s = s + '&'
    for b in sb:
        s = s + b + '*'
    return s[1:-1]


# transfer index -> a, x
def i2ax(index, A):
    return (index - 1) % (A - 1), (index - 1) // (A - 1) + 1