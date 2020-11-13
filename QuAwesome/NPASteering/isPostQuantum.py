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
from numpy import eye, zeros, shape
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
    """
    A1, A2, M1, M2, dim, Sigma = checkAssemblage(assemb)
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

    # set default result as False (Quantum realized)
    result = False
    try:
        post.solve(solver=solver)
    except Exception as e:
        if (solver == 'mosek') and (str(e) == 'Code 3: Primal solution state claimed infeasible but optimality is required (primals=True).'):
            result = True
        elif (solver == 'cvxopt') and (str(e) == 'Code 3: Primal solution state claimed empty but optimality is required (primals=True).'):
            result = True
        else:
            raise ERROR(e)
    
    # return
    if returnM:
        return result, value(moment, numpy=True)
    else:
        return result


def checkAssemblage(assemb):
    A1 = []
    A2 = []
    M1 = []
    M2 = []
    keys = assemb.keys()

    # read the index and get A1, A2, M1, M2
    for k in keys:
        try:
            [outcome, measurement] = k.split('|')
            [a, b] = outcome.split(',')
            [x, y] = measurement.split(',')
            
            a, b, x, y = int(a), int(b), int(x), int(y)
            if (a >= 0) and (b >= 0) and (x >= 1) and (y >= 1): 
                if (a not in A1): A1.append(a)
                if (b not in A2): A2.append(b)
                if (x not in M1): M1.append(x)
                if (y not in M2): M2.append(y)
            else: raise ValueError

        except ValueError:
            raise ERROR('Invalid index for assemblage: \'' + k + '\'')

    A1.sort()
    A2.sort()
    M1.sort()
    M2.sort()
    if not (A1 == list(range(len(A1)))):
        raise ERROR('The index \'a\' should start from \'0\' and continuous')
    if not (A2 == list(range(len(A2)))):
        raise ERROR('The index \'b\' should start from \'0\' and continuous')
    if not (M1 == list(range(1, len(M1) + 1, 1))):
        raise ERROR('The index \'x\' should start from \'1\' and continuous')
    if not (M2 == list(range(1, len(M2) + 1, 1))):
        raise ERROR('The index \'y\' should start from \'1\' and continuous')
    
    # check if the index is complete
    for a in A1:
        for b in A2:
            for x in M1:
                for y in M2:
                    if '{0},{1}|{2},{3}'.format(a, b, x, y) not in keys:
                        raise ERROR('Missing assemblage index: \'{0},{1}|{2},{3}\''.format(a, b, x, y))
        
    # check the dimension and create Constant-typye(for PICOS) assemblage
    Sigma = {}
    try:
        for k in keys:
            (dim_Row, dim_Col) = shape(assemb[k])

            if (dim_Row < 2) or (dim_Row != dim_Col): 
                raise ValueError
            else:
                Sigma[k] = Constant('S_' + k, assemb[k])

    except ValueError:
        raise ERROR("The dimension of input assemblage is incorrect")      

    return len(A1), len(A2), len(M1), len(M2), dim_Row, Sigma


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