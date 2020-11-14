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
from numpy import eye, zeros, shape
from picos import Constant

def readAssemblage(assemb):
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