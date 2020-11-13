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
#from QuAwesome.exceptions import QuAwesomeError as ERROR

class entry(object):
    # constructor
    def __init__(self, index, A, T='S'):
        if isinstance(index, int):
            self.__I = [[index]]

        elif isinstance(index, list):
            self.__I = index
        
        # Num of outcomes - 1
        self.__A = A - 1

        ### Types of the entry ###
        # Z: zero matrix
        # S: Sigma_ab|xy
        # V: variables
        self.__type = T


    def I(self):
        return self.__I


    def Type(self):
        return self.__type


    def copy(self):
        return entry(self.__I, self.__A, self.__type)


    def tag(self):        
        # zero matrix
        if self.__type == 'Z':
            return 'Z'

        # variable
        elif self.__type == 'V':
            s = ''
            for party in self.__I:
                for i in party:
                    s = s + str(i) + '*'
                s = s[:-1]
                s = s + '&'
            return s[:-1]
            
        # assemblage
        else:
            return self.__I[0][0], self.__I[1][0]


    def __mul__(self, e2):
        # if one of them has index_0
        if self.__I[0][0] == 0:
            return e2.copy()
        if e2.I()[0][0] == 0:
            return entry(self.__I, self.__A, self.__type)
        
        # x == x' and a == a'
        if self.__I[0][0] == e2.I()[0][0]:
            return e2.copy()
        
        # x == x' but a != a'
        elif ((e2.I()[0][0] - 1) // ( self.__A)) == ((self.__I[0][0] - 1) // ( self.__A)):
            return entry(-1, self.__A, 'Z')    # index -1 means orthogonal -> 0

        # x != x' and x != 0, x' != 0
        else: 
            self.__I[0].extend(e2.I()[0])
            return entry(self.__I, self.__A, 'V')

    def __and__(self, e2):
        # zero matrix
        if (self.__type == 'Z') or (e2.Type() == 'Z'):
            return entry(-1, self.__A, 'Z')

        # variable
        elif (self.__type == 'V') or (e2.Type() == 'V'):
            return entry(self.__I + e2.I(), self.__A, 'V')
        
        # assemblage
        else:
            return entry(self.__I + e2.I(), self.__A)