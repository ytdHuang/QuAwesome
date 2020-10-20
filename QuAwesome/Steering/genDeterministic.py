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
import numpy as np
from QuAwesome import QuAwesomeError as ERROR

def dec2base(n, b, digit):
    result = []
    for i in range(digit - 1, -1, -1):
        result.append( n // (b ** i) )
        n = n % (b ** i)
    return result

def genDeterministicArray(M, O):
    """
    generates deterministic probability distributions \n
    Inputs :
        M : the number of measurements
        O : the number of outcomes of each measurement
    """

    # number of local hidden variable
    d = O ** M

    # Deterministic Array
    DArray = np.zeros(shape=(d, M, O))

    for l in range(d):
        dec = dec2base(l, O, M)
        for x in range(M):
            for a in range(O):
                DArray[l][x][a] = (dec[x] == a)

    return DArray