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