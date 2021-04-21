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
from QuAwesome import QuAwesomeError as ERROR
from qutip import Qobj, ket2dm
from numpy import shape, ndarray

def MLE(rho, epsilon=0):
    """
    Maximum-Likelihood Estimation
    Rescale a Hermitian matrix to nearest postive semidefinite matrix.
    
    Args:
        rho: a hermitian matrix.
        epsilon: (default: 0) the threshold for setting
            eigenvalues to zero. If epsilon > 0 positive eigenvalues
            below epsilon will also be set to zero.
    
    Returns [ndarray]:
        The input matrix rescaled to have non-negative eigenvalues.
        
    References:
        J Smolin, JM Gambetta, G Smith, Phys. Rev. Lett. 108, 070502(2012). 
        Open access: arXiv:1106.5458 [quant-ph].
    """

    # Get dimension info. of rho and check if it is valid
    if isinstance(rho, ndarray):
        try:
            (N, N1) = shape(rho)

        except ValueError:
            raise ERROR("The dimension of rho is incorrect")
            
    elif isinstance(rho, Qobj):
        (N, N1) = rho.shape

    else:
        raise ERROR("Input rho should be ndarray or Qobj")
        
    if (N != N1) or (N < 2):
        raise ERROR("The dimension of rho is incorrect.\nIt should be N x N, with N >= 2")
    else:
        rho = Qobj(rho)
    
    if epsilon < 0:
        raise ERROR('epsilon must be non-negative.')
        
    # Get the eigenvalues and eigenvectors of rho
    # eigenvalues are sorted in increasing order
    [eigval, eigvec] = rho.eigenstates(sort='low')

    for j in range(N):
        if eigval[j] < epsilon:
            tmp = eigval[j]
            eigval[j] = 0.
            
            # Rescale remaining eigenvalues
            for k in range(j + 1, N):
                eigval[k] += tmp / (N - (j + 1))

    # Build positive matrix from the rescaled eigenvalues
    # and the original eigenvectors

    rho_psd = Qobj([ [0] * N ] * N)
    for j in range(N):
        if eigval[j] > 0:
            rho_psd += eigval[j] * ket2dm(eigvec[j])
        
    # renormalize the matrix
    rho_psd = rho_psd / rho_psd.tr()

    return rho_psd.full()