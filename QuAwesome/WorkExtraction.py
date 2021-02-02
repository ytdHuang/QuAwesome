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
from QuAwesome.Steering.genDeterministic import genDeterministicArray as genD
from qutip import Qobj, sigmaz
from numpy import shape, matmul
from numpy import trace as nptr, real as npreal
from picos import Problem, HermitianVariable, sum, trace, value

class WorkExtraction:
    def __init__(assemb):
        """
        Class of calculating Work Extraction \n
        Inputs:
            assemb - a 4-D array containing the assemblage members 
            Each element is an unnormalized N * N density matrix: sigma_a|X where 'a' denote outcome condition on measurement 'X'
                example : if there are totally M+1 measurements, and each has A+1 outcomes
                    assemblage = [
                        [ sigma_0|0 , sigma_1|0 , ... , sigma_A|0 ],
                        [ sigma_0|1 , sigma_1|1 , ... , sigma_A|1 ],
                            .
                            .
                            .
                        [ sigma_0|M , sigma_1|M , ... , sigma_A|M ]
                    ]

        functions:
            setAssemblage(assemb) - reset assemblage
            showAssemblage()  - print the assemblage
            Quantum()         - calculate quantum work extraction (W)
            Classical(solver) - calculate classical bound of work extraction (W_{classical})
            Witness(solver)   - calculate W - W_{classical}
        """
        self.__M = 0         # number of measurement settings
        self.__A = 0         # number of measurement outcome
        self.__N = 0         # system dimension
        self.__Nlamb = 0     # number of local hidden variable
        self.__assemb = None # assemblage
        self.__F      = None # F_a|x
        
        if assemb != None:
            self.setAssemblage(assemb)
            
    def setAssemblage(self, assemb):
        """
        Inputs:
            assemb - a 4-D array containing the assemblage members 
            Each element is an unnormalized N * N density matrix: sigma_a|X where 'a' denote outcome condition on measurement 'X'
                example : if there are totally M+1 measurements, and each has A+1 outcomes
                    assemblage = [
                        [ sigma_0|0 , sigma_1|0 , ... , sigma_A|0 ],
                        [ sigma_0|1 , sigma_1|1 , ... , sigma_A|1 ],
                            .
                            .
                            .
                        [ sigma_0|M , sigma_1|M , ... , sigma_A|M ]
                    ]
        """
        # Get dimension info. of assemblage and check if it is valid
        try:
            (self.__M, self.__A, self.__N, N1) = np.shape(assemb)

        except ValueError:
            raise ERROR("The dimension of input assemblage is incorrect")

        if (self.__N != N1) or (self.__N != 2):
            raise ERROR("The dimension of unnormalized density matrix is incorrect.\nIt should be 2 x 2")
        else:
            self.__assemb = assemb
            self.__Nlamb  = self.__A ** self.__M
            self.__genF_ax()

    def showAssemblage(self):
        """
        Print the current assemblage
        """
        print('Print Assemblage: sigma_\{a\}|\{x\}')
        for x in range(self.__M):
            for a in range(self.__A):
                print('sigma_{}|{}'.format(a, x))
                print(self.__assemb[x][a], '\n')

    def __genF_ax(self):
        sigma_z = sigmaz()
        self.__F = []
        for x in range(self.__M):
            self.__F.append([])
            for a in range(self.__A):
                E = Qobj(self.__assemb[x][a]).eigenstates(sort='high')
                
                U = (E[1][0].full()).tolist() # first column of unitary
                
                # second column
                U[0].append(E[1][1][0][0][0])
                U[1].append(E[1][1][1][0][0])
                
                U = Qobj(U)
                self.__F[x].append((U * sigma_z * U.dag() - sigma_z).full())

    def Classical(self, solver):
        """
        Calculating classical bound of work extraction
        Inputs:
            solver  - a string of solver (mosek or cvxopt)
        """
        if self.__assemb == None:
            raise ERROR("No assemblage found, please set the assemblage by calling \'setAssemblage(assemb)\' function")

        # reduce state of Bob
        sigma_B = sum([assemb[0][a] for a in range(self.__A)])

        # start to solve steerable robustness by Semidefinite program
        SP = Problem()

        # add variable (sig_lam)
        sig_lam = [HermitianVariable('sig_lam_{}'.format(l), (N, N)) for l in range(self.__Nlamb)]
        
        # generate Deterministic probability distribution array
        D = genD(M, A)
        
        # objective function 
        SP.set_objective(
            'max',
            npreal(sum([D[l][x][a] * trace(self.__F[x][a] * sig_lam[l]) for x in range(self.__M) for a in range(self.__A) for l in range(self.__Nlamb)]))
        )
        
        # add constraint
        ## sum_lam sig_lam == sigma_B
        SP.add_constraint(
            sum([sig_lam[l] for l in range(self.__Nlamb)]) == sigma_B
        )
        
        ## sig_lam >= 0
        SP.add_list_of_constraints(
            [sig_lam[l] >> 0 for l in range(self.__Nlamb)]
        )
        
        # solve the problem
        SP.solve(solver=solver)
        return SP.value / (2 * self.__M)

    def Quantum(self):
        """
        Calculate Quantum work extraction
        """
        w = 0
        for x in range(self.__M):
            for a in range(self.__A):
                w = w + nptr(matmul(self.__F[x][a], self.__assemb[x][a]))
                
        return npreal(w) / (2 * self.__M)

    def Witness(self, solver):
        """
        Calculate W - W_{classical}
        Input:
            solver - solver  - a string of solver (mosek or cvxopt) for calculating classical bound
        """
        return max(self.Quantum() - self.Classical(solver), 0)