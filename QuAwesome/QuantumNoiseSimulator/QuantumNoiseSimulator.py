from QuAwesome import ERROR, Device
import numpy as np
from qutip import Qobj, qeye, basis, sigmax, sigmaz, destroy, tensor, Options, mesolve

class QuantumNoiseSimulator:
##### Constructors #####
    # constructor
    def __init__(self, Q_min, Q_max, device):
        # check if device is an QuAwesome.Device type object
        if(not isinstance(device, Device)):
            ERROR("The device should be an QuAwesome.Device type object")

        # get info. from Device
        self.__N        = device.getN()
        self.__Q_min    = 0
        self.__Q_max    = self.__N - 1
        self.__Gamma1   = device.getGamma1_list()
        self.__Gamma2   = device.getGamma2_list()
        self.__GateTime = device.getGateTime()

        # check if Q_min and Q_max is legal, than reset Q_min, Q_max, and N
        if(Q_min > Q_max): ERROR("Q_min should be smaller than or equal Q_max")
        self.__isLegal(Q_min)
        self.__isLegal(Q_max)
        self.__Q_min = Q_min
        self.__Q_max = Q_max
        self.__N     = self.__Q_max - self.__Q_min + 1

        # set some needed variables
        self.__sx   = sigmax()
        self.__sz   = sigmaz()
        self.__sm   = destroy(2).dag()
        self.__up   = basis(2,0) * basis(2,0).dag()
        self.__down = basis(2,1) * basis(2,1).dag()
        self.__I_list    = [qeye(2) for n in range(self.__N)]
        self.__ODEoption = Options(nsteps=15000, store_states=True,rtol=1e-14,atol=1e-14)

        # set Hamiltonian as 0
        self.__H    = self.__I_list.copy()
        self.__H[0] = sigmax()
        self.__H    = 0 * tensor(self.__H)

        # set initial state (density matrix)
        self.__state = tensor( [basis(2,1) for n in range(self.__N)] )
        self.__state = self.__state * self.__state.dag()

        # set C operator list
        self.__c_op_list = []
        for n in range(self.__N):
            if self.__Gamma1[n] > 0.0:
                sm = self.__I_list.copy()
                sm[n] = self.__sm
                self.__c_op_list.append(np.sqrt(self.__Gamma1[n]) * tensor(sm))                                
            if self.__Gamma2[n] > 0.0:
                sz = self.__I_list.copy()
                sz[n] = self.__sz
                self.__c_op_list.append(np.sqrt(self.__Gamma2[n]) * tensor(sz))

##### Public Functions #####
    def getState(self): return self.__state

    def u1(self, lamb, qubit):
        # check if input is legal
        self.__isLegal(qubit)
        if(not isinstance(lamb, int) and not isinstance(lamb, float)): ERROR("Value of lambda is illegal")
        
        # set operator
        operator = self.__I_list.copy()
        operator[qubit - self.__Q_min] = (basis(2,1) * basis(2,1).dag() 
                                            + np.exp(1.0j * lamb) * basis(2,0) * basis(2,0).dag()
                                        )

        self.__ApplyGate(tensor(operator), self.__Time('u1', qubit))

    def u2(self, phi, lamb, qubit):
        # check if input is legal
        self.__isLegal(qubit)
        if(not isinstance(phi,  int) and not isinstance(phi,  float)): ERROR("Value of phi is illegal")
        if(not isinstance(lamb, int) and not isinstance(lamb, float)): ERROR("Value of lambda is illegal")
        
        # set operator
        operator = self.__I_list.copy()
        operator[qubit - self.__Q_min] = np.sqrt(0.5) * (
                                            basis(2,1) * basis(2,1).dag() 
                                            + np.exp(1.0j * (phi + lamb)) * basis(2,0) * basis(2,0).dag()
                                            - np.exp(1.0j * lamb) * basis(2,1) * basis(2,0).dag()
                                            + np.exp(1.0j * phi)  * basis(2,0) * basis(2,1).dag()
                                        )

        self.__ApplyGate(tensor(operator), self.__Time('u2', qubit))

    def u3(self, theta, phi, lamb, qubit):
        # check if input is legal
        self.__isLegal(qubit)
        if(not isinstance(theta, int) and not isinstance(theta, float)): ERROR("Value of theta is illegal")
        if(not isinstance(phi,   int) and not isinstance(phi,   float)): ERROR("Value of phi is illegal")
        if(not isinstance(lamb,  int) and not isinstance(lamb,  float)): ERROR("Value of lambda is illegal")
        
        # set operator
        operator = self.__I_list.copy()
        operator[qubit - self.__Q_min] = (np.cos(theta / 2.) * basis(2,1) * basis(2,1).dag() 
                                        + np.exp(1.0j * (phi + lamb)) * np.cos(theta / 2.) * basis(2,0) * basis(2,0).dag()
                                        - np.exp(1.0j * lamb) * np.sin(theta / 2.) * basis(2,1) * basis(2,0).dag()
                                        + np.exp(1.0j * phi)  * np.sin(theta / 2.) * basis(2,0) * basis(2,1).dag()
                                    )

        self.__ApplyGate(tensor(operator), self.__Time('u3', qubit))

    def cx(self, control, target):
        # check if input is legal
        self.__isLegal(control)
        self.__isLegal(target)
        if(control == target): ERROR('Control and target qubit should be different')

        # set operator
        down_list = self.__I_list.copy()
        down_list[control - self.__Q_min] = self.__down

        up_list = self.__I_list.copy()
        up_list[control - self.__Q_min] = self.__up

        sx_list = self.__I_list.copy()
        sx_list[target - self.__Q_min]  = self.__sx

        operator = tensor(down_list) + tensor(up_list) * tensor(sx_list)

        self.__ApplyGate(operator, self.__Time('cx', control, target))

    def iden(self, qubit):
        self.__isLegal(qubit)            # check if input is legal
        operator = self.__I_list.copy()  # set operator
        self.__ApplyGate(tensor(operator), self.__Time('id', qubit))

    def x(self, qubit):
        self.u3(np.pi, 0, np.pi, qubit)

    def y(self, qubit):
        self.u3(np.pi, np.pi / 2, np.pi / 2, qubit)

    def z(self, qubit):
        self.u3(0, 0, np.pi, qubit)

    def h(self, qubit):
        self.u3(np.pi / 2, 0, np.pi, qubit)

    def s(self, qubit):
        self.u3(0, 0, np.pi / 2, qubit)

    def sdg(self, qubit):
        self.u3(0, 0, - np.pi / 2, qubit)

    def t(self, qubit):
        self.u3(0, 0, np.pi / 4, qubit)

    def tdg(self, qubit):
        self.u3(0, 0, - np.pi / 4, qubit)

    def cu3(self, theta, phi, lamb, control, target):
        # check if input is legal
        self.__isLegal(control)
        self.__isLegal(target)
        if(control == target): ERROR('Control and target qubit should be different')
        if(not isinstance(theta, int) and not isinstance(theta, float)): ERROR("Value of theta is illegal")
        if(not isinstance(phi,   int) and not isinstance(phi,   float)): ERROR("Value of phi is illegal")
        if(not isinstance(lamb,  int) and not isinstance(lamb,  float)): ERROR("Value of lambda is illegal")

        # set operator
        down_list = self.__I_list.copy()
        down_list[control - self.__Q_min] = self.__down

        up_list = self.__I_list.copy()
        up_list[control - self.__Q_min] = self.__up

        u3_list = self.__I_list.copy()
        u3_list[target - self.__Q_min] = (np.cos(theta / 2.) * basis(2,1) * basis(2,1).dag() 
                                            + np.exp(1.0j * (phi + lamb)) * np.cos(theta / 2.) * basis(2,0) * basis(2,0).dag()
                                            - np.exp(1.0j * lamb) * np.sin(theta / 2.) * basis(2,1) * basis(2,0).dag()
                                            + np.exp(1.0j * phi)  * np.sin(theta / 2.) * basis(2,0) * basis(2,1).dag()
                                        )

        operator = tensor(down_list) + tensor(up_list) * tensor(u3_list)

        # set Gate operation time
        # since cu3 operation can be decompose into three single qubit operation and two CNOT operation
        time = 3 * self.__Time('u3', target) + 2 * self.__Time('cx', control, target)

        self.__ApplyGate(operator, time)

    def measure(self, qubit):
        qubit_list = []
        # check if qubit is legal
        if(isinstance(qubit, int)):
            self.__isLegal(qubit)
            qubit_list.append(qubit - self.__Q_min)
        elif(isinstance(qubit, list)):
            for q in qubit:
                self.__isLegal(q)
                qubit_list.append(q - self.__Q_min)
        else: ERROR('Qubit should be an integer or integer list')
        
        # set probability dict, with binary string keys
        prob = {}
        meas_state = self.__state.ptrace(qubit_list)
        for j in range(2 ** len(qubit_list)):
            prob[('{:0' + str(len(qubit_list)) + 'b}').format(j)] = meas_state[j, j]

        return prob

    #def GeneralGate(self):pass

##### Private Functions #####
    # Check if the qubit is legaltlist = np.linspace(0, self.__GateTime[gateName], self.__GateTime[gateName])
    def __isLegal(self, qubit):
        if(not isinstance(qubit, int) or (qubit < self.__Q_min) or (qubit > self.__Q_max)):
            ERROR("Qubit Index is illegal")

    def __Time(self, gateName, control, target=None):
        # single qubit gate
        if(target == None):
            gateID = gateName + '_' + str(control)
        
        # two qubit gate
        else:
            gateID = gateName + str(control) + "_" + str(target)

        # check if the gate exsist
        if gateID not in self.__GateTime:
            ERROR("The gate '" + gateID + "' doesn't exist.")
        else:
            return self.__GateTime[gateID]

    # Apply gates on states with Lindblad Master Equation
    def __ApplyGate(self, operator, time):
        tlist = np.linspace(0, time, time)
        self.__state = operator * self.__state * operator.dag()

        result = mesolve(self.__H, self.__state, tlist, self.__c_op_list, options = self.__ODEoption)
        self.__state = result.states[-1]