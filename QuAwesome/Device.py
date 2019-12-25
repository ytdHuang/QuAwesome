from QuAwesome.ERROR import ERROR
from qiskit.providers.ibmq.ibmqbackend import IBMQBackend

class Device:
    # constructor
    def __init__(self, qubitNum=1):
        if(not isinstance(qubitNum, int)): ERROR("qubitNum should be an integer at least 1")
        if(qubitNum <= 0): ERROR("qubitNum should be an integer at least 1")
        self.__N           = qubitNum
        self.__QubitConfig = []
        self.__GateTime    = {}
        self.__Gamma1      = []
        self.__Gamma2      = []

    # set IBMQ backend Info.
    def setIBMQBackend(self, backend):
        # check if backend is a real device
        if(not isinstance(backend, IBMQBackend)):
            ERROR("given backend is not a legal IBMQ real device")
        self.__N           = 0
        self.__QubitConfig = []
        self.__GateTime    = {}
        self.__Gamma1      = []
        self.__Gamma2      = []

        # Read Qubit config info. and set Gamma1, Gamma2
        for qubit in backend.properties().qubits:
            # Read Qubit config info.
            config = {}
            for j in qubit:
                # Change time unit into 'ns'
                if(j.unit == 'ms'):
                    config[j.name] = j.value * 1000000
                elif(j.unit == 'µs'):
                    config[j.name] = j.value * 1000
                else:
                    config[j.name] = j.value
            
            self.__N = self.__N + 1
            self.__QubitConfig.append(config)

            # set Gamma1 and Gamma2
            self.__Gamma1.append(1 / config['T1'])
            self.__Gamma2.append(0.5 * (1 / config['T2'] + 0.5 / config['T1']))

        # Read Gate time info.
        for i in backend.properties().gates:
            # Change time unit into 'ns'
            if(i.parameters[1].unit == 'ms'):
                self.__GateTime[i.name] = i.parameters[1].value * 1000000
            elif(i.parameters[1].unit == 'µs'):
                self.__GateTime[i.name] = i.parameters[1].value * 1000
            else:
                self.__GateTime[i.name] = i.parameters[1].value

    def getN(self): return self.__N
    
    def getQubitConfig(self): return self.__QubitConfig

    def getGamma1_list(self): return self.__Gamma1

    def getGamma2_list(self): return self.__Gamma2

    def getGateTime(self): return self.__GateTime

    def setGamma1_list(self, gamma1_list):
        # check if the input of gamma1_list is legal
        if((not isinstance(gamma1_list, list)) or (len(gamma1_list) != self.__N)):
            ERROR("Gamma1_list should be a list-type object with " + str(self.__N) + " qubits")
        else:
            self.__Gamma1 = gamma1_list
        
    def setGamma2_list(self, gamma2_list):
        # check if the input of gamma2_list is legal
        if((not isinstance(gamma2_list, list)) or (len(gamma2_list) != self.__N)):
            ERROR("Gamma2_list should be a list-type object with " + str(self.__N) + " qubits")
        else:
            self.__Gamma2 = gamma2_list

    def setGateTime(self, gateName, time, qubit_list):
        # Check if input is legal
        if((gateName != 'id') and (gateName != 'u1') and (gateName != 'u2') and (gateName != 'u3') and (gateName != 'cx')): 
            ERROR("gateName is incorrect: it can only be 'id', 'u1', 'u2', 'u3', 'cx'")
        if(not isinstance(time, int) and not isinstance(time, float)): ERROR("Value of time is illegal")
        if(not isinstance(qubit_list, list)): ERROR("qubit_list should be a list-type object")

        for q in qubit_list:
            # if gate name is 'cx'
            if(gateName == 'cx'):
                # Check if input is legal
                if(not isinstance(q, tuple) or (len(q) != 2)): 
                    ERROR("When gate name is 'cx', each element of qubit_list should be a tuple with qubit index: (control, target)")
                self.__isLegal(q[0])
                self.__isLegal(q[1])
                
                self.__GateTime[gateName + str(q[0]) + "_" + str(q[1])] = time

            # if gate name isn't 'cx'
            else:
                self.__GateTime[gateName + "_" + str(q)] = time
    
    # Check if the qubit is legal
    def __isLegal(self, qubit):
        if(not isinstance(qubit, int) or (qubit < 0) or (qubit >= self.__N)):
            ERROR("Qubit Index is illegal")