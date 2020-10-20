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
from qiskit.providers.ibmq.ibmqbackend import IBMQBackend

class Device:
    # constructor
    def __init__(self, qubitNum=1):
        if(not isinstance(qubitNum, int)): raise ERROR("qubitNum should be an integer at least 1")
        if(qubitNum <= 0): raise ERROR("qubitNum should be an integer at least 1")
        self.__N           = qubitNum
        self.__QubitConfig = []
        self.__GateTime    = {}
        self.__GateError   = {}
        self.__Gamma1      = []
        self.__Gamma2      = []

    # set IBMQ backend Info.
    def setIBMQBackend(self, backend):
        # check if backend is a real device
        if(not isinstance(backend, IBMQBackend)):
            raise ERROR("given backend is not a legal IBMQ real device")
        self.__N           = 0
        self.__QubitConfig = []
        self.__GateTime    = {}
        self.__GateError   = {}
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

        # Read Gate Error and Gate time info.
        for i in backend.properties().gates:
            # Gate Error
            self.__GateError[i.name] = i.parameters[0].value

            # Gate Time and change time unit into 'ns'
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

    def getGateTime(self):  return self.__GateTime

    def getGateError(self): return self.__GateError

    def getGateInfo(self):
        info = {}
        for key in self.__GateTime.keys():
            gate = {
                'time' : self.__GateTime[key],
                'error': self.__GateError[key]
            }
            info[key] = gate
        return info

    def setGamma1_list(self, gamma1_list):
        # check if the input of gamma1_list is legal
        if((not isinstance(gamma1_list, list)) or (len(gamma1_list) != self.__N)):
            raise ERROR("Gamma1_list should be a list-type object with " + str(self.__N) + " qubits")
        else:
            self.__Gamma1 = gamma1_list
        
    def setGamma2_list(self, gamma2_list):
        # check if the input of gamma2_list is legal
        if((not isinstance(gamma2_list, list)) or (len(gamma2_list) != self.__N)):
            raise ERROR("Gamma2_list should be a list-type object with " + str(self.__N) + " qubits")
        else:
            self.__Gamma2 = gamma2_list

    def setGateTime(self, gateName, time, qubit_list):
        # Check if input is legal
        if((gateName != 'id') and (gateName != 'u1') and (gateName != 'u2') and (gateName != 'u3') and (gateName != 'cx')): 
            raise ERROR("gateName is incorrect: it can only be 'id', 'u1', 'u2', 'u3', 'cx'")
        if(not isinstance(time, int) and not isinstance(time, float)): raise ERROR("Value of time is illegal")
        if(not isinstance(qubit_list, list)): raise ERROR("qubit_list should be a list-type object")

        for q in qubit_list:
            # if gate name is 'cx'
            if(gateName == 'cx'):
                # Check if input is legal
                if(not isinstance(q, tuple) or (len(q) != 2)): 
                    raise ERROR("When gate name is 'cx', each element of qubit_list should be a tuple with qubit index: (control, target)")
                self.__isLegal(q[0])
                self.__isLegal(q[1])
                
                self.__GateTime[gateName + str(q[0]) + "_" + str(q[1])] = time

            # if gate name isn't 'cx'
            else:
                self.__GateTime[gateName + "_" + str(q)] = time
    
    # Check if the qubit is legal
    def __isLegal(self, qubit):
        if(not isinstance(qubit, int) or (qubit < 0) or (qubit >= self.__N)):
            raise ERROR("Qubit Index is illegal")