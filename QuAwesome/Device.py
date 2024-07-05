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
from qiskit_ibm_runtime import IBMBackend

class Device:
    # constructor
    def __init__(self, qubitNum=1):
        """
        A class to store some info. of different quantum devices (computers).

        Functions
            -setIBMQBackend(backend):
                set IBMQ backend as the device

            - getN():
                return qubit number

            - getName():
                return device name

            - getQubitConfig(): 
                return a dictionary of qubits config (info.) 

            - getGamma1_list(): 
                return list of Gamma1 (with unit '1/ns')

            - getGamma2_list():
                return list of Gamma2 (with unit '1/ns')

            - getGateTime():
                return a dictionary of gate times (with unit 'ns')

            - getGateError():
                return a dictionary of gate errors

            - to_dict():
                return a dictionary contains all content

            - show():
                print all of the content
        """

        if(not isinstance(qubitNum, int)): raise ERROR("qubitNum should be an integer at least 1")
        if(qubitNum <= 0): raise ERROR("qubitNum should be an integer at least 1")
        self.__N           = qubitNum
        self.__QubitConfig = []
        self.__GateTime    = {}
        self.__GateError   = {}
        self.__Gamma1      = []
        self.__Gamma2      = []
        self.__Name        = ''

    # set IBMQ backend Info.
    def setIBMQBackend(self, backend):
        # check if backend is a real device
        if(not isinstance(backend, IBMBackend)):
            raise ERROR("given backend is not a legal IBMQ real device")
        
        self.__N           = 0
        self.__QubitConfig = []
        self.__GateTime    = {}
        self.__GateError   = {}
        self.__Gamma1      = []
        self.__Gamma2      = []
        self.__Name        = backend.name
        Property = backend.properties()

        # Read Qubit config info. and set Gamma1, Gamma2
        for qubit in Property.qubits:
            # Read Qubit config info.
            config = {}
            for nduv in qubit:
                # Change time unit into 'ns'
                if(nduv.unit == 'ms'):
                    config[nduv.name] = nduv.value * 1000000
                elif(nduv.unit == 'us'):
                    config[nduv.name] = nduv.value * 1000
                elif(nduv.unit == 'ns') or (nduv.unit == ''):
                    config[nduv.name] = nduv.value
                else:
                    config[nduv.name] = [nduv.value, nduv.unit]
            
            self.__N = self.__N + 1
            self.__QubitConfig.append(config)

            # set Gamma1 and Gamma2
            self.__Gamma1.append(1 / config['T1'])
            self.__Gamma2.append(0.5 * (1 / config['T2'] + 0.5 / config['T1']))

        # Read Gate Error and Gate time info.
        for gate in Property.gates:
            
            # set key value (gate name + qubit index)
            idx = gate.gate
            if idx != 'reset':
                for q in gate.qubits:
                    idx += ('_' + str(q))

                nduv_E, nduv_T = gate.parameters

                # Gate Error
                self.__GateError[idx] = nduv_E.value

                # Gate Time and change time unit into 'ns'
                if(nduv_T.unit == 'ms'):
                    self.__GateTime[idx] = nduv_T.value * 1000000
                elif(nduv_T.unit == 'us'):
                    self.__GateTime[idx] = nduv_T.value * 1000
                elif(nduv_T.unit == 'ns'):
                    self.__GateTime[idx] = nduv_T.value
                else:
                    raise ERROR("unknown unit for gate time " + idx)

    def getN(self): return self.__N

    def getName(self): return self.__Name

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

    def to_dict(self):
        return {
            "Name"       : self.__Name,
            "N"          : self.__N,
            "QubitConfig": self.__QubitConfig,
            "GateTime"   : self.__GateTime,
            "GateError"  : self.__GateError,
            "Gamma1"     : self.__Gamma1,
            "Gamma2"     : self.__Gamma2
        }

    def show(self):
        print(self.to_dict())

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
