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
from datetime import date
from json import dumps, loads, JSONEncoder
from numpy import integer, floating, ndarray
from qutip import Qobj

class DataProcessor:
    def __init__(self, data=[]):
        """
        ##TODO##
        """
        self.__Data = data

    def __str__(self):
        return str(self.__Data)

    def __repr__(self):
        return self.__str__()

    def add(self, data, key):
        if not isinstance(key, str):
            raise ERROR("The variable \'key\' should be in string-type.")

        if isinstance(self.__Data, dict):
            self.__Data[key] = data
        elif isinstance(self.__Data, list):
            raise ERROR("The stored data is now in list-type, please use \'append(data)\' instead.")
            
    def append(self, data):
        if isinstance(self.__Data, list):
            self.__Data.append(data)
        elif isinstance(self.__Data, dict):
            raise ERROR("The stored data is now in dictionary-type, please use \'add(data, key)\' instead.")
        else:
            self.__Data = [self.__Data, data]

    def getData(self, index=None):
        ##TODO##
        pass
    
    def removeData(self, index=None):
        ##TODO##
        pass

    def resetData(self, data):
        self.__Data = data

    def save(self, filename, dateStamp=False):
        ##TODO##
        pass

    def load(self, filename):
        ##TODO##
        pass

    def show(self):
        # PPrint
        ##TODO##
        pass

    def __Encoder(self, obj):
        """ Special json encoder for Qobj and numpy types """

        if isinstance(obj, integer):
            return int(obj)

        elif isinstance(obj, floating):
            return float(obj)

        elif isinstance(obj, complex):
            return {
                "real": obj.real,
                "imag": obj.imag
            }

        elif isinstance(obj, ndarray):
            return obj.tolist()

        elif isinstance(obj, Qobj):
            return {
                "dims": obj.dims,
                "type": obj.type,
                "data": obj.full()
            }

        return JSONEncoder.default(self, obj)
        
    def __Decoder(self, dct):
        """ Special json decoder for Qobj and numpy types """
        
        # complex number
        if ("real" in dct) or ("imag" in dct): 
            return dct["real"] + 1j*dct["imag"]
        
        # Qobj
        if ("dims" in dct) and ("type" in dct):
            return Qobj(dct["data"], dims=dct["dims"], type=dct["type"])

        return dct 