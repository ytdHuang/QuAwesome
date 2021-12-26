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
from QuAwesome import Device
from datetime import datetime
from json import dump, load
from numpy import integer, floating, ndarray
from qutip import Qobj

class DataManager:
    def __init__(self, data=[]):
        """
        The Data Manager for saving (loading) the objects of numpy-types and Qobj to (from) JSON file
        The input data type can be single-object, list-type, or dictionary-type.
        Note: Please avoid to use the key-value (such as "real", "imag", "dims", "type") for dictionary-type data

        To Create a DataManager object:
            DMObj = DataManager()           # initialize with list-type data
            DMObj = DataManager({})         # initialize with dictionary-type data
            DMObj = DataManager(exist_data) # initialize with an exist data

        Attribute
            - Data (data): return the data
                example: You can get the stored data by either one of the following commands
                    1. d = DMObj.Data
                    2. d = DMObj.data

        Functions
            - save(filename, dateStamp[optional]):
                Save data into the file in current directory
            - load(filename):
                Load data from the file
            - append(obj): 
                Append the obj into stored data (only when the stored data is not in dictionary-type)
            - keys(): 
                Returns a list containing all the keys in the stored data (only when the stored data is dictionary-type)
            - values(): 
                Returns a list of all the values in the stored data (only when the stored data is dictionary-type)
            - show():
                Print the Data

        [Example 1] Using DataManager Objects (DMObj) as List-type
        ----------------------------------------------------------------
        [1] data = ['a', 'b', 'c']
        [2] DMObj = DataManager(data) # create DataManager object
        [3] print(DMObj)              # print ['a', 'b', 'c']
        [4] DMObj.append(6)           # append data             : ['a', 'b', 'c', 6]
        [5] len(DMObj)                # length of data          : 4
        [6] DMObj[1] = 'e'            # set data with index     : ['a', 'e', 'c', 6]
        [7] del DMObj[1:3]            # delete and slicing data : ['a', 6]
        [8] for k in DMObj            # iteration of data

        [Example 2] Using DataManager Objects (DMObj) as Dictionary-type
        ----------------------------------------------------------------
        [1] data = {0: 'x', 'a': 1}
        [2] DMObj = DataManager(data) # create DataManager object
        [3] print(DMObj)              # print {0: 'x', 'a': 1}
        [4] len(DMObj)                # length of data       : 2
        [5] DMObj.keys()              # keys of data         : dict_keys([0, 'a'])
        [6] DMObj.values()            # values of data       : dict_values(['x', 1])
        [7] DMObj["b"] = 'y'          # set data with key    : {0: 'x', 'a': 1, 'b': 'y'}
        [8] del DMObj['a']            # delete data with key : {0: 'x', 'b': 'y'}
        [9] for k in DMObj            # iteration of the data key
        """
        self.__Data = data

    def __str__(self):
        return 'Data Manager Object :\n' + str(self.__Data)

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        return len(self.__Data)

    def __setitem__(self, idx, data):
        try:
            self.__Data[idx] = data

        except IndexError as e:
            raise ERROR(str(e))

        except TypeError as e:
            if isinstance(self.__Data, list):
                raise ERROR("The stored data is now in list-type, indices must be integers or slices; try to use \'append(data)\' or \'DataManagerObj[idx] = value\'.")
            elif isinstance(self.__Data, dict):
                raise ERROR("The stored data is now in dictionary-type, the type of key-value is unavailable")
            else:
                raise ERROR(str(e))
                
        except KeyError as e:
            raise ERROR('Wrong key value:', '\'' +  str(e) + '\'')
            
    def __getitem__(self, idx):
        try:
            return self.__Data[idx]

        except IndexError as e:
            raise ERROR(str(e))

        except TypeError as e:
            if isinstance(self.__Data, list):
                raise ERROR("The stored data is now in list-type, indices must be integers or slices.")
            elif isinstance(self.__Data, dict):
                raise ERROR("The stored data is now in dictionary-type, the type of key-value is unavailable")
            else:
                raise ERROR(str(e))
                
        except KeyError as e:
            raise ERROR('The key value', '\'' +  str(e) + '\'', 'does not exsist')

    def __delitem__(self, idx):
        try:
            del self.__Data[idx]
        except IndexError as e:
            raise ERROR(str(e))
        except KeyError as e:
            raise ERROR('The key value', '\'' +  str(e) + '\'', 'does not exsist')

    def __iter__(self):
        return iter(self.__Data)

    @property
    def Data(self):
        """Returns the Data"""
        return self.__Data
    
    @property
    def data(self):
        """Returns the Data"""
        return self.__Data

    def append(self, obj):
        """
        Append the obj into stored data (only when the stored data is not in dictionary-type)
        """
        if isinstance(self.__Data, list):
            self.__Data.append(obj)
        elif isinstance(self.__Data, dict):
            raise ERROR("The stored data is now in dictionary-type, try to use \'DataManagerObj[key] = value\' instead.")
        else:
            self.__Data = [self.__Data, obj]

    def keys(self):
        """
        Returns a list containing all the keys in the stored data (only when the stored data is dictionary-type)
        """
        if isinstance(self.__Data, dict):
            return self.__Data.keys()
        else:
            raise ERROR("The stored data is not in dictionary-type, no attribute \'keys\'")
    
    def values(self):
        """
        Returns a list of all the values in the stored data (only when the stored data is dictionary-type)
        """
        if isinstance(self.__Data, dict):
            return self.__Data.values()
        else:
            raise ERROR("The stored data is not in dictionary-type, no attribute \'values\'")

    def show(self):
        """
        Print the Data
        """
        print(self.__Data)

    def save(self, filename, dateStamp=False):
        """
        Save data into the file named 'filename.json' in current directory
        
        dataStamp [default as False]:
            True to add an dateStamp at the end of the file -> filename_yyyymmdd.json
        """
        if not isinstance(dateStamp, bool):
            raise ERROR("dateStamp should be True or False")
            
        if dateStamp:
            name = filename + '_' + datetime.today().strftime('%Y%m%d') + '.json'
        else:
            name = filename + '.json'

        # save to file
        with open(name, 'w') as file:
            dump(self.__Data, file, default=self.__Encoder)
        print("Save data to \'" + name + '\' (success)')

    def load(self, filename):
        """
        Load data from the file named 'filename.json' in current directory
        Return: the entire Data
        """

        # load data from file
        with open(filename + '.json', 'r') as file:
            self.__Data = load(file, object_hook=self.__Decoder)
        print("Load data from", '\'' + filename + '.json\'', '(success)', '\n')

        return self.__Data

    def __Encoder(self, obj):
        """ Special json encoder for Qobj, numpy, and Device types """

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

        elif isinstance(obj, Device):
            return obj.to_dict()
            
        else:
            return obj
        
    def __Decoder(self, dct):
        """ Special json decoder for Qobj and numpy types """
        
        # complex number
        if ("real" in dct) or ("imag" in dct): 
            return dct["real"] + 1j*dct["imag"]
        
        # Qobj
        if ("dims" in dct) and ("type" in dct):
            return Qobj(dct["data"], dims=dct["dims"], type=dct["type"])

        return dct 