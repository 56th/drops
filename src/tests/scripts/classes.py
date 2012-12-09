#!/usr/bin/python

"""
This script is part of the testing environment of the DROPS project of the
RWTH-Aachen University Scientific Computing Department

Date: 27.03.11
Version :1.0
"""

__author__ = "Oliver Fortmeier, Alin Bastea and Eva Loch"


#The class TestCL contains all the information needed to run a test
class TestCL:
    #Constructor
    def __init__(self, authorName, testType, execFileName, testName, pathExec, pathParam, numProcs, varList, stringList, status):
        self.authorName     = authorName            #name of the author of the test
        self.testType       = testType              #type of the test
        self.execFileName   = execFileName          #name of the executable
        self.testName       = testName              #name of the test
        self.pathExec       = pathExec              #path to the executable file
        self.pathParam      = pathParam             #path to the parameter file
        self.numProcs       = numProcs              #number of processors (for the serial objects this is set to 1)
        self.varList        = list(varList)         #list of the variable that need to be tested. The list is filled up                                                with objects of the class VariableCL found in classvar.py file
        self.stringList     = list(stringList)      #list of all the strings that need to be tested
        self.status         = status                #status of the run test 0 - OK, 1 - compilation failed, 2 -run failed


#Class that contains the variable information inside the tests
class VariableCL:
    #Constructor
    def __init__(self, varName, varPrefix,  varValue, testValue, varAbsError, varRelError,  varMaxValue,  varMinValue):
        self.varName         = varName               #name of variable
        self.varAppear       = varAppear             #counter for appearance
        self.varPrefix       = varPrefix                # prefix of variable
        self.varValue        = varValue              #value of variable
        self.testValue       = testValue             # value of variable in test run
        self.varAbsError     = varAbsError           #absolute error of the variable
        self.varRelError     = varRelError           #relative error of the variable in %
        self.varMinValue  = varMinValue         # maximal value
        self.varMaxValue   = varMaxValue          # minimal value
    #Constructor
    def __init__(self):
        self.varName         = None                  #name of variable
        self.varAppear       = None                  # counter for appearance
        self.varPrefix       = None                  # prefix of variable
        self.varValue        = None                  # target value of variable
        self.testValue       = None                  # value of variable in test run
        self.varAbsError     = None                  #absolute error of the variable
        self.varRelError     = None                  #relative error of the variable in %
        self.varMinValue  = None                    # maximal value
        self.varMaxValue   = None                   # minimal value
    #Setter method for the name of the variable
    def setName(self, varName):
        self.varName = varName
    #Setter method for the counter of appearance
    def setAppear(self, varAppear):
        self.varAppear = varAppear
    #Setter method for the prefix of the variable
    def setPrefix(self, varPrefix):
        self.varPrefix = varPrefix   
    #Setter method for the value of the variable
    def setValue(self, varValue):
        self.varValue = varValue
    #Setter method for the absolute error of the variable
    def setAbsError(self, varAbsError):
        self.varAbsError = varAbsError
    #Setter method for the relative error of the variable
    def setRelError(self, varRelError):
        self.varRelError(self, varRelError)
        #Setter method for the relative error of the variable
    def setMaxValue(self, varMaxValue):
        self.varMaxValuer(self, varMaxValue)
        #Setter method for the relative error of the variable
    def setMinValue(self, varMinValue):
        self.varMinValue(self, varMinValue)
    #Getter method for the name of the variable
    def getName(self):
        return self.varName
    #Getter method for the counter of appearances
    def getAppear(self):
        return self.varAppear
    #Getter method for the name of the variable
    def getPrefix(self):
        return self.varPrefix
    #Getter method for the value of the variable
    def getValue(self):
        return self.varValue
    #Getter method for the absolute error of the variable
    def getAbsError(self):
        return self.varAbsError
    #Getter method for the relative error of the variable
    def getRelError(self):
        return self.varRelError
    #Getter method for the relative error of the variable
    def getMaxValue(self):
        return self.varMaxValue
    #Getter method for the relative error of the variable
    def getMinValue(self):
        return self.varMinValue


#Class that contains the string information inside the tests
class StringCL:
    #Constructor
    def __init__(self, stringValue, find):
        self.stringValue         = stringValue               #value of string
        self.find                = find
    #Constructor
    def __init__(self):
        self.stringValue         = None                      #value of string
        self.find                = None                      #flag (0 - check if string is not found, 1 - check if string is found)
    #Setter method for the value of the string
    def setValue(self, stringValue):
        self.stringValue = stringValue
    #Getter method for the value of the string
    def getValue(self):
        return self.stringValue
    #Setter method for the find of the string
    def setValue(self, find):
        self.find = find
    #Getter method for the value of the string
    def getValue(self):
        return self.find
