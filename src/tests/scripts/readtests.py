#!/usr/bin/python

import classes
import sys, glob, os, string, re

"""
This script is part of the testing environment of the DROPS project of the
RWTH-Aachen University Scientific Computing Department

This script reads the serial tests from the folder ../specifications/serial and
returns a list containing objects of the class Test defined in classtest.py

Date: 23.03.11
Version :1.0
"""

__author__ = "Oliver Fortmeier, Alin Bastea and Eva Loch"

#path to the specification files
def readspecification( path, pattern='*.ref'):
    DBG= False
    List               = []        
    varList          = []        #intermediate list for the variables in the specification files
    stringList      = []        #intermediate list for the strings in the specification files
    relativeError = 0
    #Read all the specification files in the folder and fill up List with the data
    okCounter = 0               #Counter used for correctlly filling up the variable list
    okStringCounter = 0         #Counter used for correctlly filling up the string list. Counter is increased when reading a string value and a find flag
    lineCounter = 0
    for file in glob.glob( os.path.join(path, pattern) ):
        f = open(file,'r')
        linesOfFile = f.readlines() # lines of the specification file
        f.close()
        if (DBG):
            print("=== Parsing file " + file + " ===")
        for line in linesOfFile:
            lineCounter= lineCounter + 1
            #Eliminating the comments from the specification file
            validLine = re.sub("#.*","",line)
            #Collect test-information and read all the data needed for a test
            lengthvalidLine = len(validLine)
            if (len(validLine) > 1):
                keywordmatch =re.match('^(.*?):', validLine)
                keyword= keywordmatch.group(1)
                value= validLine[keywordmatch.end():].strip()
                if (DBG):
                    print( str(lineCounter) + ": [" + keyword + "] = <" + value +">" )
            else:
                keyword = ""
            # test keyword    
            if     keyword=="Type of test":
                testType = value
                if testType =="serial":
                  numProcs =1 
            elif  keyword=="Name of author":
                authorName = value
            elif  keyword=="Name of executable":
                execFileName = value
            elif   keyword=="Folder of executable":
                pathExec = value
            elif  keyword=="Path of parameter-file":
                pathParam = value
            elif keyword =="Number of processes":
                if (testType == "serial" and value != "1"):
                  print("number of processes for serial test case ignored")
                else:
                  numProcs = value
                # read variables
            if okCounter == 0:
                variable = classes.VariableCL()
            if okStringCounter == 0:
                string = classes.StringCL()
            if keyword== "Name of variable":
                variable.varName = re.match('"(.*)"', value).group(1)
                if (DBG):
                    print("- Reading specifications for variable " + variable.varName)
                okCounter = okCounter + 1
            elif keyword== "Prefix of variable":
                variable.varPrefix = re.match('"(.*)"', value).group(1)
                okCounter = okCounter + 1
            # value may be empty
            elif keyword== "Value of variable":
                if (len(value) > 0):
                    ValueOfVar = float(value)
                    okCounter = okCounter + 1
                else:
                    ValueOfVar = None
                variable.varValue = ValueOfVar
            elif keyword== "Relative error":
                if (len(value) > 0):
                    relativeError = float(value)
                    okCounter = okCounter + 1
                else:
                    relativeError = None
                variable.varRelError = relativeError
            elif keyword== "Absolute error":
                if (len(value) > 0):
                    absoluteError = float(value)
                    okCounter = okCounter + 1
                else:
                    absoluteError = None
                variable.varAbsError = absoluteError
            elif keyword=="Minimal value":
                if (len(value) > 0):
                    minValue = float(value)
                    okCounter = okCounter + 1
                else:
                    minValue = None
                variable.varMinValue = minValue
            elif keyword=="Maximal value":
                if (len(value) > 0):
                    maxValue = float(value)
                    okCounter = okCounter + 1
                else:
                    maxValue = None
                variable.varMaxValue = maxValue
            # strings
            elif keyword== "String value":
                string.stringValue = re.match('"(.*)"', value).group(1)
                okStringCounter = okStringCounter + 1
            elif keyword== "String find":
                if value == "0":
                    string.find = False
                else:
                    string.find = True
                okStringCounter = okStringCounter + 1
            #else:
                # warning
            if okStringCounter == 2:          #if counter is 2 it means that we read a value and a find flag.
                stringList.append(string)#append the new string to the list
                okStringCounter = 0           #reset the counter
            if okCounter == 4:          #if counter is 4 it means that we read a name a value and an error.
                varList.append(variable)#append the new variable to the list
                okCounter = 0           #reset the counter
            if (okStringCounter > 0 and okCounter > 0):
	      print(file + ":" + str(lineCounter) + ": Incomplete variable or string specification above ignored, before reading '" +keyword+ "'.")
        testName = file.split("/")[3].split(".")[0]
        List.append(classes.TestCL(authorName, testType, execFileName, testName, pathExec, pathParam, numProcs, varList, stringList, -1))#because we haven't compiled or run the test yet status is 0
        #Empty the intermediate lists of variables and strings
        varList=[]
        stringList=[]
    return List
    
def serial(pattern='*.ref'):
    return readspecification('../specifications/serial/',pattern)
    
def parallel(pattern='*.ref'):
    return readspecification('../specifications/parallel/',pattern)
    
