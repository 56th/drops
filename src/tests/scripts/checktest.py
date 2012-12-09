#!/usr/bin/python

import sys, glob, os, string, re

"""
This script is part of the testing environment of the DROPS project of the
RWTH-Aachen University Scientific Computing Department

Wrapper script for the testing logic.

Date: 23.03.11
Version :1.0
"""

__author__ = "Oliver Fortmeier, Alin Bastea and Eva Loch"

# line is a string in which the number should be found
# see python-docs for module re for this regex
def findNumber(line):
    rem = re.match('[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?', line)
    if rem == None:
        return None
    else:
        number = float(rem.group(0))
    return number

#Test - parameter of type TestCL
def checkNumbers(Test,  linesOutput):
    failedNumbersList = []
    for var in Test.varList:
        var.varAppear = 0
        found = False           #we presume that the variable will not be found and we search for it and its value
        fail  = False           #we presume that the variable will pass the test
        for outputline in linesOutput:
            # compare output to  variable.varPrefix
            if re.search(var.varPrefix,  outputline) != None:
              found = True
              tmpline =  outputline[(re.search(var.varPrefix,  outputline)).end():]
              helpstring = tmpline.strip()
              actualValue = findNumber(helpstring)
              if actualValue == None: #false alarm
                  found = False
                  continue
              
              var.varAppear = var.varAppear + 1
              var.testValue = actualValue
              if (var.varRelError != None):
                    #we deal with a relative error logic
                  Err = abs(actualValue-var.varValue)
                  if Err > var.varRelError*abs(var.varValue):
                      fail =True
              elif (var.varAbsError != None ):
                    #we deal with an absolute error logic
                  absErr =abs(actualValue-var.varValue)
                  if absErr > var.varAbsError:
                      fail = True
              elif (var.varValue == None):
                  if actualValue < var.varMinValue:
                        fail = True
                  if actualValue > var.varMaxValue:
                        fail = True
            if fail:
                break
        if fail:
            failedNumbersList.append(var)
        if found == False:
            failedNumbersList.append(var)
    return failedNumbersList

#The class Test contains all the information needed to run a test
def checkStrings(Test,  linesOutput):
    failedStringsList = []
    for teststring in Test.stringList:
        fail = True
        found = False #suppose we don't find the string
        for outputline in linesOutput:
            if (teststring.stringValue in outputline): #if we find it
                found = True    
                
        if found == True and teststring.find == True:#if string was found and we HAD to find it
            fail = False
        elif found == False and teststring.find == False:#if we don't find it and we were NOT suppose to find it
            fail = False
        # in all other cases the test failed: fail = true
            
        if fail == True:
            failedStringsList.append(teststring)
    return failedStringsList
    
#Test - parameter of type TestCL
def main(Test):
    outputFile = string.join(['../output/',Test.testName,'.out'],'')
    #read the output
    outFile = open(outputFile)
    linesOutput = outFile.readlines()
    outFile.close()
    return checkNumbers(Test, linesOutput),  checkStrings(Test, linesOutput)
if __name__ == "__main__":
    main()
