#!/usr/bin/python

import sys, glob, os, string, re

"""
This script is part of the testing environment of the DROPS project of the
RWTH-Aachen University Scientific Computing Department

Script that wites the report of a test

Date: 23.03.11
Version :1.0
"""

__author__ = "Oliver Fortmeier, Alin Bastea and Eva Loch"

#Test - parameter of type TestCL
#resultOfTest - parameter list of failed variables of Test
def main(Test, resultOfTest):
    reportfile = '../report/reportfile.txt'
    reportFile = open(reportfile, "a")
    line = '#==================================================================\n'
    reportFile.writelines(line)
    line = 'Test: ' + Test.testName + '\n'
    reportFile.writelines(line)
    line = 'Author: ' + Test.authorName + '\n'
    reportFile.writelines(line)
    if (Test.status == 0):
        if len(resultOfTest[0]) > 0:
            for result in resultOfTest[0]:
                line = 'THE VARIABLE : | ' + result.varName + ' | FAILED!\nvalue in test run '
                if result.varAppear != None and result.varAppear != 1:
                    line = line + '[ appearance ' + str(result.varAppear) + ' ] '
                line = line + ' is ' + str(result.testValue) + ', \tshould be ' + str(result.varValue) + '\n'
                reportFile.writelines(line)
        if len(resultOfTest[1]) > 0:
            for result in resultOfTest[1]:
                if result.find == True: # string should have been found, but was not
                    line = 'THE STRING : | ' + result.stringValue + ' | WAS NOT FOUND !\n'
                else: # string should not have been found, but was found
                    line = 'THE STRING : | ' + result.stringValue + ' | WAS FOUND!\n'
                reportFile.writelines(line)
    elif (Test.status == 1):
                line = 'THE TEST : | ' + Test.testName + ' | FAILED TO COMPILE!\n'
                reportFile.writelines(line)
    elif (Test.status == 2):
                line = 'THE TEST : | ' + Test.testName + ' | FAILED TO RUN!\n'
                reportFile.writelines(line)
    line = '#==================================================================\n'
    reportFile.writelines(line)
    reportFile.close()
    return 0

if __name__ == "__main__":
    main()
