#!/usr/bin/python

import report, checktest,  setup,  readtests
import sys, glob, os, string, re

"""
This script is part of the testing environment of the DROPS project of the
RWTH-Aachen University Scientific Computing Department

This script calls all the apropiate scripts for performing all the tests in the specifications files

Valid arguments: parallel serial
If no arguments are given, parallel and serial tests are run.

Date: 23.03.11
Version :1.1
"""

__author__ = "Oliver Fortmeier Alin Bastea and Eva Loch"

def runtest(Test):
    command = "cd ../../"+Test.pathExec +"; "
    if Test.testType == "serial":
        command = command + "./"
    else:
        command = command + "mpirun -np " + str(Test.numProcs) + " "
    command = command + Test.execFileName + " " + Test.pathParam +" &> ../tests/output/" + Test.testName +".out"
    print("Running test " +Test.testName+ "...")
    retCode = os.system(command)
    if (retCode != 0):
        Test.status = 2
    else:
        Test.status = 0
    return 0

def runtests(testlist):
    for test in testlist:
        failedtests= []
        setup.compile(test)
        #Run and gather the output of the tests
        runtest(test)
        #perform the testing logic
        resultOfTest=checktest.main(test)
        report.main(test, resultOfTest)
        if (len(resultOfTest[0]) > 0 or len(resultOfTest[1]) > 0):
            failedtests.append(test.testName)
    return failedtests

def main(argv=sys.argv):
    command = "cd .. ; mkdir -p report; mkdir -p output"
    retcode = os.system(command)
    failedTests = []                    #list for failed tests
    reportfile = '../report/reportfile.txt'    #path to the report file
    reportFile = open(reportfile, "w")  #open the report-file for writing the header of the report
    line = '#=======================TESTING DROPS REPORT=======================\n'
    reportFile.writelines(line)
    reportFile.close()
    #generate the lists of objects containing the tests ----> see Test class in classtest.py
    parallelList = readtests.parallel()
    serialList = readtests.serial()
    if ('parallel' in sys.argv or ('parallel' not in sys.argv and 'serial' not in sys.argv)):
        #Set up DROPS for parallel testing
        setup.parallel()
        #Compile the parallel tests
        failedTests = failedTests + runtests(parallelList)
    if ('serial' in sys.argv or ('parallel' not in sys.argv and 'serial' not in sys.argv)):
        #Set up DROPS for serial testing
        setup.serial()
        #Compile the serial tests
        failedTests = failedTests + runtests(serialList)
    reportFile = open(reportfile, "a")
    line = '#=========================FAILED TESTS==============================\n'
    reportFile.writelines(line)
    for i in range (0, len(failedTests)):
        reportFile.writelines(failedTests[i]+'\n')
    reportFile.close()
    return 0



if __name__ == "__main__":
    main()
