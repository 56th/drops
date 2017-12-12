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
    retCode = os.system(command)
    if (retCode != 0):
        Test.status = 2
    else:
        Test.status = 0
    return 0

def runtests(testlist, compileproc=1):
    failedtests= []
    i=1
    n_test=len(testlist)
    for test in testlist:
        setup.compile(test, compileproc)
        #Run and gather the output of the tests
        print "Running test " +test.testName+ "...\t(",i,"of",n_test,")"
        runtest(test)
        i=i+1
        #perform the testing logic
        resultOfTest=checktest.main(test)
        report.main(test, resultOfTest)
        if (len(resultOfTest[0]) > 0 or len(resultOfTest[1]) > 0 or test.status > 0):
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


    pattern = '*.ref' # default
    nondefaultpattern = False
    compileproc = "1" # default

    # parse command line arguments
    for arg in sys.argv:
        if re.match('--select=',arg):
            [dummy,pattern] = arg.split('=',1)
            nondefaultpattern = True
        if re.match('--compileproc=',arg):
            [dummy,compileproc] = arg.split('=',1)
        if re.match('--help',arg):
            print "Usage: python testDROPS.py [OPTIONS]"
            print "Run DROPS test suite.\n"
            print "  serial\t\tconsider only serial test cases"
            print "  parallel\t\tconsider only parallel test cases"
            print "  --select=PATTERN\ttest specifications to use, default pattern is *.ref"
            print "  --compileproc=P\tnumber of processes used for compilation"
            print "  --help\t\tprint this message and exit"
            return 0
    print "chosen pattern is ", pattern
    print "number of processes for compilation ", compileproc

    parallelList = readtests.parallel(pattern)
    serialList = readtests.serial(pattern)

    if ('parallel' in sys.argv or ('parallel' not in sys.argv and 'serial' not in sys.argv)):
        if nondefaultpattern:
            print "parallelList =", [t.testName for t in parallelList]
        if (parallelList != []):
            #Set up DROPS for parallel testing
            setup.parallel()
            #Compile the parallel tests
            failedTests = failedTests + runtests(parallelList,compileproc)
    if ('serial' in sys.argv or ('parallel' not in sys.argv and 'serial' not in sys.argv)):
        if nondefaultpattern:
            print "serialList =", [t.testName for t in serialList]
        if (serialList != []):
            #Set up DROPS for serial testing
            setup.serial()
            #Compile the serial tests
            failedTests = failedTests + runtests(serialList,compileproc)
    reportFile = open(reportfile, "a")
    if (len(failedTests) == 0):
        line = '#=======================ALL TESTS PASSED===========================\n'
        reportFile.writelines(line)
    else:
        line = '#=========================FAILED TESTS=============================\n'
        reportFile.writelines(line)
        for i in range (0, len(failedTests)):
            reportFile.writelines(failedTests[i]+'\n')
    reportFile.close()
    return 0



if __name__ == "__main__":
    main()
