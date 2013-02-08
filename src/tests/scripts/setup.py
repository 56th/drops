#!/usr/bin/python

import sys, glob, os, string, re

"""
This script is part of the testing environment of the DROPS project of the
RWTH-Aachen University Scientific Computing Department

This script sets up DROPS for serial computation

Date: 23.03.11
Version :1.0
"""

__author__ = "Oliver Fortmeier, Alin Bastea and Eva Loch"

def replacearch( a,  b):
    confFile = '../../drops.conf'               #configuration file for DROPS
    #read the configuration file
    file = open(confFile)
    confFileLines = file.readlines()
    file.close()
    #modify the architecture to strategy b
    for i,line in enumerate(confFileLines):
        confFileLines[i] = line.replace(a,b)
    #rewrite the configuration file
    file = open(confFile, "w")
    file.writelines(confFileLines)
    file.close()

def serial():
    #clean DROPS and enable serial build
    command = "cd ../../; make clean; cmake -DMPI=0"
    return os.system(command)
    
def parallel():
    #clean DROPS and enable parallel build
    command = "cd ../../; make clean; cmake -DMPI=1"
    return os.system(command)

#Test - parameter, object of class TestCL found in classtest.py file
def compile(Test, compileproc=1):
    #create and then run the command for compiling the executable files
    command = string.join(["cd ../../",Test.pathExec," ; make -j",compileproc," ",Test.execFileName],"")
    retCode = os.system(command)
    if (retCode != 0):#if the return code of the compilation is not 0
        Test.status = 1#set the object's status to 1
    else:
        Test.status = 0#set the object's status to 0
    return 0
