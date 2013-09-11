#!/usr/bin/python

from numpy import *

import os
import sys

import re

def openfile(filename):
    #adjust file path in case we're running on fucking windows
    filename = os.path.normcase(filename)
	
    #open file
    try:
        #relative path?
        f = open(filename,'r')
        return f
    except:
        #probably not. let us try a full path
        filename = os.path.join(os.path.dirname(sys.argv[0]), filename)
        try:
            f = open(filename,'r')
            return f
        except:
            #ok this thing cannot be read
            raise IOError('Could not read file')

def readfile(filename):

    #create empty list
    a = []
    try:
        f = openfile(filename)
    except IOError:
        raise IOError

    #we need to check if a line is actually useful
    num_tab_num = re.compile('^[0-9]+((\.){1}[0-9]+)?\\t[0-9]+((\.){1}[0-9]+)?.*[\\r]?\\n$')

    #read file, skip comment lines
    for line in f:
        #no comments
        if not line.startswith('#'):
            #only number tabulator number
            if num_tab_num.match(line):
                #strip newline and split by tabulator and append to array
                a.append(line.strip('\r\n').split('\t'))

    #convert list a to float array
    data = array(a,dtype = float)

    if len(data) == 0:
        raise IOError('File did not contain any valid lines')

    #close file
    f.close()

    return data

