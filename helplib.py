#!/usr/bin/python

import os
import sys

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
