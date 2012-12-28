#!/usr/bin/python

#library for logging stuff

import os
import sys
import datetime

class Log:
    """This class provides logging functions"""
    def __init__(self, filename = 'output/output.log'):

        #let's assume we can create a logfile
        self.nologfile = False

	    #adjust file path in case we're running on fucking windows
        self.filename = os.path.normcase(filename)
                
        #make an absolute path for fucking windows
        self.filename = os.path.join(os.path.dirname(sys.argv[0]), self.filename)
        
        #see if the given or standard logfile can be opened
        try:
            self.fh = open(self.filename, 'w')
        except IOError:
            #maybe it is supposed to be in a subdirectory and we have to the create said directory
            directory = os.path.dirname(self.filename)
            if not os.path.exists(directory):
                os.makedirs(directory)

            try:
                self.fh = open(self.filename, 'w')
            except:
                self.nologfile = True
                
        #make a timestamp
        self.starttime = datetime.datetime.now()

        if self.nologfile is not True:
            self.write('New logfile created. It is now: %s' % self.starttime)
        else:
            self.write('Warning! Could not create logfile! Giving all output to stdout!')

    def write(self, text):
        if self.nologfile is True:
            print text + '\r\n'
        else:
            self.fh.write(text + '\r\n')

    def ioerror(self, filename):
        self.write('File %s could not be read.' % filename)

    def stop(self):
        self.stoptime = datetime.datetime.now()
        timedelta = self.stoptime - self.starttime
        self.write('Finished. Processing took %s' % timedelta)
        self.fh.close()

    def AE_fit_p(self, params):
        if len(params) == 3:
            self.write('AE: %f (Alpha fixed)' % params[2])
        elif len(params) == 4:
            self.write('Alpha: %s, AE: %s' % (params[2], params[3]))
