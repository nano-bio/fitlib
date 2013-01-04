#!/usr/bin/python

#library for logging stuff

import os
import sys
import datetime

class Log:
    """This class provides logging functions"""
    def __init__(self, filename = None):
    
        #make a timestamp
        self.starttime = datetime.datetime.now()
        
        #create a default
        if filename is None:
            filename = ('output/output_%s.log' % self.starttime.strftime('%d_%m_%Y_%Hh%Mm%S'))
        
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
        if self.nologfile is not True:
            self.fh.close()

    def AE_fit_p(self, params, alpha):
        if alpha is not None:
            self.write('AE: %f (Alpha fixed to %f)' % (params[1], alpha))
        else:
            self.write('Alpha: %s, AE: %s' % (params[3], params[1]))
            
    def printargs(self):
        if self.cmdargs.filename is not None:
            self.write('AE.py is in filename-mode.')
        elif self.cmdargs.folder is not None:
            self.write('AE.py is in folder-mode.')
        elif self.cmdargs.filelist is not None:
            self.write('AE.py is in filelist-mode.')
            
        if self.cmdargs.alpha is not None:
            self.write('Alpha was set in the command line to %s.' % self.cmdargs.alpha)
            
        if self.cmdargs.sigma is not None:
            self.write('Energy resolution was set in the command line to %s eV.' % self.cmdargs.sigma)
            
        if self.cmdargs.linearbackground is True:
            self.write('AE.py was set to fit a linear background (non-constant) from command line.')
              
        if self.cmdargs.noshow is True:
            self.write('AE.py was set to not show any plots from command line.')
            
        if self.cmdargs.nosave is True:
            self.write('AE.py was set to not save any plots from command line.')
            
    def setargs(self, cmdargs, printargs = True):
        self.cmdargs = cmdargs
        if printargs is True:
            self.printargs()
