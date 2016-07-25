#!/usr/bin/python

#library for logging stuff

import os
import sys
import datetime

class Log:
    """This class provides logging functions"""
    def __init__(self, filename = None, outputfolder = 'output', tovariable = False, timestamp = False):
    
        #make a timestamp
        self.starttime = datetime.datetime.now()
        
        if timestamp == True:
            self.timestamp = True
        else:
            self.timestamp = False
        
        #shall we log everything to a variable?
        if tovariable == True:
            #we obviously want to write to a variable
            self.logcontent = ''
            self.nologfile = True
            self.tovariable = True

        elif tovariable == False:
            self.outputfolder = outputfolder
            
            #create a default
            if filename is None:
                filename = (self.outputfolder + '/output_%s.log' % self.starttime.strftime('%d_%m_%Y_%Hh%Mm%S'))
        
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
        #we might require a timestamp
        if self.timestamp == True:
            now = str(datetime.datetime.now()) + ': '
        else:
            now = ''
            
        if self.nologfile is True and self.tovariable is False:
            print(now + text + '\r\n')
        elif self.nologfile is  True and self.tovariable is True:
            self.logcontent = self.logcontent + now + text + '\r\n'
        else:
            self.fh.write(now + text + '\r\n')

    def ioerror(self, filename):
        self.write('File %s could not be read.' % filename)

    def stop(self):
        self.emptyline()
        self.stoptime = datetime.datetime.now()
        timedelta = self.stoptime - self.starttime
        self.write('Finished. Processing took %s' % timedelta)
        if self.nologfile is not True:
            self.fh.close()
            
    def setargs(self, cmdargs, printargs = True):
        self.cmdargs = cmdargs
        if printargs is True:
            self.printargs()
            
    def emptyline(self):
        self.write('')
