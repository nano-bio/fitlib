#!/usr/bin/python

#library for logging stuff

import os

class Log:
    """This class provides logging functions"""
    def __init__(self, filename = 'output/output.log'):

        #let's assume we can create a logfile
        self.nologfile = False

	    #adjust file path in case we're running on fucking windows
        self.filename = os.path.normcase(filename)

        #see if the given or standard logfile can be opened
        try:
            self.fh = open(filename, 'w')
        except IOError:
            #maybe it is supposed to be in a subdirectory and we have to the create said directory
            directory = os.path.dirname(filename)
            if not os.path.exists(directory):
                os.makedirs(directory)

            try:
                self.fh = open(filename, 'w')
            #still not working
            except IOError:
                self.nologfile = True

        if self.nologfile is not True:
            self.write('New logfile created')
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
        self.write('Finished.')
        self.fh.close()

    def AE_fit_p(self, params):
        if len(params) == 3:
            self.write('AE: %f (Alpha fixed)' % params[2])
        elif len(params) == 4:
            self.write('Alpha: %s, AE: %s' % (params[2], params[3]))
