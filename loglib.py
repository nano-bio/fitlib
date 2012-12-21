#!/usr/bin/python

#library for logging stuff

class Log:
    """This class provides logging functions"""
    def __init__(self, filename = 'output.log'):
        self.fh = open(filename, 'w')
        self.write('New logfile created')

    def write(self, text):
        self.fh.write(text + '\r\n')
