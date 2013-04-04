#!/usr/bin/python
import fitlib as fl
import loglib

import os
import re

import matplotlib.pyplot as plt

def do_SF6_calibration(filelist, showplots = True):

    #instanciate a log object
    log = loglib.Log()

    #define empty list of peaks used for calibration
    calpeaks = []

    #pattern to match something like nES_SF6_bla.ee
    fn_pattern = re.compile('^nES_([A-Z]{1}[a-z]?)+[1-9]*_.*\.ee$')

    #set information for typical calibration peaks
    #dictionary with index fragmentname
    #data-structure: [[listofpeaksvalues], [actualpeakvalues], [fit parameters]]
    #actualpeakvalues are set to 0 in the beginning
    frag_info = {'SF6': [[0], [0], []], 'SF5': [[0], [0], []], 'F': [[5.5,9,11.5], [0,0,0], []], 'F2': [[4,6], [0,0], []]}

    #for plot numbering
    i = 1

    fig1 = plt.figure()

    for file in filelist:
        #is it an expected filename?
        if fn_pattern.search(file):
            #we only want the second part of the filename (e.g. SF6)
            filenameparts = file.split('_')
            fragment = filenameparts[1]

            #read the file
            data = fl.readfile(file)

            #guess peaks
            peaks = fl.guess_ES_peaks(data, len(frag_info[fragment][0]))
            #use guessed values to fit; return peak list to frag_info and add fit function parameters to frag_info
            frag_info[fragment][1], frag_info[fragment][2] = fl.fitES(data, peaks)

            #only plot if we actually had the file (fit parameter is still of type 'list' if this is the case)
            if type(frag_info[fragment][2]) is not list:
                plt.subplot(3,2,i)
                fl.plotES(data, fragment)
                # the /3 stems from the fact that the functions takes the number of peaks, not parameters
                fl.plot_fit(data, fl.gaussfunctions(len(frag_info[fragment][2])/3), frag_info[fragment][2])
                i = i + 1

    plt.tight_layout()
    if showplots is True:
        plt.show()
    #plt.savefig('ES_plots/bla.png')


if __name__ == "__main__":

    #read path and create filelist
    filepath = 'ES_data/cal/'
    filelist = os.listdir(filepath)
    newfilelist = []
    #create a list of complete paths
    for file in filelist:
        newfilelist.append(filepath + file)
    #actually call the main routine
    do_SF6_calibration(filelist, False)
