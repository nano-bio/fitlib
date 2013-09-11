#!/usr/bin/python
import loglib

import os
import re

from numpy import *
from scipy import *

def do_SF6_calibration(filelist, showplots = True, quadratic = True, outputfile = 'SF6_cal.pdf'):
    #this is tricky. we don't need to set the backed to PDF here before importing, because fitlib does that
    import matplotlib
    import matplotlib.pyplot as plt
    import fitlib as fl

    #instanciate a log object
    log = loglib.Log(tovariable = True)

    #define empty list of peaks used for calibration
    calpeaks = []

    #pattern to match something like nES_SF6_bla.ee
    fn_pattern = re.compile('/nES_([A-Z]{1}[a-z]?)+[1-9]*_.*\.ee$')

    #set information for typical calibration peaks
    #dictionary with index fragmentname
    #data-structure: [[listofpeaksvalues], [actualpeakvalues], [fit parameters], [peak search offset], [peak search end]]
    #actualpeakvalues are set to 0 in the beginning
    frag_info = {'SF6': [[0.0], [0], [], None, [7.0]], 'SF5': [[0.1], [0], [], None, [7.0]], 'F': [[5.5, 9.0, 11.5], [0, 0, 0], [], [5.0], [16.0]], 'F2': [[4.625, 11.49], [0,0], [], None, [15.0]]}

    #for plot numbering
    i = 1

    fitdata = []
    fig1 = plt.figure()

    for file in filelist:
        if fn_pattern.search(file):
            badfile = False
            log.write(file)
            #we only want the second but last part of the filename (e.g. SF6)
            filenameparts = file.split('_')
            fragment = filenameparts[1]

            log.write('Now dealing with fragment: ' + fragment)

            #read the file
            try:
                data = fl.readfile(file)
            except IOError:
                badfile = True
                log.write('Could not read the file: ' + file)

            if badfile is False:

                #guess peaks
                peaks = fl.guess_ES_peaks(data, len(frag_info[fragment][0]), offset = frag_info[fragment][3], limit = frag_info[fragment][4])
                log.write('Result of peak guessing: ' + str(peaks))
                #use guessed values to fit; return peak list to frag_info and add fit function parameters to frag_info
                frag_info[fragment][1], frag_info[fragment][2] = fl.fitES(data, peaks)
                n = 0
                for peakpos in frag_info[fragment][1]:
                    fitdata.append([peakpos, frag_info[fragment][0][n]])
                    n += 1

                #only plot if we actually had the file (fit parameter is still of type 'list' if this is the case)
                if type(frag_info[fragment][2]) is not list:
                    plt.subplot(3,2,i)
                    fl.plotES(data, fragment)
                    # the /3 stems from the fact that the functions takes the number of peaks, not parameters
                    fl.plot_fit(data, fl.gaussfunctions(len(frag_info[fragment][2])/3), frag_info[fragment][2])
                    i = i + 1

    #now we fitted all the calibration scans and have an array with peaks in frag_info

    log.emptyline()
    log.write('Raw data of the fits: ' + str(frag_info))
    log.emptyline()

    #convert to numpy array

    fitdata = array(fitdata, dtype = float)

    #quadratic or linear?

    if quadratic is True:
        fitfunc = lambda p, x: p[0] + p[1]*x + p[2]*x**2
        p = [1]*3
    else:
        fitfunc = lambda p, x: p[0] + p[1]*x
        p = [0]*2
        p[0] = 2
        p[1] = 1

    #we fit the peak positions as they are with where they should be

    print fitdata
    parameters = fl.fit_function_to_data(fitdata, fitfunc, p)

    log.emptyline()
    log.write('Fitted and received the following parameters: ' + str(parameters))
    
    if quadratic is True:
        log.write('Fit function is therefore y = %f + %f*x + %f*x^2' % (parameters[0],  parameters[1], parameters[2]))
    else:
        log.write('Fit function is therefore y = %f + %f*x' % (parameters[0], parameters[1]))

    plt.subplot(3, 2, 5)
    plt.plot(fitdata[:,0], fitdata[:,1], 'bo')
    plt.xlabel('Measured')
    plt.ylabel('Corrected')
    plt.grid(True)
    plt.title('Calibration function')
    fl.plot_fit(fitdata, fitfunc, parameters)

#    plt.tight_layout()
    if showplots is True:
        plt.show()

    plt.savefig(outputfile)

    return parameters, log.logcontent


if __name__ == "__main__":

    #read path and create filelist
    filepath = 'ES_data/cal/'
    filelist = os.listdir(filepath)
    newfilelist = []
    #create a list of complete paths
    for file in filelist:
        newfilelist.append(filepath + file)

    #actually call the main routine
    parameters, logcontent = do_SF6_calibration(newfilelist, showplots = False, quadratic = False)

    print parameters, logcontent
