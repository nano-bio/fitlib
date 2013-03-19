#!/usr/bin/python

from numpy import *
from scipy import *

from scipy.special import *
from scipy import optimize

import loglib
import helplib as hl

import argparse
import os, sys
import re

#to compile on windows
if os.name == 'nt':
    from scipy.sparse.csgraph import _validation

#using the argparse module to make use of command line options
parser = argparse.ArgumentParser(description="Fit Appearance Energy measurements and extract the AE")

#we need at least one thing to fit
filegroup = parser.add_mutually_exclusive_group(required=True)
filegroup.add_argument("--filename", help="Specify a filename to be fitted")
filegroup.add_argument("--folder", help="Specify a folder. All files contained in the folder will be fitted")
filegroup.add_argument("--filelist", help="Specify a file that includes a list of filenames to be fitted")

#those are optional
parser.add_argument("--alpha", "-a", help="Specify the exponent for the fit function. If not specified, it will be fitted as well", type = float)
parser.add_argument("--sigma", "-s", help="Specify the FWHM-resolution in eV. If not specified, 1 eV will be assumed", type = float)
parser.add_argument("--linearbackground", help="Set this, if you want to fit a linear (non-constant) background.", action = 'store_true')
parser.add_argument("--noshow", help="Do not show the plot windows.", action = 'store_true')
parser.add_argument("--nosave", help="Do not save the plots.", action = 'store_true')
parser.add_argument("--outputfolder", help="This option can be used to output files to a specific directory.")
parser.add_argument('--version', action='version', version='r1')

#parse it
args = parser.parse_args()

#now we can assign variables from the given arguments
sigma = args.sigma
alpha = args.alpha
linearbackground = args.linearbackground

#we need an output folder. default is 'output'
outputfolder = 'output'
if args.outputfolder is not None:
    outputfolder = args.outputfolder.rstrip('/')
    if not os.path.exists(outputfolder):
        os.makedirs(outputfolder)
    
#instanciate a log object (now that we know where)
log = loglib.Log(outputfolder = outputfolder)

#set arguments in the log object
log.setargs(args)

#this is tricky. if we don't plot, we need to use a different backend for matplotlib, as the normal one crashes if plot is not called
import matplotlib
if (args.noshow is True) and (args.nosave is False):
    matplotlib.use('PDF')

#now we can import fitlib and matplotlib completely, because the backend is now set
import matplotlib.pyplot as plt
import fitlib as fl

#we need this to make a filelist list, that contains filenames
filelist = []

#here we go through the 3 cases: filename, folder, filelist

#this is the easiest case - add the filename to the filelist
if args.filename is not None:
    filelist.append([args.filename, os.path.basename(args.filename)])

elif args.folder is not None:

    #adjust file path in case we're running on fucking windows
    args.folder = os.path.normcase(args.folder)

    #lets go through that folder and add every filename to the filelist
    caldirlist = os.listdir(args.folder)

    for file in caldirlist:
        filelist.append([os.path.join(args.folder, file), file])

elif args.filelist is not None:

    #in this we have to read the list of filenames in a file
    f = hl.openfile(args.filelist)
    
    #we can also compile a regex for the optional arguments
    argre = re.compile('^[a-z]+=(([0-9]+(\.)?[0-9]?)|(True|False))+$')
    
    for line in f:
        #we split the array by whitespaces or tabs (split tries both)
        line = line.strip('\r\n').split()
        if len(line) > 0:
            #first argument should be the filename + path, therefore appending it to a temp array together with the filename
            linearray = [line[0], os.path.basename(line[0])]
            #first argument of the list out
            del line[0]
            #append rest (could also be zero length, doesn't matter)
            linearray = linearray + line
            #append the whole list to the filelist
            filelist.append(linearray)

#if there are too many plots, we shouldn't show them all
if len(filelist) > 5:
    args.noshow = True
    log.write('Not showing any plots, because there are more than 5 files to deal with.')
    
#one empty line in the logfile
log.emptyline()

#let's walk our filelist
for file in filelist:
    #this variable is set to false if we encounter a non-readable file
    usefulfile = True

    try:
        data = fl.readfile(file[0])
    except IOError:
        usefulfile = False
        log.ioerror(file[0])

    if usefulfile is True:

        #in case we don't want any graphical output whatsoever, we skip this, in order to avoid errors from plt
        if (args.noshow is False) or (args.nosave is False):
            fig1 = plt.figure()
            fl.plotES(data, file[1])
        
        #default values for initial guesses
        offset = 10
        ea = 8
        sigma = float64(1.0)

        #by default we don't cut away data
        minfit = None
        maxfit = None

        
        #if there were arguments specified in the filelist, we set them here
        if len(file) > 2:
            #reset other values that are read from the command line
            alpha = None
            linearbackground = False

            #loop through all given arguments
            for arg in file:
                #are they matching "arg=value" where value is a float or int
                if argre.match(arg):
                    #split them up
                    arg = arg.split('=', 2)
                    if arg[0] == 'min':
                        minfit = arg[1]
                    if arg[0] == 'max':
                        maxfit = arg[1]
                    if arg[0] == 'offset':
                        offset = float64(arg[1])
                    if arg[0] == 'ea':
                        ea = float64(arg[1])
                    if arg[0] == 'alpha':
                        alpha = float64(arg[1])
                    if arg[0] == 'sigma':
                        sigma = float64(arg[1])
                    if arg[0] == 'linearbackground':
                        if arg[1] == 'True':
                            linearbackground = True
                        elif arg[1] == 'False':
                            linearbackground = False
                        
            #doesn't do anything if minfit and maxfit are None (and they are by default)
            data = fl.cutarray(data, lowerlim = minfit, upperlim = maxfit)
            
        #if the standard guess of ea is outside the data range, we set it to the middle of the range
        datamin = float64(data[:,0].min())
        datamax = float64(data[:,0].max())
        if (ea < datamin) or (ea > datamax):
            ea = (datamax - datamin) / 2 + datamin
            
        #if alpha or linearbackground were set from commandline, we overwrite the settings from the file
        if args.alpha is not None:
            alpha = args.alpha
            log.write('Overwriting alpha from command line!')
            
        if args.linearbackground is True:
            linearbackground = args.linearbackground
            log.write('Overwriting linear background from command line!')
            
        if args.sigma is not None:
            sigma = args.sigma
            log.write('Overwriting sigma from command line!')

        #depending on the situation of alpha and the lin background we need different amounts of params
        if (alpha is None) and (linearbackground is False):
            p0 = [0]*4
            p0[0] = offset
            p0[1] = ea
            p0[2] = 1
            p0[3] = 1
        elif (alpha is not None) and (linearbackground is False):
            p0 = [0]*3
            p0[0] = offset
            p0[1] = ea
            p0[2] = 1
        elif (alpha is None) and (linearbackground is True):
            p0 = [0]*5
            p0[0] = offset
            p0[1] = ea
            p0[2] = 1
            p0[3] = 1
            p0[4] = 1
        elif (alpha is not None) and (linearbackground is True):
            p0 = [0]*4
            p0[0] = offset
            p0[1] = ea
            p0[2] = 1
            p0[3] = 1
       
        #retrieve function for Appearance Energy - the alpha is None if not specified, hence returning a function with alpha fit-able
        ae_func = fl.get_AE_func(sigma, alpha, linearbackground)

        #actually fit
        p1 = fl.fit_function_to_data(data, ae_func, p0)

        #log success
        if p1 is not None:
            log.write('============================')
            log.write('Fitted file %s with success.' % file[0])
            log.AE_fit_p(p1, alpha, minfit, maxfit, linearbackground, sigma)
        else:
            log.write('Failed with fitting of file %s.' % file[0])

        #we don't even need to plot if we neither save nor show
        if (args.noshow is False) or (args.nosave is False):
            fl.plot_fit(data, ae_func, p1)
        if args.nosave is False:
            #we need to create a more speaking filename
            additions = ''
            
            additions += '_sigma=%s' % sigma
                        
            if alpha is not None:
                additions += '_alpha=%s' % alpha
            
            if minfit is not None:
                additions += '_min=%s' % minfit
                
            if maxfit is not None:
                additions += '_max=%s' % maxfit
                
            if linearbackground is True:
                additions += '_linearbackground'
                
            plotfile = os.path.normcase(os.path.join(os.path.dirname(sys.argv[0]), outputfolder + '/' + file[1] + additions + '.pdf'))
            plt.savefig(plotfile, format='pdf')
            log.write('Plotted to file: %s' % plotfile)

#done
log.stop()

#showing the plot happens in the end, so to show all windows at the same time and not block the script execution
if args.noshow is False:
    plt.show()
