#!/usr/bin/python

from numpy import *
from scipy import *

import numpy as np

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
parser.add_argument("--fwhm", "-s", help="Specify the FWHM-resolution in eV. If not specified, 1 eV will be assumed", type = float)
parser.add_argument("--offset", help="Specify the offset. If not specified, an offset of 10 will be assumed", type = float)
parser.add_argument("--energyshift", help="Specifiy a shift for the AE in eV. If not specified, 0 eV will be assumed", type = float)
parser.add_argument("--ea", help="Specify the EA in eV. If not specified, 8 eV will be assumed", type = float)
parser.add_argument("--linearbackground", help="Set this, if you want to fit a linear (non-constant) background.", action = 'store_true')
parser.add_argument("--noshow", help="Do not show the plot windows.", action = 'store_true')
parser.add_argument("--nosave", help="Do not save the plots.", action = 'store_true')
parser.add_argument("--writetoplot", "-w", help="Write AE to plot and mark it.", action = 'store_true')
parser.add_argument("--writefit", help="Write a file with an array that contains the x- and y-values of the fit.", action = 'store_true')
parser.add_argument("--outputfolder", help="This option can be used to output files to a specific directory.")
parser.add_argument('--version', action='version', version='r1')

#parse it
args = parser.parse_args()

#now we can assign variables from the given arguments
fwhm = args.fwhm
alpha = args.alpha
linearbackground = args.linearbackground
energyshift = args.energyshift
offset = args.offset
ea = args.ea

#we need an output folder. default is 'output'
outputfolder = 'output'
if args.outputfolder is not None:
    outputfolder = args.outputfolder.rstrip('/')
    if not os.path.exists(outputfolder):
        os.makedirs(outputfolder)

# derive the log class from the loglib and add a function to write AE parameters
class AELog(loglib.Log):
    def AE_fit_p(self, params, alpha, min, max, linearbackground, energyshift, fwhm, offsetfixed):
        if alpha is not None:
            self.write('AE: %f (Alpha fixed to %f)' % (params[1], alpha))
        else:
            self.write('Alpha: %s, AE: %s' % (params[3], params[1]))

        if offsetfixed is not None:
            self.write('Offset fixed at: %s, Slope: %s' % (offsetfixed, params[4]))
        else:
            self.write('Offset: %s, Slope: %s' % (params[0], params[4]))

        self.write('Energy Resolution was set to %s eV FWHM' % fwhm)

        if energyshift is not None:
            self.write('Energies were shifted by %s eV.' % energyshift)

        if min is not None:
            self.write('Fit was started at %s eV.' % min)

        if max is not None:
            self.write('Fit was ended at %s eV.' % max)

        if linearbackground is True:
            self.write('A linear background (non-constant) was used.')
    def printargs(self):
        if self.cmdargs.filename is not None:
            self.write('AE.py is in filename-mode.')
        elif self.cmdargs.folder is not None:
            self.write('AE.py is in folder-mode.')
        elif self.cmdargs.filelist is not None:
            self.write('AE.py is in filelist-mode.')

        if self.cmdargs.alpha is not None:
            self.write('Alpha was set in the command line to %s.' % self.cmdargs.alpha)

        if self.cmdargs.fwhm is not None:
            self.write('Energy resolution was set in the command line to %s eV FWHM.' % self.cmdargs.fwhm)

        if self.cmdargs.linearbackground is True:
            self.write('AE.py was set to fit a linear background (non-constant) from command line.')

        if self.cmdargs.noshow is True:
            self.write('AE.py was set to not show any plots from command line.')

        if self.cmdargs.nosave is True:
            self.write('AE.py was set to not save any plots from command line.')


log = AELog(outputfolder = outputfolder)

# set arguments in the log object
log.setargs(args)

# this is tricky. if we don't plot, we need to use a different backend for matplotlib, as the normal one crashes if plot is not called
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

# one empty line in the logfile
log.emptyline()

# let's walk our filelist
for file in filelist:
    #this variable is set to false if we encounter a non-readable file
    usefulfile = True

    try:
        data = hl.readfile(file[0])
    except IOError:
        usefulfile = False
        log.ioerror(file[0])

    if usefulfile is True:
        #by default we don't cut away data
        minfit = None
        maxfit = None

        #by default we don't have the offset fixed
        offsetfixed = None

        # we set the default values for initial guesses, if they weren't set in the arguments
        if offset is None:
            offset = float64(10)
        if ea is None:
          ea = float64(8)
        if fwhm is None:
            fwhm = float64(1.0)
        if energyshift is None:
            energyshift = float64(0)

        #if there were arguments specified in the filelist, we set them here
        if len(file) > 2:
            #reset other values that are read from the command line
            alpha = None
            linearbackground = False

            #loop through all given arguments
            for arg in file:
                # remove brackets
                arg = arg.replace('[', '').replace(']', '')
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
                    if arg[0] == 'offsetfixed':
                        offsetfixed = float64(arg[1])
                    if arg[0] == 'ea':
                        ea = float64(arg[1])
                    if arg[0] == 'alpha':
                        alpha = float64(arg[1])
                    if arg[0] == 'fwhm':
                        fwhm = float64(arg[1])
                    if arg[0] == 'linearbackground':
                        if arg[1] == 'True':
                            linearbackground = True
                        elif arg[1] == 'False':
                            linearbackground = False
                    if arg[0] == 'energyshift':
                        energyshift = float64(arg[1])

            # adapt the energyshift to minfit and maxfit
            minfit += energyshift
            maxfit += energyshift

        # if the energies have to be shifted...
        if energyshift != 0:
          # ...we first shift all of the values
          for set in data:
            set[0] += energyshift

        # we want to plot the complete data (and not a subset), so we copy it here
        # in case it gets cut later
        complete_data = data

        #doesn't do anything if minfit and maxfit are None (and they are by default)
        data = fl.cutarray(data, lowerlim = minfit, upperlim = maxfit)

        #if the standard guess of ea is outside the data range, we set it to the middle of the range
        datamin = float64(data[:,0].min())
        datamax = float64(data[:,0].max())
        if (ea < datamin) or (ea > datamax):
            ea = (datamax - datamin) / 2 + datamin

        #if any parameters were set from commandline, we overwrite the settings from the file
        if args.alpha is not None:
            alpha = args.alpha
            log.write('Overwriting alpha from command line!')

        if args.linearbackground is True:
            linearbackground = args.linearbackground
            log.write('Overwriting linear background from command line!')

        if args.fwhm is not None:
            fwhm = args.fwhm
            log.write('Overwriting FWHM from command line!')

        if args.offset is not None:
            offset = args.offset
            log.write('Overwriting offset from command line!')

        if args.ea is not None:
            ea = args.ea
            log.write('Overwriting EA from command line!')

        if args.energyshift is not None:
            energyshift = args.energyshift
            log.write('Overwriting energyshift from command line!')

        #depending on the situation of alpha and the lin background we need different amounts of params
        """
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
        """
        p = [0]*5

        p[0] = offset
        p[1] = ea
        p[2] = 1
        p[3] = 1
        p[4] = 0

        #retrieve function for Appearance Energy - the alpha is None if not specified, hence returning a function with alpha fit-able
        #ae_func = fl.get_AE_func(sigma, alpha, linearbackground)

        sigma = fwhm / (2*sqrt(2*np.log(2)))
        ae_func = fl.AE_func(alpha, offsetfixed, linearbackground)
        ae_func = eval(ae_func)

        #actually fit
        errmsg, p1 = fl.fit_function_to_data(data, ae_func, p)

        #log success
        if p1 is not None:
            log.write('============================')
            log.write('Fitted file %s with success.' % file[0])
            log.AE_fit_p(p1, alpha, minfit, maxfit, linearbackground, energyshift, fwhm, offsetfixed)
        else:
            log.write('Failed with fitting of file %s. Error: %s' % (file[0], errmsg))

        #we need to create a more speaking filename
        additions = ''

        additions += '_fwhm=%s' % str(fwhm)

        if alpha is not None:
            additions += '_alpha=%s' % alpha

        if minfit is not None:
            additions += '_min=%s' % minfit

        if maxfit is not None:
            additions += '_max=%s' % maxfit

        if linearbackground is True:
            additions += '_linearbackground'

        if energyshift is not None:
            additions += '_energyshift=%s' % energyshift

        if offsetfixed is not None:
            additions += '_offsetfixed=%s' % offsetfixed

        # we offer an option to write an array of x- and y-values to a file
        if args.writefit is True:
            fitdata = fl.data_from_fit_and_parameters(data, ae_func, p1)

            arrayfilename = os.path.normcase(os.path.join(os.path.dirname(sys.argv[0]), outputfolder + '/' + file[1] + additions + '_fitarray.txt'))
            hl.writearray(fitdata, arrayfilename)

            log.write('Wrote array to %s' % arrayfilename)

        #we don't even need to plot if we neither save nor show
        if (args.noshow is False) or (args.nosave is False):
            fig1 = plt.figure()
            ae_x = p1[1]
            fl.plotES(complete_data, file[1] + ' / AE = %.2f' % ae_x)
            fl.plot_fit(data, ae_func, p1)

            #we annotate the AE in the plot
            if args.writetoplot is True:
                # which alpha?
                if alpha is not None:
                    alpha_value = alpha
                else:
                    alpha_value = p1[3]

                annotate_string = 'AE = %.2f\n$\\alpha$ = %.3f' %(ae_x, alpha_value)
                fig1.text(0.15, 0.85, annotate_string,
                        verticalalignment='top', horizontalalignment='left',
                        color='green', fontsize=15)

        if args.nosave is False:
            plotfile = os.path.normcase(os.path.join(os.path.dirname(sys.argv[0]), outputfolder + '/' + file[1] + additions + '.pdf'))
            plt.savefig(plotfile, format='pdf')
            log.write('Plotted to file: %s' % plotfile)

#done
log.stop()

#showing the plot happens in the end, so to show all windows at the same time and not block the script execution
if args.noshow is False:
    plt.show()
