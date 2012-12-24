#!/usr/bin/python

from numpy import *
from scipy import *

from scipy.special import *
from scipy import optimize

import loglib

import argparse
import os, sys

#instanciate a log object
log = loglib.Log()

#using the argparse module to make use of command line options
parser = argparse.ArgumentParser(description="Fit Appearance Energy measurements and extract the AE")

#we need at least one thing to fit
filegroup = parser.add_mutually_exclusive_group(required=True)
filegroup.add_argument("--filename", help="Specify a filename to be fitted")
filegroup.add_argument("--folder", help="Specify a folder. All files contained in the folder will be fitted")
filegroup.add_argument("--filelist", help="Specify a file that includes a list of filenames to be fitted")

#those are optional
parser.add_argument("--alpha", "-a", help="Specify the exponent for the fit function. If not specified, it will be fitted as well", type = float)
parser.add_argument("--sigma", "-s", help="Specify the FWHM-resolution in eV. If not specified, 1 eV will be assumed", default = 1.0, type = float)
parser.add_argument("--noshow", help="Do not show the plot windows.", action = 'store_true')
parser.add_argument("--nosave", help="Do not save the plots.", action = 'store_true')

#parse it
args = parser.parse_args()

#now we can assign variables from the given arguments
sigma = args.sigma
alpha = args.alpha

#this is tricky. if we don't plot, we need to use a different backend for matplotlib, as the normal one crashes if plot is not called
import matplotlib
if (args.noshow is True) and (args.nosave is False):
    matplotlib.use('PDF')

#now we can import fitlib and matplotlib completely, because the backend is now set
import matplotlib.pyplot as plt
import fitlib as fl

#retrieve function for Appearance Energy - the alpha is None if not specified, hence returning a function with alpha fit-able
ae_func = fl.get_AE_func(sigma, alpha)

#Log if we have a fixed alpha
if args.alpha is not None:
    log.write('Using fixed Alpha: %s' % alpha)

#we need this to make a filelist list, that contains filenames
filelist = []

#this is the easiest case - add the filename to the filelist
if args.filename is not None:
    filelist.append([args.filename, os.path.basename(args.filename)])
elif args.folder is not None:
    #lets go through that folder and add every filename to the filelist
    caldirlist = os.listdir(args.folder)

    for file in caldirlist:
        filelist.append([os.path.join(args.folder,file), file])

#if there are to many plots, we shouldn't show them all
if len(filelist) > 5:
    args.noshow = True

#let's walk our filelist
for file in filelist:
    try:
        data = fl.readfile(file[0])
    except IOError:
        log.ioerror(file[0])

    #in case we don't want any graphical outpot whatsoever, we skip this, in order to avoid errors from plt
    if (args.noshow is False) or (args.nosave is False):
        fig1 = plt.figure()
        fl.plotES(data, file[1])

    #if we fit alpha, we need 4 fit parameters with (very rough) guesses
    if args.alpha is None:
        p0 = [0]*4
        p0[0] = 50
        p0[1] = 1
        p0[2] = 1
        p0[3] = 20
    else:
        p0 = [0]*3
        p0[0] = 50
        p0[1] = 1
        p0[2] = 20

    #actually fit
    p1 = fl.fit_function_to_data(data, ae_func, p0)

    #log success
    if p1 is not None:
        log.write('Fitted file %s with success.' % file[0])
        log.AE_fit_p(p1)
    else:
        log.write('Failed with fitting of file %s.' % file[0])

    #we don't even need to plot if we neither save nor show
    if (args.noshow is False) or (args.nosave is False):
        fl.plot_fit(data, ae_func, p1)
    if args.nosave is False:
        plt.savefig('output/' + file[1] + '.pdf', format='pdf')

#showing the plot happens in the end, so to show all windows at the same time and not block the script execution
if args.noshow is False:
    plt.show()

#done
log.stop()
