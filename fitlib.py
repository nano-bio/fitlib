#!/usr/bin/python
from numpy import *
from scipy import *
from scipy.special import *

from scipy import optimize
import matplotlib.pyplot as plt

from operator import itemgetter

import re
import os

import helplib as hl

#library for plotting and fitting

def pbdv_fa(x,y):
    # we need this, because b is the derivative of a, which is not needed in fits and annoying when defining fit functions
    a, b = pbdv(x, y)
    return a

def get_AE_func(sigma, alpha = None, linearbackground = False):
    #this function defines fit functions for appearance energies
    #they are of the form b + (x-AE)^a convoluted with a gaussian
    #see file docs/AE_conv.pdf for details
    
    sigma = sigma / 2.35482
    
    #0 = offset, 1 is EA, 2 is a factor, 3 is alpha or the linear bg, 4 is the linear background
    
    if linearbackground is False:
        if alpha is not None:
            fitfunc = lambda p, x: p[0] + p[2]*sigma**alpha*gamma(alpha+1)*exp(-1.0/(4.0*sigma**2)*(p[1]-x)**2)*pbdv_fa(-(alpha+1), (p[1]-x)/sigma)
        else:
            fitfunc = lambda p, x: p[0] + p[2]*sigma**p[3]*gamma(p[3]+1)*exp(-1.0/(4.0*sigma**2)*(p[1]-x)**2)*pbdv_fa(-(p[3]+1), (p[1]-x)/sigma)
    elif linearbackground is True:
        #this one supports a linear background. needs one more parameter
        if alpha is not None:
            fitfunc = lambda p, x: p[0] + p[3]*x + p[2]*sigma**alpha*gamma(alpha+1)*exp(-1.0/(4.0*sigma**2)*(p[1]-x)**2)*pbdv_fa(-(alpha+1), (p[1]-x)/sigma)
        else:
            fitfunc = lambda p, x: p[0] + p[4]*x + p[2]*sigma**p[3]*gamma(p[3]+1)*exp(-1.0/(4.0*sigma**2)*(p[1]-x)**2)*pbdv_fa(-(p[3]+1), (p[1]-x)/sigma)   

    return fitfunc

def gaussfunctions(numpeaks):
    #this function defines gauss-shapes for (numpeaks) peaks
    if numpeaks == 1:
        # fit function for one peak
        fitfunc = lambda p, x: p[0]*exp(-(p[1]-x)**2/p[2])
    elif numpeaks == 2:
        # fit function for two peaks
        fitfunc = lambda p, x: p[0]*exp(-(p[1]-x)**2/p[2]) + p[3]*exp(-(p[4]-x)**2/p[5])
    elif numpeaks == 3:
        # fit function for three peaks
        fitfunc = lambda p, x: p[0]*exp(-(p[1]-x)**2/p[2]) + p[3]*exp(-(p[4]-x)**2/p[5]) + p[6]*exp(-(p[7]-x)**2/p[8])
    elif numpeaks == 4:
        # fit function for four peaks
        fitfunc = lambda p, x: p[0]*exp(-(p[1]-x)**2/p[2]) + p[3]*exp(-(p[4]-x)**2/p[5]) + p[6]*exp(-(p[7]-x)**2/p[8]) + p[9]*exp(-(p[10]-x)**2/p[11])
    elif numpeaks == 5:
        # fit function for five peaks
        fitfunc = lambda p, x: p[0]*exp(-(p[1]-x)**2/p[2]) + p[3]*exp(-(p[4]-x)**2/p[5]) + p[6]*exp(-(p[7]-x)**2/p[8]) + p[9]*exp(-(p[10]-x)**2/p[11]) + p[12]*exp(-(p[13]-x)**2/p[14])

    return fitfunc

def readfile(filename):

    #create empty list
    a = []
    try:
        f = hl.openfile(filename)
    except IOError:
        raise IOError
        
    #we need to check if a line is actually useful
    num_tab_num = re.compile('^[0-9]+((\.){1}[0-9]+)?\\t[0-9]+((\.){1}[0-9]+)?.*[\\r]?\\n$')

    #read file, skip comment lines
    for line in f:
        #no comments
        if not line.startswith('#'):
            #only number tabulator number
            if num_tab_num.match(line):
                #strip newline and split by tabulator and append to array
                a.append(line.strip('\r\n').split('\t'))

    #convert list a to float array
    data = array(a,dtype = float)

    if len(data) == 0:
        raise IOError('File did not contain any valid lines')

    #close file
    f.close()

    return data
    

def fitES(data, peaks):

    #we need a place to put our newly found peaks
    peaksfound = peaks

    i = 0
    p0 = [0]*len(peaks)*3

    for peak in peaks:
        #Initial guess for the parameters

        p0[i*3] = peak[1]
        p0[i*3+1] = peak[0]
        p0[i*3+2] = 0.5

        i += 1

    #define fitfunction with n peaks
    fitfunc = gaussfunctions(len(peaks))

    #fit
    p1 = fit_function_to_data(data, fitfunc, p0)

    #success?
    if p1 is not None:
        i = 0
        for peakfound in peaksfound:
            peaksfound[i] = p1[i*3+1]
            i += 1

    return peaksfound, p1

def guess_ES_peaks(data,numberofpeaks):
    #set some variables
    localmax = 0
    datalength = len(data)
    i = 0
    counterleft = 0
    counterright = 0
    nm = numberofpeaks

    #maxima found = mf
    mf = []

    #now we actually go through the data
    for point in data:
        i += 1

        #we go here, if the new point is higher, than our highest so far
        if point[1] > localmax:
            #this is the new local maximum for the last (counterleft) data points
            counterleft += 1
            if counterleft > (datalength/25):
                #ok, this has been the highest point for the last 5% of the file
                localmax = point[1]
                localmaxx = point[0]
                counterright = 0
        elif point[1] < localmax:
            #we are descending again
            counterright += 1
        if counterright > (datalength/25):
            #ok this was a maximum for 5% of the file in each direction
            mf.append([localmaxx, localmax])
            localmax = 0
            counterleft = 0
            counterright = 0

    #sort them by their size
    mf.sort(key=itemgetter(1))

    #empty array to return
    mf_final = []
    i = 0

    #return the (numberofpeaks) highest peaks
    while i < nm:
        i += 1
        mf_final.append(mf.pop())
    
    #sort them by energy
    mf_final.sort(key=itemgetter(0))

    return mf_final

def plotES(data, title):
    #this function plots an ES
    plt.plot(data[:,0],data[:,1],'b-')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Counts (1/s)')
    plt.grid(True)
    plt.title(title)

def plot_fit(data, fitfunc, parameters):
    #create equidistant x points for plotting the fit
    fitx = linspace(data[:,0].min(), data[:,0].max(), len(data[:,0]))

    #plot
    plt.plot(fitx, fitfunc(parameters, fitx), 'r--', linewidth=3)

def fit_function_to_data(data, fitfunc, initial_parameters):
    #data has to be a numpy array
    #fits the function to data[:,0] (as x) and data[:,1] (as y) using the initial_parameters
    #returns an array of parameters of the same size as initial_parameters in case of success
    #returns None if data couldn't be fitted

    # Distance to the target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y

    #fit
    p1, success = optimize.leastsq(errfunc, initial_parameters[:], args=(data[:,0],data[:,1]))

    if success:
        return p1
    else:
        return None
    
def cutarray(data, lowerlim = None, upperlim = None):
    #this function cuts an array and returns it
    #if lowerlim or upperlim are not defined, the maximum is assumed
        
    if lowerlim is None:
        lowerlim = data[:,0].min()
        
    if upperlim is None:
        upperlim = data[:,0].max()
        
    lowerlim = float64(lowerlim)
    upperlim = float64(upperlim)
    
    newdata = []

    for point in data:
        if (point[0] > lowerlim) and (point[0] < upperlim):
            newdata.append(point)
            
    data = array(newdata, dtype = float)
    
    return data