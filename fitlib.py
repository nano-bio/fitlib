#!/usr/bin/python
from numpy import *
from scipy import *
from scipy.special import *

from scipy import optimize
from scipy import signal

from operator import itemgetter

import re
import os

import helplib as hl

from math import sqrt

#now for the interesing part: we only load the matplotlib.pyplot if we are not called externally
#reason for this is the problem with setting the backend after loading the module

if __name__ == '__main__':
    import sys
    sys.exit('This program does not have a main routine is not supposed to be called from command line')
elif __name__ == 'fitlib':
    import matplotlib.pyplot as plt
elif __name__ == 'fitlib.fitlib':
    import matplotlib
    matplotlib.use('PDF')

    #now we can import matplotlib completely, because the backend is now set
    import matplotlib.pyplot as plt


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
    expr_list = []

    # 9 gaussians ought to be enough for anybody.
    if numpeaks < 10:
        for n in range(0, numpeaks):
            expr_list.append('p[%s]*exp(-(p[%s]-x)**2/p[%s])' % (n*3, n*3+1, n*3+2))
    else:
        raise ValueError('Maximum of 9 Gaussians are allowed.')

    complete_expr = ' + '.join(expr_list)

    #now define one single lambda
    fitfunc = lambda p, x: eval(complete_expr)

    return fitfunc

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

def guess_ES_peaks(data, numberofpeaks, offset = None, limit = None):
    #first thing: shall we search from a certain offset?
    if offset is not None:
        data = cutarray(data, lowerlim = offset, upperlim = limit)

    #retrieve amount of entries in array
    datalength = len(data[: ,1])

    #use the CWT method implemented in the signal package of scipy
    #minimal signal to noise ratio for peaks of 2 seems to be a good choice
    peakindices = signal.find_peaks_cwt(data[:, 1], arange(1, datalength / 6), min_snr = 2.5, noise_perc = 25)

    #create array of all maxima found (mf)
    mf = []
    for peakindex in peakindices:
        mf.append([data[peakindex, 0], data[peakindex, 1]])

    #sort them by their size
    mf.sort(key=itemgetter(1))

    #empty array to return
    mf_final = []
    i = 0
    
    peaksfound = len(mf)

    #return the (numberofpeaks) highest peaks
    while i < numberofpeaks:
        #do we even have that many peaks?
        if i < peaksfound:
            mf_final.append(mf.pop())
        #we don't. return zeros.
        else:
            mf_final.append([0.0, 0.0])
        i += 1
    
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
    
def data_from_fit_and_parameters(data, fitfunc, parameters):
    minimumpoints = 20

    #check, whether there is a reasonable amount of points to plot:
    if len(data[:, 0]) > minimumpoints:
        #create empy float array in the size of data
        a = empty_like(data, dtype = float)
    else:
        a = empty((minimumpoints, 2), dtype = float)
    
    #create equidistant x points
    a[:, 0] = linspace(data[:, 0].min(), data[:, 0].max(), len(a[:, 0]))
    #calculate corresponding y values
    a[:, 1] = fitfunc(parameters, a[:, 0])
    
    return a

def plot_fit(data, fitfunc, parameters):
    #use above function create equidistant x points for plotting the fit
    data = data_from_fit_and_parameters(data, fitfunc, parameters)

    #plot
    plt.plot(data[:, 0],data[:, 1], 'r--', linewidth=3)

def fit_function_to_data(data, fitfunc, initial_parameters):
    #data has to be a numpy array
    #fits the function to data[:,0] (as x) and data[:,1] (as y) using the initial_parameters
    #returns an array of parameters of the same size as initial_parameters in case of success
    #returns None if data couldn't be fitted

    # Distance to the target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y
    
    #calculate the sqrt of the data values for the fit weights
    vecsqrt = vectorize(sqrt)
    vecabs = vectorize(abs)
    weights = vecsqrt(vecabs(data[:,1] + 1.0))

    #fit
    p1, success = optimize.leastsq(errfunc, initial_parameters[:], args = (data[:,0],data[:,1]), diag = weights)
    
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
