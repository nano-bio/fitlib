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
from math import log

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


def AE_func(alpha = None, offsetfixed = None, linearbackground = False):
    #this function defines fit functions for appearance energies
    #they are of the form b + (x-AE)^a convoluted with a gaussian
    #see file docs/AE_conv.pdf for details

    #0 = offset, 1 is EA, 2 is a factor, 3 is alpha, 4 is the linear background

    base_func = 'lambda p, x: p[2]*sigma**alpha*gamma(alpha+1)*exp(-1.0/(4.0*sigma**2)*(p[1]-x)**2)*fl.pbdv_fa(-(alpha+1), (p[1]-x)/sigma)'

    if alpha is None:
        base_func = base_func.replace('alpha', 'p[3]')
    
    if offsetfixed is None:
        base_func = base_func + ' + p[0]'
    else:
        base_func = base_func + ' + %s' % (offsetfixed)

    if linearbackground is True:
        base_func = base_func + ' + p[4]*x'

    return base_func


# DEPRECATED!
def get_AE_func(sigma, alpha = None, linearbackground = False):
    #this function defines fit functions for appearance energies
    #they are of the form b + (x-AE)^a convoluted with a gaussian
    #see file docs/AE_conv.pdf for details
    
    sigma = sigma / 2*sqrt(2*log(2))
    
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

def get_multiple_AE_func(numonsets, alpha = None):
    #this function defines AE-functions for (numonsets) onsets.
    #see function above
    expr_list = []
    
    '''
    This is slightly different than above. We have to evaluate
    the complete expression outside of this function in the main
    program, because we need p, sigma and alpha.
    Therefore this function just returns the string.
    
    the fittable parameters for the functions are as follows:
    alpha fixed:
    p[0] ... offset or slope (first function) for linear addition (all others)
    p[1] ... AE
    p[2] ... constant (see docs)
    
    alpha fitted:
    as above, but with
    p[3] ... alpha
    '''

    # 10 gaussians ought to be enough for anybody.
    if numonsets < 10:
        for n in range(0, numonsets):
            if alpha is not None:
                if n == 0:
                    expr_list.append('p[%s] + p[%s]*sigma**alpha*gamma(alpha+1)*exp(-1.0/(4.0*sigma**2)*(p[%s]-x)**2)*pbdv_fa(-(alpha+1), (p[%s]-x)/sigma)' % (n*3, n*3 + 2, n*3 + 1, n*3 + 1))
                else:
                    expr_list.append('p[%s]*x + p[%s]*sigma**alpha*gamma(alpha+1)*exp(-1.0/(4.0*sigma**2)*(p[%s]-x)**2)*pbdv_fa(-(alpha+1), (p[%s]-x)/sigma)' % (n*3, n*3 + 2, n*3 + 1, n*3 + 1))
            else:
                if n == 0:
                    expr_list.append('p[%s] + p[%s]*sigma**p[%s]*gamma(p[%s]+1)*exp(-1.0/(4.0*sigma**2)*(p[%s]-x)**2)*pbdv_fa(-(p[%s]+1), (p[%s]-x)/sigma)' % (n*4, n*4 + 2, n*4 + 3, n*4 + 3, n*4 + 1, n*4 + 3, n*4 + 1))
                else:
                    expr_list.append('p[%s]*x + p[%s]*sigma**p[%s]*gamma(p[%s]+1)*exp(-1.0/(4.0*sigma**2)*(p[1%sx)**2)*pbdv_fa(-(p[%s]+1), (p[%s]-x)/sigma)' % (n*4, n*4 + 2, n*4 + 3, n*4 + 3, n*4 + 1, n*4 + 3, n*4 + 1))
    else:
        raise ValueError('Maximum of 10 Onsets are allowed.')

    complete_expr = ' + '.join(expr_list)

    return complete_expr

def gaussfunctions(numpeaks, linear_addition = False):
    #this function defines gauss-shapes for (numpeaks) peaks
    expr_list = []

    # 9 gaussians ought to be enough for anybody.
    if numpeaks < 10:
        for n in range(0, numpeaks):
            expr_list.append('p[%s]*exp(-(p[%s]-x)**2/p[%s])' % (n*3, n*3 + 1, n*3 + 2))
    else:
        raise ValueError('Maximum of 9 Gaussians are allowed.')
        
    # in some cases, one needs a linear addition (ion pair formation)
    if linear_addition is True:
        expr_list.append('p[%s]*x' % (numpeaks*3))

    complete_expr = ' + '.join(expr_list)

    #now define one single lambda
    fitfunc = lambda p, x: eval(complete_expr)

    return fitfunc

def fitES(data, peaks, linear_addition = False):

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

    # in case we want a linear addition (+ p[n]*x)
    # we just guess 1 as the slope
    if linear_addition is True:
        p0.append(1)

    #define fitfunction with n peaks
    fitfunc = gaussfunctions(len(peaks), linear_addition)

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
    #minimal signal to noise ratio for peaks of 2.5 seems to be a good choice
    cwt_peakindices = signal.find_peaks_cwt(data[:, 1], arange(1, datalength / 6), min_snr = 2.5, noise_perc = 25)

    #create array of all maxima found with the CWT method
    cwt_maxima = []
    for peakindex in cwt_peakindices:
        cwt_maxima.append([data[peakindex, 0], data[peakindex, 1]])

    #sort them by their size (signal)
    cwt_maxima.sort(key=itemgetter(1))
    
    peaksfound = len(cwt_maxima)
    
    #have we found the requested number of maxima with the CWT method?
    if peaksfound < numberofpeaks:
        #we haven't found enough. as a backup we use a simple algorithm:
        #search for maxima that dominate the neighbouring 15th of the area
        rel_peakindices = signal.argrelmax(data, order = len(data)/15)
        
        rel_maxima = []
        last = 0
        #now we pick the first of a group that spans a 13th of the whole range
        for peakindex in rel_peakindices[0]:
            if (data[peakindex, 0] - last) > data[:, 0].max()/13:
                last = data[peakindex, 0]
                rel_maxima.append([data[peakindex, 0], data[peakindex, 1]])
        
        #sort them by their size
        rel_maxima.sort(key=itemgetter(1))
        
        #combine the two arrays (rel_maxima, cwt_maxima) to mf giving cwt the preference
        mf = []
        # we go through cwt_maxima and remove all rel_maxima from their list,
        # if they are within one 1 eV
        for cwt_maximum in cwt_maxima:
            for rel_maximum in rel_maxima:
                if abs(rel_maximum[0] - cwt_maximum[0]) < 1:
                    rel_maxima.remove(rel_maximum)
            mf.append(cwt_maximum)
        
        #now we fill up with the largest rel_maxima
        while len(mf) < numberofpeaks:
            #only works if rel_maxima search was successful.
            try:
                mf.append(rel_maxima.pop())
            except:
                break
    else:
        #we didn't use the simple algorithm. all maxima found by CWT
        mf = cwt_maxima
    
    #sort again and define amount of found peaks again
    mf.sort(key=itemgetter(1))
    peaksfound = len(mf)

    #empty array to return
    mf_final = []
    i = 0
    
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
    absolutedatavalues = vecabs(data[:,1])
    weights = vecsqrt(absolutedatavalues)
    # note that we set all zero values to 1, in order to have useful weights
    place(weights, weights==0, 1)

    #fit
    res = optimize.leastsq(errfunc, initial_parameters[:], args = (data[:,0],data[:,1]), diag = weights, full_output = True)
    (p1, pcov, infodict, errmsg, ier) = res

    # bad luck, didn't converge
    if ier not in [1, 2, 3, 4]:
        return None
    else:
        return p1
    
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
