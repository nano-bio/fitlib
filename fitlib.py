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

# now we can import matplotlib completely, because the backend is now set
import matplotlib.pyplot as plt


def pbdv_fa(x,y):
  # we need this, because b is the derivative of a, which is not needed in fits and annoying when defining fit functions
  a, b = pbdv(x, y)

  return a


def AE_func(alpha = None, offsetfixed = None, linearbackground = False):
   # this function defines fit functions for appearance energies
   # they are of the form b + (x-AE)^a convoluted with a gaussian
   # see file docs/AE_conv.pdf for details

   # 0 = offset, 1 is EA, 2 is a factor, 3 is alpha, 4 is the linear background

   #deprecated, not compliant with cuve_fit method
   #base_func = 'lambda x, p, n: p[2]*sigma**alpha*gamma(alpha+1)*exp(-1.0/(4.0*sigma**2)*(p[1]-x)**2)*fl.pbdv_fa(-(alpha+1), (p[1]-x)/sigma)'


   base_func = 'lambda x, p0, p1, p2, p3, p4: p2*(sigma**alpha)*gamma(alpha+1)*exp(-1.0/(4.0*sigma**2)*(p1-x)**2)*fl.pbdv_fa(-(alpha+1), (p1-x)/sigma)' # with convolution

   #base_func ='lambda x, p0, p1, p2, p3, p4: p2*(x - p1)**p3 * heaviside(x-p3,1)' # no convolution


   if alpha is None:
       base_func = base_func.replace('alpha', 'p3')
    
   if offsetfixed is None:
       base_func = base_func + ' + p0'
   else:
       base_func = base_func + ' + %s' % (offsetfixed)

   if linearbackground is True:
       base_func = base_func + ' + p4*x'

   return base_func

def gaussfunctions(numpeaks, linear_addition = False):
    # this function defines gauss-shapes for (numpeaks) peaks
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
    #a[:, 1] = fitfunc(a[:, 0], parameters) #to be deprecated
    a[:, 1] = fitfunc(a[:, 0], parameters[0], parameters[1], parameters[2], parameters[3],parameters[4])
    
    return a

def plot_fit(data, fitfunc, parameters):
    #use above function create equidistant x points for plotting the fit
    data = data_from_fit_and_parameters(data, fitfunc, parameters)

    #plot
    plt.plot(data[:, 0],data[:, 1], 'r--', linewidth=3)
    #mark AE
    plt.plot(parameters[1],0,'g|', linewidth=10)

def plot_fit_with_testinfo(data, fitfunc, parameters, upperfitlimit, lowerfitlimit):
    #use above function create equidistant x points for plotting the fit
    data = data_from_fit_and_parameters(data, fitfunc, parameters)

    #plot
    plt.plot(data[:, 0],data[:, 1], 'r--', linewidth=3)
    #mark AE
    plt.plot(parameters[1],0,'g|', linewidth=10)

    plt.plot(upperfitlimit, 0, 'y>', linewidth=10)
    plt.plot(lowerfitlimit, 0, 'y<', linewidth=10)

#!! weights has not been implemented
def do_the_fit(fitfunc, data, initial_parameters, weights=[]):
    res = optimize.curve_fit(fitfunc,
                             data[:, 0],
                             data[:, 1],
                             initial_parameters[:],
                             method='dogbox',
                             absolute_sigma=False,
                             bounds=([0, 0, -np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf, np.inf]),
                             #sigma=weights,
                             check_finite=True)

    #return initial_parameters and stddev for the parameters
    return res[0], np.sqrt(np.diag(res[1]))


def fit_function_to_data(data, fitfunc, initial_parameters):
    #data has to be a numpy array
    #fits the function to data[:,0] (as x) and data[:,1] (as y) using the initial_parameters
    #returns an array of parameters of the same size as initial_parameters in case of success
    #returns None if data couldn't be fitted


    # calculate the sqrt of the data values for the fit weights

    vecsqrt = vectorize(lambda x: x**(1/2))
    vecabs = vectorize(abs)
    absolutedatavalues = vecabs((data[:,1]))
    weights = vecsqrt(absolutedatavalues)

    # note that we set all zero values to 1, in order to have useful weights
    place(weights, weights==0, 1)

    #fit
    initial_parameters, stddevs = do_the_fit(fitfunc, data, initial_parameters, weights)

    # copy data for cutting down to significant area.
    cutdata = data
    newae = initial_parameters[1]
    #n = int(sqrt(len(data))/2)
    n = 3

    #fit another n times while closing in on the EA
    for iteration in range(0, n):
        #find the position in data[] where the found ae lies
        ae_pos = 0
        for energy in cutdata[:,0]:
            if energy < newae:
              ae_pos += 1
            else:
              break

        cdl = len(cutdata)

        #We'll cut data to energy entries within some % of the fitted AE in order to get a better fit around that point
        #where are we in the following loop
        e_pos = 0

        #a buffer to save the cut data to
        buff_cutdata = []

         # cut data down to 90% of its data points, where the amount left and right of AE stay proportionate
#        for energy in cutdata[:,0]:
#            if e_pos > 0.1 * ae_pos and e_pos < (0.1 * ae_pos + 0.9*cdl):
#                buff_cutdata.append(cutdata[e_pos,:])
#            e_pos += 1

        cutpercent = 0.90
        cutpercent = cutpercent/2

        # cut data down to cut_percent of its data points, where the amount left and right gets cut down symmectrically around ae_pos
        for energy in cutdata[:,0]:
            #if abs(e_pos - ae_pos)<cdl*cutpercent/2:
            #if e_pos > ae_pos - cdl*cutpercent/2 and e_pos < ae_pos + cdl*cutpercent/2:
            if e_pos >= ae_pos - cdl * cutpercent  and e_pos <= ae_pos + cdl * cutpercent:
                buff_cutdata.append(cutdata[e_pos,:])
            e_pos += 1

        cutdata = np.array(buff_cutdata)

        #fit again to compare initial_parameters[1] before and after
        initial_parameters, stddevs = do_the_fit(fitfunc, cutdata, initial_parameters, weights)

        newae = initial_parameters[1]

    cutdata, initial_parameters, stdevs = cut_to_ev(fitfunc, cutdata, 1, initial_parameters)

    # bad luck, didn't converge
#    if ier not in [1, 2, 3, 4]:
#        return errmsg, None
#    else:
    return None, initial_parameters, stddevs[1], cutdata[0][0], cutdata[len(cutdata)-1][0]
    
def cut_to_ev(fitfunc, cutdata, ev, p):
    # cut data down to +- ev around ae_pos
#!! definitely improve list and numpy.array usage
    e_pos = 0
    ae = p[1]
    ret_cutdata = []

    for energy in cutdata[:, 0]:

        if abs(energy - ae) < ev:
            ret_cutdata.append(np.array(cutdata[e_pos, :]))
        e_pos += 1

    ret_cutdata = np.array(ret_cutdata)

    initial_parameters, stdevs = do_the_fit(fitfunc, cutdata, p)

    return ret_cutdata, initial_parameters, stdevs


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
