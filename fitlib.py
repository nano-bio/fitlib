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

# now for the interesing part: we only load the matplotlib.pyplot if we are not called externally
# reason for this is the problem with setting the backend after loading the module

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

def pbdv_fa(x, y):
    # we need this, because b is the derivative of a, which is not needed in fits and annoying when defining fit functions
    a, b = pbdv(x, y)

    return a

def AE_func(alpha=None, offsetfixed=None, linearbackground=False):
    # this function defines fit functions for appearance energies
    # they are of the form b + (x-AE)^a convoluted with a gaussian
    # see file docs/AE_conv.pdf for details

    # 0 = offset, 1 is EA, 2 is a factor, 3 is alpha, 4 is the linear background

    # deprecated, not compliant with cuve_fit method
    # base_func = 'lambda x, p, n: p[2]*sigma**alpha*gamma(alpha+1)*exp(-1.0/(4.0*sigma**2)*(p[1]-x)**2)*fl.pbdv_fa(-(alpha+1), (p[1]-x)/sigma)'

    base_func = 'lambda x, p0, p1, p2, p3, p4: p2*(sigma**alpha)*gamma(alpha+1)*exp(-1.0/(4.0*sigma**2)*(p1-x)**2)*fl.pbdv_fa(-(alpha+1), (p1-x)/sigma)'  # with convolution

    # base_func ='lambda x, p0, p1, p2, p3, p4: p2*(x - p1)**p3 * heaviside(x-p3,1)' # no convolution

    if alpha is None:
        base_func = base_func.replace('alpha', 'p3')

    if offsetfixed is None:
        base_func = base_func + ' + p0'
    else:
        base_func = base_func + ' + %s' % (offsetfixed)

    if linearbackground is True:
        base_func = base_func + ' + p4*x'

    return base_func

def gaussfunctions(numpeaks, linear_addition=False):
    # this function defines gauss-shapes for (numpeaks) peaks
    expr_list = []

    # 9 gaussians ought to be enough for anybody.
    if numpeaks < 10:
        for n in range(0, numpeaks):
            expr_list.append('p[%s]*exp(-(p[%s]-x)**2/p[%s])' % (n * 3, n * 3 + 1, n * 3 + 2))
    else:
        raise ValueError('Maximum of 9 Gaussians are allowed.')

    # in some cases, one needs a linear addition (ion pair formation)
    if linear_addition is True:
        expr_list.append('p[%s]*x' % (numpeaks * 3))

    complete_expr = ' + '.join(expr_list)

    # now define one single lambda
    fitfunc = lambda p, x: eval(complete_expr)

    return fitfunc

#not in use
def fitES(data, peaks, linear_addition=False):
    # we need a place to put our newly found peaks
    peaksfound = peaks

    i = 0
    p0 = [0] * len(peaks) * 3

    for peak in peaks:
        # Initial guess for the parameters

        p0[i * 3] = peak[1]
        p0[i * 3 + 1] = peak[0]
        p0[i * 3 + 2] = 0.5

        i += 1

    # in case we want a linear addition (+ p[n]*x)
    # we just guess 1 as the slope
    if linear_addition is True:
        p0.append(1)

    # define fitfunction with n peaks
    fitfunc = gaussfunctions(len(peaks), linear_addition)

    # fit
    p1 = fit_function_to_data(data, fitfunc, p01)

    # success?
    if p1 is not None:
        i = 0
        for peakfound in peaksfound:
            peaksfound[i] = p1[i * 3 + 1]
            i += 1

    return peaksfound, p1

def guess_ES_peaks(data, numberofpeaks, offset=None, limit=None):
    # first thing: shall we search from a certain offset?
    if offset is not None:
        data = cutarray(data, lowerlim=offset, upperlim=limit)

    # retrieve amount of entries in array
    datalength = len(data[:, 1])

    # use the CWT method implemented in the signal package of scipy
    # minimal signal to noise ratio for peaks of 2.5 seems to be a good choice
    cwt_peakindices = signal.find_peaks_cwt(data[:, 1], arange(1, datalength / 6), min_snr=2.5, noise_perc=25)

    # create array of all maxima found with the CWT method
    cwt_maxima = []
    for peakindex in cwt_peakindices:
        cwt_maxima.append([data[peakindex, 0], data[peakindex, 1]])

    # sort them by their size (signal)
    cwt_maxima.sort(key=itemgetter(1))

    peaksfound = len(cwt_maxima)

    # have we found the requested number of maxima with the CWT method?
    if peaksfound < numberofpeaks:
        # we haven't found enough. as a backup we use a simple algorithm:
        # search for maxima that dominate the neighbouring 15th of the area
        rel_peakindices = signal.argrelmax(data, order=len(data) / 15)

        rel_maxima = []
        last = 0
        # now we pick the first of a group that spans a 13th of the whole range
        for peakindex in rel_peakindices[0]:
            if (data[peakindex, 0] - last) > data[:, 0].max() / 13:
                last = data[peakindex, 0]
                rel_maxima.append([data[peakindex, 0], data[peakindex, 1]])

        # sort them by their size
        rel_maxima.sort(key=itemgetter(1))

        # combine the two arrays (rel_maxima, cwt_maxima) to mf giving cwt the preference
        mf = []
        # we go through cwt_maxima and remove all rel_maxima from their list,
        # if they are within one 1 eV
        for cwt_maximum in cwt_maxima:
            for rel_maximum in rel_maxima:
                if abs(rel_maximum[0] - cwt_maximum[0]) < 1:
                    rel_maxima.remove(rel_maximum)
            mf.append(cwt_maximum)

        # now we fill up with the largest rel_maxima
        while len(mf) < numberofpeaks:
            # only works if rel_maxima search was successful.
            try:
                mf.append(rel_maxima.pop())
            except:
                break
    else:
        # we didn't use the simple algorithm. all maxima found by CWT
        mf = cwt_maxima

    # sort again and define amount of found peaks again
    mf.sort(key=itemgetter(1))
    peaksfound = len(mf)

    # empty array to return
    mf_final = []
    i = 0

    # return the (numberofpeaks) highest peaks
    while i < numberofpeaks:
        # do we even have that many peaks?
        if i < peaksfound:
            mf_final.append(mf.pop())
        # we don't. return zeros.
        else:
            mf_final.append([0.0, 0.0])
        i += 1

    # sort them by energy
    mf_final.sort(key=itemgetter(0))

    return mf_final

def write_fits(data, p, filename, err, lowerfitbound=-1, upperfitbound=-1, fwhm=-1):
    # use above function create equidistant x points for plotting the fit

    fitfunc = eval_fit_function(fwhm)

    fitdata = data_from_fit_and_parameters(data, fitfunc, p)

    f1 = open(filename.replace('.txt','') + '_AE_%.2f_' % (p[1]) + 'fit_1.txt', "w+")

    for i in range(0, len(data)):
        f1.write("%f\t%f\t%f\r\n" % (data[i][0], data[i][1], fitdata[i][1]))


    # show another zoomed in fit, option to hide/show it might be interesting
    if lowerfitbound != -1 and upperfitbound != -1:

        if lowerfitbound == upperfitbound:
            lowerfitbound = fitdata[0, 0]
            upperfitbound = fitdata[len(fitdata)-1, 0]

            lb = p[1]-1.5
            ub = p[1]+1.5
        else:
            lb = lowerfitbound
            ub = upperfitbound

        # second subplot at the fit relevant area

        f2 = open(filename.replace('.txt', '') + '_AE_%.2f_' % (p[1]) + 'fit_2.txt', "w+")

        cutdata = cutarray(data, lb, ub)
        cut_fit_data = cutarray(fitdata, lb, ub)

        for i in range(0, len(cutdata)):
            f2.write("%f\t%f\t%f\r\n" % (cutdata[i][0], cutdata[i][1], cut_fit_data[i][1]))

        #third subplot around +- 1 eV around AE

        #show 1 eV around AE
        lb = p[1]-0.5
        ub = p[1]+0.5

        if lb > 0 and ub > lb:
          cutdata = cutarray(data, lb, ub)
          cut_fit_data = cutarray(fitdata, lb, ub)

          f3 = open(filename.replace('.txt', '') + '_AE_%.2f_' % (p[1]) + 'fit_3.txt', "w+")

          for i in range(0, len(cutdata)):
              f3.write("%f\t%f\t%f\r\n" % (cutdata[i][0], cutdata[i][1], cut_fit_data[i][1]))

    return f1, f2, f3

#obsolete, plot_fit used
def plotES(data, title):
    # this function plots an ES
    plt.plot(data[:, 0], data[:, 1], 'b.')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Counts (1/s)')
    plt.grid(True)
    plt.title(title)

def plot_fit(data, fitfunc, p, filename, err, lowerfitbound=-1, upperfitbound=-1, fwhm=-1):
    # use above function create equidistant x points for plotting the fit
    fitdata = data_from_fit_and_parameters(data, fitfunc, p)

    #testing
    fitdata = difference_data_fit(data, fitfunc, p)

    fig, axs = plt.subplots(3, 1, constrained_layout=True, figsize=(6, 8))

    c_ax = axs.flatten()[0]

    c_ax.plot(data[:, 0], data[:, 1], 'o', markersize=3, markeredgecolor='b', markerfacecolor='none')

    # plot
    c_ax.plot(fitdata[:, 0], fitdata[:, 1], 'r--', linewidth=2)

    # show another zoomed in fit, option to hide/show it might be interesting
    if lowerfitbound != -1 and upperfitbound != -1:

        if lowerfitbound == upperfitbound:
            lowerfitbound = fitdata[0, 0]
            upperfitbound = fitdata[len(fitdata)-1, 0]

            lb = p[1]-1.5
            ub = p[1]+1.5
        else:
            lb = lowerfitbound
            ub = upperfitbound

        title = filename + ' / AE = %.2f, ERR=%.5f\n' % (p[1], err) + 'relevant fitdata: (%.2f ; %.2f), ' % (
        lowerfitbound, upperfitbound) + 'fwhm = %.2f' % fwhm  # + '\n#iterations = %i' % iterations

        c_ax.set(title=title)

        mark_ae_data_in_plot(c_ax, p, lowerfitbound, upperfitbound, fwhm)

        # second subplot at the fit relevant area

        c_ax = axs.flatten()[1]

        cutdata = cutarray(data, lb, ub)

        if len(cutdata) > 0:
            c_ax.plot(cutdata[:,0], cutdata[:,1], 'o', markersize=3, markeredgecolor='b', markerfacecolor='none')
        c_ax.set(ylabel='Counts (1/s)')

        cut_fit_data = cutarray(fitdata, lb, ub)
        if len(cut_fit_data) > 0:
            c_ax.plot(cut_fit_data[:, 0], cut_fit_data[:, 1], 'r--', linewidth=2)

        mark_ae_data_in_plot(c_ax,p, lowerfitbound, upperfitbound, fwhm)

        #third subplot around +- 1 eV around AE

        c_ax = axs.flatten()[2]

        #show 1 eV around AE
        lb = p[1]-0.5
        ub = p[1]+0.5

        if lb > 0 and ub > lb:
          cutdata = cutarray(data, lb, ub)

          if len(cutdata)>0:
              c_ax.plot(cutdata[:,0], cutdata[:,1], 'o', markersize=3, markeredgecolor='b', markerfacecolor='none')


          c_ax.set(xlabel='Energy (eV)')

          cut_fit_data = cutarray(fitdata, lb, ub)
          if len(cut_fit_data)>0:
              c_ax.plot(cut_fit_data[:, 0], cut_fit_data[:, 1], 'r--', linewidth=2)

          mark_ae_data_in_plot(c_ax,p, lowerfitbound, upperfitbound, fwhm)

    return fig

def mark_ae_data_in_plot(ax, p, lowerfitbound, upperfitbound, fwhm):
    # mark AE
    ax.axvline(x=p[1], color='g', linestyle='-', linewidth='1')

    x_min = ax.viewLim.intervalx[0]
    x_max = ax.viewLim.intervalx[1]

    if lowerfitbound < x_min:
        lowerfitbound = x_min

    if upperfitbound > x_max:
        upperfitbound = x_max

    lb_fwhm = p[1] - 0.5*fwhm
    ub_fwhm = p[1] + 0.5*fwhm

    if lb_fwhm < x_min:
        lb_fwhm = x_min

    if ub_fwhm > x_max:
        ub_fwhm = x_max

    # mark relevant fit area

    if lowerfitbound != -1 and upperfitbound != -1:
        ax.axvspan(lowerfitbound, upperfitbound, facecolor='y', alpha=0.5)

    if fwhm > 0:
        ax.axvspan(lb_fwhm, ub_fwhm, facecolor='c', alpha=0.5)

def data_from_fit_and_parameters(data, fitfunc, p):
    minimumpoints = 20

    # check, whether there is a reasonable amount of points to plot:
    if len(data[:, 0]) > minimumpoints:
        # create empty float array in the size of data
        a = empty_like(data, dtype=float)
    else:
        a = empty((minimumpoints, 2), dtype=float)

    # created a problem, as x-points can slightly differ from the x-points in data, which is necessary for writing plot
    # to file and/or plotting the figures
    # Note: this also will make the check for minimumpoints not working for now, which might be left like this for now,
    # as data with less than 20 data-points seems rare
    # # create equidistant x points
    # a[:, 0] = linspace(data[:, 0].min(), data[:, 0].max(), len(a[:, 0]))

    a[:, 0] = data[:, 0]

    # calculate corresponding y values
    a[:, 1] = fitfunc(a[:, 0], p[0], p[1], p[2], p[3], p[4])

    return a

def eval_fit_function(fwhm):
    sigma = fwhm / (2 * sqrt(2 * np.log(2)))

    base_func = 'lambda x, p0, p1, p2, p3, p4: p2*(%s**p3)*gamma(p3+1)*exp(-1.0/(4.0*%s**2)*(p1-x)**2)*pbdv_fa(-(p3+1), (p1-x)/%s)+p0' % (
        sigma, sigma, sigma)

    return eval(base_func)

# !! weights has not been implemented, so far not needed
def do_the_fit(data, initial_parameters, fwhm, weights=[]):
    # res = optimize.curve_fit(fitfunc,
    #                         data[:, 0],
    #                         data[:, 1],
    #                         initial_parameters[:],
    #                         method='dogbox',
    #                         absolute_sigma=False,
    #                         bounds=([0, 0, -np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf, np.inf]),
    #                         #sigma=weights,
    #                         check_finite=True)

    # experimental: tighten sigma on the go

    base_func = eval_fit_function(fwhm)

    res = optimize.curve_fit(base_func,
                             data[:, 0],
                             data[:, 1],
                             initial_parameters[:],
                             method='trf',
                             absolute_sigma=False,
                             bounds=([0, 0, -np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf, np.inf]),
                             # sigma=weights,
                             check_finite=True)

    # return initial_parameters and stddev for the parameters
    return res[0], np.sqrt(np.diag(res[1]))[1]

def find_fwhm(starting_fwhm, ending_fwhm, iteration, n):

    a = starting_fwhm
    b = ending_fwhm
    c = n-1
    x = iteration

    return (b - a) / (c**2) * (x**2) + a #quadratic, nonconvs 74, 3, 0, 0, 0, 0, 0
    #return (b - a) / c * x + a #linear, nonconvs 63, 6, 0, 0, 0, 0, 0
    #return (a - b) / sqrt(c) * sqrt(-x + c) + b  # sqrt, nonconvs 50, 4, 0, 0, 0, 0, 0
    #return (a - b) * log(-x + c + 1)/log(c + 1) + b #logarithmic, nonconvs 33, 4, 0, 0, 0, 0, 0
    #return (a - b)*(1 - exp(x - c)) + b #exponential,nonconvs: 31, 3, 0, 0, 0, 0, 0
    # concave functions seem to work better in general

def find_cut_numbers(data_length, data_length_ev, minspan):

    n = 1
    cutpercent = 0

    steps = float64(data_length-1)

    if minspan > 0:
        # set a minimum of datapoints that are needed to do a decent fit (experimental value) (still not sure/experimenting)
        # not pretty but seems to be working well so far with tested data
        min_data_points = int64(25)

        energy_step = float64(data_length_ev)/steps

        min_data_points_energy = int64(minspan/energy_step)

        if min_data_points < min_data_points_energy:
            min_data_points = min_data_points_energy

        a = 2 #number of cut datapoints at last cut, experimental value

        n = int64((log(min_data_points)-log(data_length))/(log(min_data_points)-log(a+min_data_points)))

        # n is equal to the number of iterations that it will take to get the data points down to its minimum, when it is
        # cut by cutpercent*100 % every iteration.
        # the int() cast will floor the number, which will most likely lead to having more than min_data_points
        # data points left to fit.

        cutpercent = pow(min_data_points/data_length, 1.0/float64(n-1))
        #(n-1)st root, because the data won't be cut for the first iteration, and therefore only n-1 times

        if cutpercent > 1.0:
            cutpercent = 1.0
            n = 1

    return n, cutpercent

#Heart of the fit process
def fit_function_to_data(data, p, fwhm, minspan):
    # data has to be a numpy array
    # fits the function to data[:,0] (as x) and data[:,1] (as y) using the initial_parameters
    # returns an array of parameters of the same size as initial_parameters in case of success
    # returns None if data couldn't be fitted

    message = None

    #weights = set_fit_weights(data)

    # find number of iterations to take place where the data is cut down by x %
    n, cutpercent = find_cut_numbers(len(data), data[len(data)-1][0]-data[0][0], minspan)

    cutdata = data

    starting_fwhm = fwhm*3
    ending_fwhm = fwhm

    # fit another n times while closing in on the EA
    # The following conditions apply:
    # a) there have to be less than n iterations, to guarantee enough data points to fit AND
    # b) the energy span of the cut data has to be enough to guarantee a sensible fit

    iteration = 0
    stddev = -1
    nonconv = 0 #testing variable, not relevant for fitting
    while iteration < n:

        fwhm = find_fwhm(starting_fwhm, ending_fwhm, iteration, n)

        # try to fit
        try:
            p, stddev = do_the_fit(cutdata, p, fwhm)
            message = 'fit succeeded.'
        except Exception as error:
            if type(error) is ValueError:
              if error.args[0] == 'Residuals are not finite in the initial point.':
                  message = 'fit doesn\'t converge.'
            elif type(error) is RuntimeError:
                if error.args[0] == 'Optimal parameters not found: The maximum number of function evaluations is exceeded.':
                  message = 'fit doesn\'t converge'
            else:
                message = "unhandled error: " + error.args[0]

        # get the position of EA in the data
        ae_pos = find_ev_position(data, p[1])

        if ae_pos == -1:
            ae_pos = int(len(cutdata)/2)
            #set ae_pos to the center since it wasn't in the data
            #this will lead to cutting data nevertheless and could lead to a ae_pos in the center, if the fit never works.

            nonconv += 1
            #n += 1
            #increment n by one in order to ensure another run with a different fwhm to try to get a sensible ae
            #would probably bear problems with the accuracy in the play of minspan, cutpercent, and n

        # ae_pos has to be valid to cut the data.
        # At the end of the last run, it is not needed to cut again, and would even falsify data
        if iteration < n-1:
            # use that position, to cut the data down to cutpercent*100% of its size
            cutdata = cut_relatively_equal(data, ae_pos, cutpercent ** (iteration + 1))
        iteration += 1

    # bad luck, didn't converge
    #    if ier not in [1, 2, 3, 4]:
    #        return errmsg, None
    #    else:
    #return message, p, stddev, cutdata[0][0], cutdata[len(cutdata) - 1][0], fwhm, iteration
    return message, p, float(n)+ float(nonconv/100), cutdata[0][0], cutdata[len(cutdata) - 1][0], fwhm, iteration

def find_ev_position(data, ev):
    # find the position in data[] where the found ae lies
    e_pos = -1
    for energy in data[:, 0]:
        if energy < ev:
            e_pos += 1
        else:
            break

    if e_pos == len(data):
        e_pos = -1

    return e_pos

def cut_relatively_equal(cutdata, ae_pos, cutpercent):
    # cuts data down to cutpercent*100 % of its size, with EA in the middle
    # if however there are not enough data points on one of the sides to guarantee cutpercent*100 % of the data
    # more of the other side will be added.

    cdl = len(cutdata)
    buff_cutdata = []

    buff_cutdata.append(cutdata[ae_pos, :])

    # deviation from ae
    deviation = 1

    while len(buff_cutdata) / cdl <= cutpercent:

        current_pos = ae_pos - deviation

        if current_pos >= 0:
            buff_cutdata.insert(0, cutdata[current_pos, :])

        current_pos = ae_pos + deviation

        if current_pos < cdl:
            buff_cutdata.append(cutdata[current_pos, :])

        deviation += 1

    return np.array(buff_cutdata)

def cut_to_ev(cutdata, ae_pos, ev):
    # cut data down to +- ev around ae_pos
    # !! definitely improve list and numpy.array usage
    e_pos = 0
    ae = cutdata[ae_pos, 0]
    ret_cutdata = []

    for energy in cutdata[:, 0]:
        if abs(energy - ae) < ev:
            ret_cutdata.append(np.array(cutdata[e_pos, :]))
        e_pos += 1

    return np.array(ret_cutdata)

def difference_data_fit(data, fitfunc, p):
    fitdata = data_from_fit_and_parameters(data, fitfunc, p)
    newdata = []
    i = 0

    for datapoint in fitdata:
        newdata.append([data[i][0], data[i][1] - fitdata[i][1]])
        i += 1

    return np.array(newdata)

def set_fit_weights(data):
    # calculate the sqrt of the data values for the fit weights

    vecsqrt = vectorize(lambda x: x ** (1 / 2))
    vecabs = vectorize(abs)
    absolutedatavalues = vecabs((data[:, 1]))
    weights = vecsqrt(absolutedatavalues)

    # note that we set all zero values to 1, in order to have useful weights
    place(weights, weights == 0, 1)

    return weights

def cutarray(data, lowerlim=None, upperlim=None):
    # this function cuts an array and returns it
    # if lowerlim or upperlim are not defined, the maximum is assumed

    if lowerlim is None:
        lowerlim = data[:, 0].min()

    if upperlim is None:
        upperlim = data[:, 0].max()

    lowerlim = float64(lowerlim)
    upperlim = float64(upperlim)

    newdata = []

    for point in data:
        if (point[0] >= lowerlim) and (point[0] <= upperlim):
            newdata.append(point)

    data = array(newdata, dtype=float)

    return data

def create_plot_figure(data, filename, p, lowerfitbound, upperfitbound, fwhm, iterations, err):

    fig = plot_fit(data, eval_fit_function(fwhm), p, filename, err, lowerfitbound, upperfitbound, fwhm)  # plot fit

    return fig