#!/usr/bin/python

from numpy import *
from scipy import *

from scipy.special import *
from scipy import optimize

import matplotlib.pyplot as plt

import fitlib as fl

def pbdv_fa(x,y):
    a, b = pbdv(x, y)
    return a

ae_func = lambda p, x: p[0] + p[1]*exp(-1.0/(4.0*sigma**2)*(p[2]-x)**2)*pbdv_fa(-(alpha+1), (p[2]-x)/sigma)
errfunc = lambda p, x, y: ae_func(p, x) - y

data = fl.readfile('AE_data/AE_HMX_H2O.txt')

sigma = 0.9

alpha = 0.0

p0 = [0]*3
p0[0] = 50
p0[1] = 10
p0[2] = 22

i = 0.0
j = 0

runs = 300

results = []

fig1 = plt.figure()


while i < runs:
    i = i + 1
    alpha = i/100
    p1, success = optimize.leastsq(errfunc, p0[:], args=(data[:,0],data[:,1]))

    results.append([alpha, p1[2]])
    j = j + 1
    print [alpha, p1[2]]

results = array(results,dtype = float)

plt.plot(results[:,0], results[:,1], 'bo')

plt.show()


