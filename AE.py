#!/usr/bin/python

from numpy import *
from scipy import *

from scipy.special import *
from scipy import optimize

import matplotlib.pyplot as plt

import fitlib as fl

sigma = 1.0
alpha = 1.27

#retrieve function for Appearance Energy with zero fixed parameters (hence the '')
ae_func = fl.get_AE_func(sigma, '')

#fixed alpha
#ae_func = lambda p, x: p[0] + p[1]*exp(-1.0/(4.0*sigma**2)*(p[2]-x)**2)*pbdv_fa(-(alpha+1), (p[2]-x)/sigma)

data = fl.readfile('AE_data/AE_HMX_H2O.txt')

#t1 = arange(1, 5, 0.1)

#y1 = zeros(t1.shape)
#y2 = zeros(t1.shape)

#j = 0

#for i in t1:
#    y1[j] = 10*pbdv_fa(i, i+3) + rand(1,1)*0.1
#    j += 1

#print t1.shape, y1.shape

fig1 = plt.figure()

fl.plotES(plt,data,'bla')

p0 = [0]*4
p0[0] = 50
p0[1] = 1
p0[2] = 1
p0[3] = 22

#fitfunc = lambda p, x: p[0]*pbdv_fa(x, x+3)
errfunc = lambda p, x, y: ae_func(p, x) - y

p1, success = optimize.leastsq(errfunc, p0[:], args=(data[:,0],data[:,1]))

print p1

plt.plot(data[:,0], ae_func(p1, data[:,0]), 'r--', linewidth=3)

plt.show()



