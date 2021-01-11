# version d6df05b15079004c7f95062fb8228580
from __future__ import division, print_function
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import make_interp_spline
from numpy.polynomial import polynomial

# these values are obtained from a fit (explained in https://github.com/philippeller/angular_acceptance/blob/master/Angular_acceptance.ipynb)
support_x = np.array([-1.  , -0.5 , -0.2 ,  0.35,  0.65,  0.95,  1.05])
components = np.array([[ 0.        , -0.00545327,  0.01655248, -0.1366884 , -0.07822519,
         0.51390016,  0.84308956],
       [-0.        , -0.45087175, -0.51799428,  0.24076603,  0.60613246,
         0.30925652, -0.08597765]])
mean = np.array([0.        , 0.14623077, 0.26576604, 0.4832101 , 0.63038745,
       0.57494938, 0.48991044])
n_components = 2

def ang(params, values):
    """
    New angular acceptance function

    params : list / array
        the parameters p0, p1, ...
    values : float, list, array
        the eta values to compute the angular acceptance for in (-1, 1)
    """
    # sanity check
    assert np.all(np.logical_and(np.greater_equal(values, -1), np.less_equal(values, 1))), 'values must be in range -1, 1'
    p = np.zeros(n_components)
    p[:len(params)] = params
    # inverse PCA transform
    transformed_params = np.dot(p, components) + mean
    # construct spline
    f = make_interp_spline(support_x, transformed_params, bc_type=([(1, 0.)], [(1, 0.)]))
    # make sure we're positive everywhere
    positive_f = lambda x : np.clip(f(x), 0., None)
    # normalize
    norm = quad(positive_f, -1, 1)
    out = positive_f(values)
    #normalize....why 0.68?
    out *= 0.68 / norm[0]
    return out

def angsens_poly(params):
    """
    Return standard IceCube format

    params : list/array
        the parameters p0, p1, ...
    """

    x = np.linspace(-1,1,101)
    sampled = ang(params, x)
    coeffs = np.zeros(12)
    coeffs[0] = np.max(sampled)
    coeffs[1:] = polynomial.polyfit(x, sampled, 10)
    return coeffs

from sys import argv
from os.path import join
args = argv[1:]

# get unified parameter input
p0 = args[0]
p1 = args[1]
p  = [p0, p1]

# get dima parameterization
coeffs = angsens_poly(p)

# write dima to as.dat
out = args[2]
f   = open(out, 'w')
for c in coeffs:
    f.write('%.5f\n' %(c))
f.close()
