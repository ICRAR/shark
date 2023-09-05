#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2018
# Copyright by UWA (in the framework of the ICRAR)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

import numpy as np
import scipy.optimize as so

def wmedians_2sigma(x=None, y=None, xbins=None):

    nbins = len(xbins)
    #define size of bins, assuming bins are all equally spaced.
    dx = xbins[1] - xbins[0]
    result = np.zeros(shape = (3, nbins))

    for i in range (0,nbins):
        xlow = xbins[i]-dx/2.0
        xup  = xbins[i]+dx/2.0
        ind  = np.where((x > xlow) & (x< xup))
        if(len(x[ind]) > 9):

            obj_bin = len(x[ind])
            ybin    = y[ind]
            result[0, i] = np.median(ybin)
            #sort array on 1/y because we want it to sort from the smallest to the largest item, and the default of argsort is to order from the largest to the smallest.
            IDs = np.argsort(ybin,kind='quicksort')
            ID16th = int(np.floor(obj_bin*0.025))+1   #take the lower edge.
            ID84th = int(np.floor(obj_bin*0.975))-1   #take the upper edge.
            result[1, i] = np.abs(result[0, i] - ybin[IDs[ID16th]])
            result[2, i] = np.abs(ybin[IDs[ID84th]] - result[0, i])

    return result

def gpercentiles(x=None):

    result = np.zeros(shape = (3))

    ind = np.where(x != 0)
    if (len(x[ind]) > 0):    
       x = x[ind]
       result[0] = np.median(x)
       if(len(x) > 9):
          IDs = np.argsort(x,kind='quicksort')
          obj_bin = len(x)
          ID16th = int(np.floor(obj_bin*0.16))+1   #take the lower edge.
          ID84th = int(np.floor(obj_bin*0.84))-1   #take the upper edge.
          result[1] = np.abs(result[0] - x[IDs[ID16th]])
          result[2] = np.abs(x[IDs[ID84th]] - result[0])
       else:
          result[1] = np.abs(result[0] - min(x))
          result[2] = np.abs(max(x) - result[0])

    return result

def wmedians(x=None, y=None, xbins=None, low_numbers=False, nmin=10):

    nbins = len(xbins)
    #define size of bins, assuming bins are all equally spaced.
    dx = xbins[1] - xbins[0]
    result = np.zeros(shape = (3, nbins))

    for i in range (0,nbins):
        xlow = xbins[i]-dx/2.0
        xup  = xbins[i]+dx/2.0
        ind  = np.where((x > xlow) & (x< xup))
        if(len(x[ind]) >= nmin):

            obj_bin = len(x[ind])
            ybin    = y[ind]
            result[0, i] = np.median(ybin)
            #sort array on 1/y because we want it to sort from the smallest to the largest item, and the default of argsort is to order from the largest to the smallest.
            IDs = np.argsort(ybin,kind='quicksort')
            ID16th = int(np.floor(obj_bin*0.16))+1   #take the lower edge.
            ID84th = int(np.floor(obj_bin*0.84))-1   #take the upper edge.
            result[1, i] = np.abs(result[0, i] - ybin[IDs[ID16th]])
            result[2, i] = np.abs(ybin[IDs[ID84th]] - result[0, i])
        elif(low_numbers and len(x[ind]) > 0):
            ybin    = y[ind]
            result[0, i] = np.median(ybin)
            result[1, i] = np.abs(result[0, i] - np.min(y[ind]))
            result[2, i] = np.abs(np.max(y[ind]) - result[0, i])

    return result

def wmedians_2575(x=None, y=None, xbins=None, low_numbers=False, nmin=10):

    nbins = len(xbins)
    #define size of bins, assuming bins are all equally spaced.
    dx = xbins[1] - xbins[0]
    result = np.zeros(shape = (3, nbins))

    for i in range (0,nbins):
        xlow = xbins[i]-dx/2.0
        xup  = xbins[i]+dx/2.0
        ind  = np.where((x > xlow) & (x< xup))
        if(len(x[ind]) >= nmin):

            obj_bin = len(x[ind])
            ybin    = y[ind]
            result[0, i] = np.median(ybin)
            #sort array on 1/y because we want it to sort from the smallest to the largest item, and the default of argsort is to order from the largest to the smallest.
            IDs = np.argsort(ybin,kind='quicksort')
            ID16th = int(np.floor(obj_bin*0.25))+1   #take the lower edge.
            ID84th = int(np.floor(obj_bin*0.75))-1   #take the upper edge.
            result[1, i] = np.abs(result[0, i] - ybin[IDs[ID16th]])
            result[2, i] = np.abs(ybin[IDs[ID84th]] - result[0, i])
        elif(low_numbers and len(x[ind]) > 0):
            ybin    = y[ind]
            result[0, i] = np.median(ybin)
            result[1, i] = np.abs(result[0, i] - np.min(y[ind]))
            result[2, i] = np.abs(np.max(y[ind]) - result[0, i])

    return result


def stacking(x=None, y=None, xbins=None, low_numbers=False):

    nbins = len(xbins)
    #define size of bins, assuming bins are all equally spaced.
    dx = xbins[1] - xbins[0]
    result = np.zeros(shape = (2,nbins))

    for i in range (0,nbins):
        xlow = xbins[i]-dx/2.0
        xup  = xbins[i]+dx/2.0
        ind  = np.where((x > xlow) & (x<= xup))
        if(len(x[ind]) > 0):
            ybin    = y[ind]
            xbin    = x[ind]
            result[0,i] = np.log10(np.mean(10**xbin))
            result[1,i] = np.log10(np.mean(ybin))

    return result


def fractions(x=None, y=None, xbins=None, ythresh=None):

    nbins = len(xbins)
    #define size of bins, assuming bins are all equally spaced.
    dx = xbins[1]-xbins[0]
    result = np.zeros(shape = (nbins))

    for i in range (0,nbins):
        xlow = xbins[i]-dx/2.0
        xup  = xbins[i]+dx/2.0
        ind  = np.where((x > xlow) & (x< xup))
        if(len(x[ind]) > 9):
            ngalaxies = len(x[ind])
            ybin      = y[ind]
            above     = np.where(ybin > ythresh)
            nabove    = len(ybin[above])
            result[i] = (nabove+ 0.0)/(ngalaxies+0.0)
        else:
            result[i] = -1

    return result

def fractional_contribution(x=None, y=None, xbins=None):

    nbins = len(xbins)
    #define size of bins, assuming bins are all equally spaced.
    dx = xbins[1]-xbins[0]
    result = np.zeros(shape = (nbins))

    for i in range (0,nbins):
        xlow = xbins[i]-dx/2.0
        xup  = xbins[i]+dx/2.0
        ind  = np.where((x > xlow) & (x< xup))
        if(len(x[ind]) > 4):
            mtotal    = sum(x[ind])
            mbulge_tot= sum(x[ind]*y[ind])
            result[i] = mbulge_tot/mtotal
        else:
            result[i] = -1

    return result



def find_confidence_interval(x, pdf, confidence_level):
    return pdf[pdf > x].sum() - confidence_level



def density_contour(ax, xdata, ydata, nbins_x, nbins_y, cmap = 'viridis'):
    """ Create a density contour plot.
    Parameters
    ----------
    ax : matplotlib.Axes
        Plot the contour to this axis
    xdata : numpy.ndarray
    ydata : numpy.ndarray
    nbins_x : int
        Number of bins along x dimension
    nbins_y : int
        Number of bins along y dimension
    contour_kwargs : dict
        kwargs to be passed to pyplot.contour()
    """

    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), normed=True)
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

    pdf = (H*(x_bin_sizes*y_bin_sizes))

    low_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.01))
    twenty_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.2))
    thirty_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.3))
    forty_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.4))
    fifty_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.5))
    one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
    eighty_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.8))
    ninty_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.9))
    two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
    three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
    levels = [three_sigma, two_sigma, ninty_sigma, eighty_sigma, one_sigma, fifty_sigma, forty_sigma, thirty_sigma, twenty_sigma, low_sigma]

    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T

    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.colors as col

    # The viridis colormap is only available since mpl 1.5
    extra_args = {}
    if tuple(mpl.__version__.split('.')) >= ('1', '5'):
        extra_args['cmap'] = plt.get_cmap(cmap)

    return ax.contourf(X, Y, Z, levels=levels, origin="lower", alpha=0.75,
                      norm=col.Normalize(vmin=0, vmax=0.01), **extra_args)

def density_contour_reduced_col(ax, xdata, ydata, nbins_x, nbins_y, cmap = 'viridis'):
    """ Create a density contour plot.
    Parameters
    ----------
    ax : matplotlib.Axes
        Plot the contour to this axis
    xdata : numpy.ndarray
    ydata : numpy.ndarray
    nbins_x : int
        Number of bins along x dimension
    nbins_y : int
        Number of bins along y dimension
    contour_kwargs : dict
        kwargs to be passed to pyplot.contour()
    """

    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), normed=True)
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

    pdf = (H*(x_bin_sizes*y_bin_sizes))

    thirty_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.5))
    one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
    two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
    three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
    levels = [three_sigma, two_sigma, one_sigma, thirty_sigma]


    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T

    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.colors as col

    # The viridis colormap is only available since mpl 1.5
    extra_args = {}
    if tuple(mpl.__version__.split('.')) >= ('1', '5'):
        extra_args['cmap'] = plt.get_cmap(cmap)

    return ax.contourf(X, Y, Z, levels=levels, origin="lower", alpha=0.75,
                      norm=col.Normalize(vmin=0, vmax=0.01), **extra_args)


def density_contour_reduced(ax, xdata, ydata, nbins_x, nbins_y, cmap = 'viridis'):
    """ Create a density contour plot.
    Parameters
    ----------
    ax : matplotlib.Axes
        Plot the contour to this axis
    xdata : numpy.ndarray
    ydata : numpy.ndarray
    nbins_x : int
        Number of bins along x dimension
    nbins_y : int
        Number of bins along y dimension
    contour_kwargs : dict
        kwargs to be passed to pyplot.contour()
    """

    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), normed=True)
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

    pdf = (H*(x_bin_sizes*y_bin_sizes))

    thirty_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.5))
    one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
    two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
    three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
    levels = [three_sigma, two_sigma, one_sigma, thirty_sigma]

    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T

    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.colors as col

    return ax.contour(X, Y, Z, levels=levels, origin="lower", 
                      norm=col.Normalize(vmin=0, vmax=0.01), colors='Maroon')




def look_back_time(z, h=0.6751, omegam=0.3121, omegal=0.6879):

	"""Calculates the look back time of an array of redshifts
	Parameters
	---------
	z: array of redshifts
	h: hubble parameter
	omegam: omega matter
	omegal: omega lambda
	"""

	#define some constants:
	H0100=100.0
	KM2M=1.0e3
	GYR2S=3.15576e16
	MPC2M=3.0856775807e22
	H0100PGYR=H0100*KM2M*GYR2S/MPC2M

	#calculate the expansion parameters
	a = 1.0 / (1.0 + z)

	#The Hubble time for H_0=100km/s/Mpc
	Hubble_Time=1.0/H0100PGYR
	t0= Hubble_Time*(2/(3*h*np.sqrt(1-omegam)))*np.arcsinh(np.sqrt((1.0/omegam-1.0)*1.0)*1.0)
	t = Hubble_Time*(2/(3*h*np.sqrt(1-omegam)))*np.arcsinh(np.sqrt((1.0/omegam-1.0)*a)*a)

	return t0-t


def redshift(lbt, h=0.6751, omegam=0.3121, omegal=0.6879):

	"""Calculates the look back time of an array of redshifts
	Parameters
	---------
	z: array of redshifts
	h: hubble parameter
	omegam: omega matter
	omegal: omega lambda
	"""

	#define some constants:
	H0100=100.0
	KM2M=1.0e3
	GYR2S=3.15576e16
	MPC2M=3.0856775807e22
	H0100PGYR=H0100*KM2M*GYR2S/MPC2M

	#calculate the expansion parameters
	#a = 1.0 / (1.0 + z)

	#The Hubble time for H_0=100km/s/Mpc
	Hubble_Time=1.0/H0100PGYR
	t0= Hubble_Time*(2/(3*h*np.sqrt(1-omegam)))*np.arcsinh(np.sqrt((1.0/omegam-1.0)*1.0)*1.0)
	age = t0 - lbt

	a = pow(np.sinh(age/(Hubble_Time*(2/(3*h*np.sqrt(1-omegam))))) / np.sqrt(1.0/omegam-1.0), 2.0/3.0)
	
	z = 1.0 / a - 1.0
	for i in range (0,len(z)):
		z[i] = round(z[i], 2)

	return z


def bootstrap_error(x=None, y=None, xbins = None, iterations=50):

    nbins = len(xbins)
    len_subsample = 100 # same as used in observational data

    # initialise array, one row per mass bin
    medians = np.zeros(shape = [nbins,iterations])

    #define size of bins, assuming bins are all equally spaced.
    dx = xbins[1] - xbins[0]

    for i in range (0,nbins):
        xlow = xbins[i]-dx/2.0
        xup  = xbins[i]+dx/2.0
        ind  = np.where((x > xlow) & (x< xup))
        #if(len(x[ind]) > 9):
         #   ybin    = y[ind]
        ybin = y[ind]
        if len(ybin) < 2:
            ybin = [0,0] # sample function behaves weirdly when len(ybin) =1, allocate sd = 0 for low number cases
        for j in range(iterations):
             # calculate median size 30 times
             subsample = np.random.choice(ybin, size=len_subsample)
             medians[i,j] = np.median(subsample)

    # return stddev of sample medians for each mass bin
    return np.std(medians, axis=1)
