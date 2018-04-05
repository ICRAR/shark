#
#    ICRAR - International Centre for Radio Astronomy Research
#    (c) UWA - The University of Western Australia, 2018
#    Copyright by UWA (in the framework of the ICRAR)
#    All rights reserved
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#    MA 02111-1307  USA
#

import numpy as np
import scipy.optimize as so

def wmedians(x=None, y=None, xbins=None):

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
            ID16th = int(np.floor(obj_bin*0.16))+1   #take the lower edge.
            ID84th = int(np.floor(obj_bin*0.84))-1   #take the upper edge.
            result[1, i] = np.abs(result[0, i] - ybin[IDs[ID16th]])
            result[2, i] = np.abs(ybin[IDs[ID84th]] - result[0, i])
        elif(len(x[ind]) > 0 & len(x[ind]) < 9):
            result[0, i] = np.median(y[ind])
            result[1, i] = np.abs(result[0, i]-np.min(y[ind]))
            result[2, i] = np.abs(np.max(y[ind])-result[0, i])

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
        if(len(x[ind]) > 4):
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



def density_contour(xdata, ydata, nbins_x, nbins_y, ax=None):
    """ Create a density contour plot.
    Parameters
    ----------
    xdata : numpy.ndarray
    ydata : numpy.ndarray
    nbins_x : int
        Number of bins along x dimension
    nbins_y : int
        Number of bins along y dimension
    ax : matplotlib.Axes (optional)
        If supplied, plot the contour to this axis. Otherwise, open a new figure
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

    import matplotlib.pyplot as plt
    import matplotlib.colors as col

    if ax == None:
        contour = plt.contour(X, Y, Z, levels=levels, origin="lower")
    else:
        contour = ax.contourf(X, Y, Z, levels=levels, origin="lower", alpha=0.75,norm=col.Normalize(vmin=0,vmax=0.01),cmap=plt.get_cmap('viridis'))

    return contour