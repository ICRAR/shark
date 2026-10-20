import common
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
import numpy as np
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import (FixedLocator, MaxNLocator,
                                             DictFormatter)
import matplotlib.pyplot as plt
import scipy.optimize as so
from matplotlib.image import NonUniformImage


def find_confidence_interval(x, pdf, confidence_level):
    return pdf[pdf > x].sum() - confidence_level

def density_contour(ax, xdata, ydata, nbins_x, nbins_y, cmap = 'grey'):
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

    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), density=True)
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

    return ax.contour(X, Y, Z, levels=levels, origin="lower",
                      norm=col.Normalize(vmin=0, vmax=0.025), **extra_args)


def plot_3d_dist(plt, outdir, dat, zout):
    fig = plt.figure(figsize=(5,4))
    ax = fig.add_subplot(111)
    
    x = dat[:,0]
    y = dat[:,1]
    ms = dat[:,2]
    lbol = dat[:,3]
    xmin, xmax, ymin, ymax = min(x), max(x), min(y), max(y)
    print(xmin, xmax, ymin, ymax)
    ytit = 'X/cMpc'
    xtit = 'Y/cMpc'
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(30, 30, 30, 30))

    #ax.hexbin(x, y, gridsize = (100, 100), cmap = 'Greys')
    #density_contour(ax, x, y, 200, 200, cmap = 'Greys')

    im = NonUniformImage(ax, interpolation='bilinear', cmap='magma_r')
    H, xedges, yedges = np.histogram2d(x, y, bins=(300,300))
    H = H.T
    xcenters = (xedges[:-1] + xedges[1:]) / 2
    ycenters = (yedges[:-1] + yedges[1:]) / 2
    im.set_data(xcenters, ycenters, H)
    ax.add_image(im)

    common.savefig(outdir, fig, "projected_density_nogalaxies_" + str(zout) + ".pdf")


    if(zout >= 3):
       fig = plt.figure(figsize=(5,4))
       ax = fig.add_subplot(111)
       common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(30, 30, 30, 30))
       im = NonUniformImage(ax, interpolation='bilinear', cmap='magma_r')
       H, xedges, yedges = np.histogram2d(x, y, bins=(300,300))
       H = H.T
       xcenters = (xedges[:-1] + xedges[1:]) / 2
       ycenters = (yedges[:-1] + yedges[1:]) / 2
       im.set_data(xcenters, ycenters, H)
       ax.add_image(im)
   
       ind = np.where((lbol > 43) & (lbol < 44.5))
       ax.plot(x[ind], y[ind], marker='o', markersize=1, color='darkgreen', fillstyle='full', alpha=0.5, linestyle='None')
       ind = np.where((lbol > 45.8) & (lbol < 50))
       print(lbol[ind])
       ax.plot(x[ind], y[ind], marker='*', markersize=7, color='black', fillstyle='full', alpha=0.5, linestyle='None')
   
       common.savefig(outdir, fig, "projected_density_withgalaxies_" + str(zout) + ".pdf")

def prepare_data(hdf5_data):

    (h0, _, mstard, mstarb, lbol, typeg, xg, yg, zg) = hdf5_data

    mstart = (mstard + mstarb)/h0

    ind = np.where(zg/h0 < 30)
    data = np.zeros(shape = (len(xg[ind]),4))
    data[:,0] = xg[ind]/h0
    data[:,1] = yg[ind]/h0
    data[:,2] = np.log10(mstart[ind])
    data[:,3] = np.log10(lbol[ind] * 1e40 / h0) #bolometric luminosity in erg/s

    return data

def main(modeldir, outdir, redshift_table, subvols, obsdir):

    zlist = [3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    plt = common.load_matplotlib()

    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'bolometric_luminosity_agn', 'type', 'position_x', 'position_y', 'position_z')}

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        data = prepare_data(hdf5_data)
        plot_3d_dist(plt, outdir, data, zlist[index])

if __name__ == '__main__':
    main(*common.parse_args())

