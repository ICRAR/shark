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
import functools

import numpy as np

import common
import utilities_statistics as us

##################################
# Constants
mlow = 6.0
mupp = 12.0
dm = 0.2
mbins = np.arange(mlow, mupp, dm)
xmf = mbins + dm/2.0


def add_observations_to_plot(obsdir, fname, ax, marker, label, color='k'):
    fname = '%s/Gas/%s' % (obsdir, fname)
    x, y, yerr_up, yerr_down = common.load_observation(obsdir, fname, (0, 1, 2, 3))
    common.errorbars(ax, x, y, yerr_down, yerr_up, color, marker, label=label, err_absolute=False)

def prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit):
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)
    xleg = xmax - 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)
    ax.text(xleg, yleg, 'z=0')

def prepare_data(hdf5_data):

    bin_it = functools.partial(us.wmedians, xbins=xmf)

    # Unpack data
    (h0, _, typeg, mdisk, mbulge, _, _, mHI, mH2, mgas,
     mHI_bulge, mH2_bulge, mgas_bulge) = hdf5_data

    XH = 0.72
    h0log = np.log10(float(h0))

    n_typeg = len(typeg)
    mh2_gals = np.zeros(shape = (2, n_typeg))
    mh1_gals = np.zeros(shape = (2, n_typeg))
    mgas_gals = np.zeros(shape = (2, n_typeg))

    # Constrains
    ind = np.where((mdisk + mbulge > 9) & (mgas + mgas_bulge > 0) & (mH2 + mH2_bulge > 0))

    # Data we'll use later
    mass = mdisk[ind] + mbulge[ind]
    mgas_gals[0,ind] = mh1_gals[0,ind] = mh2_gals[0,ind] = np.log10(mass) - h0log
    mgas_gals[1,ind] = np.log10(XH * (mgas[ind] + mgas_bulge[ind]) / mass)
    mh1_gals[1,ind] = np.log10(XH * (mHI[ind] + mHI_bulge[ind]) / mass)
    mh2_gals[1,ind] = np.log10(XH * (mH2[ind] + mH2_bulge[ind]) / (mass))

    # Binned relations
    mgas_relation = bin_it(x=mgas_gals[0, ind], y=mgas_gals[1, ind])
    mh1_relation = bin_it(x=mh1_gals[0, ind], y=mh1_gals[1, ind])
    mh2_relation = bin_it(x=mh2_gals[0, ind], y=mh2_gals[1, ind])
    mhr_relation = bin_it(x=np.log10(mdisk[ind]+mbulge[ind]) - h0log,
                        y=np.log10((mH2[ind] + mH2_bulge[ind]) / (mHI[ind] + mHI_bulge[ind])))

    ind = np.where((mdisk+mbulge > 0) & (typeg == 0))
    mass_central = mdisk[ind] + mbulge[ind]
    mgas_relation_cen = bin_it(x=np.log10(mass_central) - h0log,
                               y=np.log10(XH * (mgas[ind] + mgas_bulge[ind]) / (mdisk[ind] + mbulge[ind] + XH*mgas[ind] + XH*mgas_bulge[ind])))
    mhr_relation_cen = bin_it(x=np.log10(mass_central) - h0log,
                              y=np.log10((mH2[ind] + mH2_bulge[ind]) / (mHI[ind] + mHI_bulge[ind])))

    ind = np.where((mdisk+mbulge > 0) & (typeg > 0))
    mass_sat = np.log10(mdisk[ind] + mbulge[ind]) - h0log
    mgas_relation_sat = bin_it(x=mass_sat,
                               y=np.log10(XH * (mgas[ind] + mgas_bulge[ind]) / (mdisk[ind] + mbulge[ind] + XH*mgas[ind] + XH*mgas_bulge[ind])))
    mhr_relation_sat = bin_it(x=mass_sat,
                              y=np.log10((mH2[ind] + mH2_bulge[ind]) / (mHI[ind] + mHI_bulge[ind])))

    return (mgas_relation, mgas_relation_cen, mgas_relation_sat,
            mh2_gals, mh1_gals, mgas_gals,
            mh2_relation, mh1_relation, mhr_relation, mhr_relation_cen, mhr_relation_sat)

def plot_cold_gas_fraction(plt, output_dir, obs_dir, mgas_relation, mgas_relation_cen, mgas_relation_sat):

    ###################################
    #   Plots global mass densities
    fig = plt.figure(figsize=(5,5))

    xtit="$\\rm log_{10} (\\rm M_{\\rm star}/M_{\odot})$"
    ytit="$\\rm log_{10}(M_{\\rm cold}/M_{\\rm star})$"

    ax = fig.add_subplot(111)
    prepare_ax(ax, 8, 12, -3, 0.1, xtit, ytit)

    #Predicted SMHM
    ind = np.where(mgas_relation[0,:] != 0)
    xplot = xmf[ind]
    yplot = mgas_relation[0,ind]
    errdn = mgas_relation[1,ind]
    errup = mgas_relation[2,ind]

    ax.errorbar(xplot,yplot[0],color='k', label="all galaxies")
    ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'k', mec='k',marker='+',markersize=2)

    ind = np.where(mgas_relation_cen[0,:] != 0)
    xplot = xmf[ind]
    yplot = mgas_relation_cen[0,ind]
    ax.errorbar(xplot,yplot[0],color='b',linestyle='dotted', label="centrals")

    ind = np.where(mgas_relation_sat[0,:] != 0)
    xplot = xmf[ind]
    yplot = mgas_relation_sat[0,ind]
    ax.errorbar(xplot,yplot[0],color='r',linestyle='dashed', label="satelites")

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    add_observations_to_plot(obs_dir, 'NeutralGasFraction_NonDetEQUpperLimits.dat', ax, 'v', "GASS+COLDGASS", color='grey')
    add_observations_to_plot(obs_dir, 'NeutralGasFraction_NonDetEQZero.dat', ax, '^', "GASS+COLDGASS", color='grey')

    common.prepare_legend(ax, ['k','b','r','grey','grey'])
    common.savefig(output_dir, fig, "cold_gas_fraction.pdf")


def plot_molecular_gas_fraction(plt, output_dir, obs_dir, mgas_gals, mgas_relation, mh1_gals, mh1_relation, mh2_gals, mh2_relation):

    xmin, xmax, ymin, ymax = 9, 12, -3, 1
    fig = plt.figure(figsize=(5,12.5))

    # First subplot
    ax = fig.add_subplot(311)
    xtit="$\\rm log_{10} (\\rm M_{\\rm star}/M_{\odot})$"
    ytit="$\\rm log_{10}(M_{\\rm HI+H_2}/M_{\\rm star})$"
    prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)

    #Predicted relation
    ind = np.where((mgas_gals[0,:] > 0) & (mgas_gals[1,:] != 0) )
    xdata = mgas_gals[0,ind]
    ydata = mgas_gals[1,ind]
    us.density_contour(xdata[0], ydata[0], 30, 30, ax=ax) #, **contour_kwargs)
    ind = np.where(mgas_relation[0,:] != 0)
    xplot = xmf[ind]
    yplot = mgas_relation[0,ind]
    ax.plot(xplot,yplot[0],color='k', label="SHArk all galaxies")

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    add_observations_to_plot(obs_dir, 'NeutralGasRatio_NonDetEQUpperLimits.dat', ax, 'v', "COLDGAS+GASS")
    add_observations_to_plot(obs_dir, 'NeutralGasRatio_NonDetEQZero.dat', ax, '^', "COLDGAS+GASS")

    common.prepare_legend(ax, ['k','k','k'])

    # Second subplot
    ax = fig.add_subplot(312)
    xtit="$\\rm log_{10} (\\rm M_{\\rm star}/M_{\odot})$"
    ytit="$\\rm log_{10}(M_{\\rm HI}/M_{\\rm star})$"
    prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)

    #Predicted relation
    ind = np.where((mh1_gals[0,:] > 0) & (mh1_gals[1,:] != 0) )
    xdata = mh1_gals[0,ind]
    ydata = mh1_gals[1,ind]
    us.density_contour(xdata[0], ydata[0], 30, 30, ax=ax) #, **contour_kwargs)
    ind = np.where(mh1_relation[0,:] != 0)
    xplot = xmf[ind]
    yplot = mh1_relation[0,ind]
    ax.plot(xplot,yplot[0],color='k', label="SHArk all galaxies")

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    add_observations_to_plot(obs_dir, 'HIGasRatio_NonDetEQUpperLimits.dat', ax, 'v', "GASS")
    add_observations_to_plot(obs_dir, 'HIGasRatio_NonDetEQZero.dat', ax, '^', "GASS")

    # Legend
    common.prepare_legend(ax, ['k','k','k'])

    # Third subplot
    ax = fig.add_subplot(313)
    xtit="$\\rm log_{10} (\\rm M_{\\rm star}/M_{\odot})$"
    ytit="$\\rm log_{10}(M_{\\rm H_2}/M_{\\rm star})$"
    prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)

    #Predicted relation
    ind = np.where((mh2_gals[0,:] > 0) & (mh2_gals[1,:] != 0) )
    xdata = mh2_gals[0,ind]
    ydata = mh2_gals[1,ind]
    us.density_contour(xdata[0], ydata[0], 30, 30, ax=ax) #, **contour_kwargs)
    ind = np.where(mh2_relation[0,:] != 0)
    xplot = xmf[ind]
    yplot = mh2_relation[0,ind]
    ax.plot(xplot,yplot[0],color='k', label="SHArk all galaxies")

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    add_observations_to_plot(obs_dir, 'MolecularGasRatio_NonDetEQUpperLimits.dat', ax, 'v', "COLDGASS")
    add_observations_to_plot(obs_dir, 'MolecularGasRatio_NonDetEQZero.dat', ax, '^', "COLDGASS")

    common.prepare_legend(ax, ['k','k','k'])
    common.savefig(output_dir, fig, "molecular_gas_fraction.pdf")

def plot_h1h2_gas_fraction(plt, output_dir, mhr_relation, mhr_relation_cen, mhr_relation_sat):

    fig = plt.figure(figsize=(5,5))

    ax = fig.add_subplot(111)
    xtit="$\\rm log_{10} (\\rm M_{\\rm star}/M_{\odot})$"
    ytit="$\\rm log_{10}(M_{\\rm H_2}/M_{\\rm HI})$"
    prepare_ax(ax, 8, 12, -3, 1.0, xtit, ytit)

    # Predicted SMHM
    ind = np.where(mhr_relation[0,:] != 0)
    xplot = xmf[ind]
    yplot = mhr_relation[0,ind]
    errdn = mhr_relation[1,ind]
    errup = mhr_relation[2,ind]

    ax.errorbar(xplot,yplot[0],color='k', label="all galaxies")
    ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'k', mec='k',marker='+',markersize=2)

    ind = np.where(mhr_relation_cen[0,:] != 0)
    xplot = xmf[ind]
    yplot = mhr_relation_cen[0,ind]
    ax.errorbar(xplot,yplot[0],color='b',linestyle='dotted', label="centrals")

    ind = np.where(mhr_relation_sat[0,:] != 0)
    xplot = xmf[ind]
    yplot = mhr_relation_sat[0,ind]
    ax.errorbar(xplot,yplot[0],color='r',linestyle='dashed', label="satelites")

    common.prepare_legend(ax, ['k','b','r','grey','grey'])
    common.savefig(output_dir, fig, "HIH2_gas_fraction.pdf")


def main():

    plt = common.load_matplotlib()

    model_dir, output_dir, obs_dir, snapshot = common.parse_args()

    fields = {'Galaxies': ('type', 'mstars_disk', 'mstars_bulge',
                           'rdisk', 'mBH', 'matom_disk', 'mmol_disk', 'mgas_disk',
                           'matom_bulge', 'mmol_bulge', 'mgas_bulge')}
    hdf5_data = common.read_data(model_dir, snapshot, fields)

    (mgas_relation, mgas_relation_cen, mgas_relation_sat,
     mh2_gals, mh1_gals, mgas_gals,
     mh2_relation, mh1_relation, mhr_relation, mhr_relation_cen, mhr_relation_sat) = prepare_data(hdf5_data)

    plot_cold_gas_fraction(plt, output_dir, obs_dir, mgas_relation, mgas_relation_cen, mgas_relation_sat)
    plot_molecular_gas_fraction(plt, output_dir, obs_dir, mgas_gals, mgas_relation, mh1_gals, mh1_relation, mh2_gals, mh2_relation)
    plot_h1h2_gas_fraction(plt, output_dir, mhr_relation, mhr_relation_cen, mhr_relation_sat)

if __name__ == '__main__':
    main()