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


def add_observations_to_plot(obsdir, fname, ax, marker, label, color='k', err_absolute=False):
    fname = '%s/Gas/%s' % (obsdir, fname)
    x, y, yerr_down, yerr_up = common.load_observation(obsdir, fname, (0, 1, 2, 3))
    common.errorbars(ax, x, y, yerr_down, yerr_up, color, marker, label=label, err_absolute=err_absolute)

def prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit):
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)
    xleg = xmax - 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)
    #ax.text(xleg, yleg, 'z=0')

def prepare_data(hdf5_data):

    bin_it = functools.partial(us.wmedians, xbins=xmf)

    # Unpack data
    (h0, _, typeg, mdisk, mbulge, _, _, mHI, mH2, mgas,
     mHI_bulge, mH2_bulge, mgas_bulge) = hdf5_data

    XH = 0.72
    h0log = np.log10(float(h0))

    n_typeg = len(typeg)
    morpho_type = np.zeros(shape = (n_typeg))
    morpho_type_stellar = np.zeros(shape = (n_typeg))

    # Early-type galaxies criterion of Khochfar et al. (2011).
    ind = np.where((mbulge + mgas_bulge)/(mdisk + mbulge + mgas + mgas_bulge) > 0.5)
    morpho_type[ind] = 1.0

    # Early-type galaxies criterion based on stellar mass alone.
    ind = np.where((mbulge)/(mdisk + mbulge) > 0.5)
    morpho_type_stellar[ind] = 1.0

    mh2_gals = np.zeros(shape = (2, n_typeg))
    mh1_gals = np.zeros(shape = (2, n_typeg))
    mgas_gals = np.zeros(shape = (2, n_typeg))

    mh2_gals_ltg = np.zeros(shape = (2, n_typeg))
    mh1_gals_ltg = np.zeros(shape = (2, n_typeg))
    mgas_gals_ltg = np.zeros(shape = (2, n_typeg))

    mh2_gals_etg = np.zeros(shape = (2, n_typeg))
    mh1_gals_etg = np.zeros(shape = (2, n_typeg))
    mgas_gals_etg = np.zeros(shape = (2, n_typeg))

    #Gas scaling relations based on morphological criterion calculated using total baryon mass of disk and bulge
    ind = np.where((mdisk + mbulge > 0) & (mgas + mgas_bulge > 0) & (mH2 + mH2_bulge > 0) & (morpho_type == 0))
    # Data we'll use later
    mass = mdisk[ind] + mbulge[ind]
    mgas_gals_ltg[0,ind] = mh1_gals_ltg[0,ind] = mh2_gals_ltg[0,ind] = np.log10(mass) - h0log
    mgas_gals_ltg[1,ind] = np.log10(XH * (mgas[ind] + mgas_bulge[ind]) / mass)
    mh1_gals_ltg[1,ind] = np.log10(XH * (mHI[ind] + mHI_bulge[ind]) / mass)
    mh2_gals_ltg[1,ind] = np.log10(XH * (mH2[ind] + mH2_bulge[ind]) / (mass))

    mgas_relation_ltg = bin_it(x=mgas_gals_ltg[0, ind], y=mgas_gals_ltg[1, ind])
    mh1_relation_ltg = bin_it(x=mh1_gals_ltg[0, ind], y=mh1_gals_ltg[1, ind])
    mh2_relation_ltg = bin_it(x=mh2_gals_ltg[0, ind], y=mh2_gals_ltg[1, ind])

    ind = np.where((mdisk + mbulge > 0) & (mgas + mgas_bulge > 0) & (mH2 + mH2_bulge > 0) & (morpho_type == 1))
    # Data we'll use later
    mass = mdisk[ind] + mbulge[ind]
    mgas_gals_etg[0,ind] = mh1_gals_etg[0,ind] = mh2_gals_etg[0,ind] = np.log10(mass) - h0log
    mgas_gals_etg[1,ind] = np.log10(XH * (mgas[ind] + mgas_bulge[ind]) / mass)
    mh1_gals_etg[1,ind] = np.log10(XH * (mHI[ind] + mHI_bulge[ind]) / mass)
    mh2_gals_etg[1,ind] = np.log10(XH * (mH2[ind] + mH2_bulge[ind]) / (mass))

    mgas_relation_etg = bin_it(x=mgas_gals_etg[0, ind], y=mgas_gals_etg[1, ind])
    mh1_relation_etg = bin_it(x=mh1_gals_etg[0, ind], y=mh1_gals_etg[1, ind])
    mh2_relation_etg = bin_it(x=mh2_gals_etg[0, ind], y=mh2_gals_etg[1, ind])

    #Gas scaling relations based on morphological criterion calculated using stellar mass of disk and bulge
    ind = np.where((mdisk + mbulge > 0) & (mgas + mgas_bulge > 0) & (mH2 + mH2_bulge > 0) & (morpho_type_stellar == 0))
    # Data we'll use later
    mass = mdisk[ind] + mbulge[ind]
    mgas_gals_ltg[0,ind] = mh1_gals_ltg[0,ind] = mh2_gals_ltg[0,ind] = np.log10(mass) - h0log
    mgas_gals_ltg[1,ind] = np.log10(XH * (mgas[ind] + mgas_bulge[ind]) / mass)
    mh1_gals_ltg[1,ind] = np.log10(XH * (mHI[ind] + mHI_bulge[ind]) / mass)
    mh2_gals_ltg[1,ind] = np.log10(XH * (mH2[ind] + mH2_bulge[ind]) / (mass))

    mgas_ms_relation_ltg = bin_it(x=mgas_gals_ltg[0, ind], y=mgas_gals_ltg[1, ind])
    mh1_ms_relation_ltg = bin_it(x=mh1_gals_ltg[0, ind], y=mh1_gals_ltg[1, ind])
    mh2_ms_relation_ltg = bin_it(x=mh2_gals_ltg[0, ind], y=mh2_gals_ltg[1, ind])

    ind = np.where((mdisk + mbulge > 0) & (mgas + mgas_bulge > 0) & (mH2 + mH2_bulge > 0) & (morpho_type_stellar == 1))
    # Data we'll use later
    mass = mdisk[ind] + mbulge[ind]
    mgas_gals_etg[0,ind] = mh1_gals_etg[0,ind] = mh2_gals_etg[0,ind] = np.log10(mass) - h0log
    mgas_gals_etg[1,ind] = np.log10(XH * (mgas[ind] + mgas_bulge[ind]) / mass)
    mh1_gals_etg[1,ind] = np.log10(XH * (mHI[ind] + mHI_bulge[ind]) / mass)
    mh2_gals_etg[1,ind] = np.log10(XH * (mH2[ind] + mH2_bulge[ind]) / (mass))

    mgas_ms_relation_etg = bin_it(x=mgas_gals_etg[0, ind], y=mgas_gals_etg[1, ind])
    mh1_ms_relation_etg = bin_it(x=mh1_gals_etg[0, ind], y=mh1_gals_etg[1, ind])
    mh2_ms_relation_etg = bin_it(x=mh2_gals_etg[0, ind], y=mh2_gals_etg[1, ind])




    # Constrains
    ind = np.where((mdisk + mbulge > 0) & (mgas + mgas_bulge > 0) & (mH2 + mH2_bulge > 0))

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
            mh2_relation, mh1_relation, mhr_relation, mhr_relation_cen, mhr_relation_sat, 
            mgas_relation_ltg, mh2_relation_ltg, mh1_relation_ltg, mgas_relation_etg, mh2_relation_etg, 
	    mh1_relation_etg, mgas_ms_relation_ltg, mh2_ms_relation_ltg, mh1_ms_relation_ltg, 
	    mgas_ms_relation_etg, mh2_ms_relation_etg, mh1_ms_relation_etg)

def plot_cold_gas_fraction(plt, output_dir, obs_dir, mgas_relation, mgas_relation_cen, mgas_relation_sat):

    ###################################
    #   Plots global mass densities
    fig = plt.figure(figsize=(5,4.5))

    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm log_{10}(M_{\\rm cold}/M_{\\star})$"

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

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
    add_observations_to_plot(obs_dir, 'NeutralGasRatio_NonDetEQUpperLimits.dat', ax, 'v', "GASS+COLDGASS", color='grey')
    add_observations_to_plot(obs_dir, 'NeutralGasRatio_NonDetEQZero.dat', ax, '^', "GASS+COLDGASS", color='grey')

    common.prepare_legend(ax, ['k','b','r','grey','grey'])
    common.savefig(output_dir, fig, "cold_gas_fraction.pdf")


def plot_molecular_gas_fraction(plt, output_dir, obs_dir, mgas_gals, mgas_relation, mh1_gals, mh1_relation, mh2_gals, mh2_relation, 
    mgas_relation_ltg, mh2_relation_ltg, mh1_relation_ltg, mgas_relation_etg, mh2_relation_etg, mh1_relation_etg, 
    mgas_ms_relation_ltg, mh2_ms_relation_ltg, mh1_ms_relation_ltg, mgas_ms_relation_etg, mh2_ms_relation_etg, 
    mh1_ms_relation_etg):

    xmin, xmax, ymin, ymax = 9, 12, -3, 1
    fig = plt.figure(figsize=(11,11))

    # First subplot
    ax = fig.add_subplot(321)
    plt.subplots_adjust(left=0.15)

    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm log_{10}(M_{\\rm HI+H_2}/M_{\\star})$"
    prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)

    #Predicted relation for all galaxies
    ind = np.where((mgas_gals[0,:] > 0) & (mgas_gals[1,:] != 0) )
    xdata = mgas_gals[0,ind]
    ydata = mgas_gals[1,ind]
    us.density_contour(ax, xdata[0], ydata[0], 30, 30) #, **contour_kwargs)

    def plot_mrelation(mrelation, color, label=None, linestyle=None):
        ind = np.where(mrelation[0,:] != 0)
        xplot = xmf[ind]
        yplot = mrelation[0,ind]
        linestyle = linestyle or ''
        ax.plot(xplot,yplot[0], color=color, label=label, linestyle=linestyle)

    def plot_mrelation_fill(mrelation, color, colorfill, label=None, linestyle=None):
        ind = np.where(mrelation[0,:] != 0)
        xplot = xmf[ind]
        yplot = mrelation[0,ind]
        errdn = mrelation[1,ind]
        errup = mrelation[2,ind]

        ax.plot(xplot,yplot[0], color=color, label=label, linestyle=linestyle)
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=colorfill, alpha=0.2,interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=colorfill, alpha=0.2,interpolate=True)


    plot_mrelation(mgas_relation, 'k', label="Shark all galaxies")

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    add_observations_to_plot(obs_dir, 'NeutralGasRatio_NonDetEQUpperLimits.dat', ax, 'v', "xCOLDGAS+xGASS")
    add_observations_to_plot(obs_dir, 'NeutralGasRatio_NonDetEQZero.dat', ax, '^', "xCOLDGAS+xGASS")

    common.prepare_legend(ax, ['k','k','k'])

    # Second subplot
    ax = fig.add_subplot(322)
    xmin, xmax, ymin, ymax = 9, 12, -4.5, 1

    prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)

    plot_mrelation_fill(mgas_relation_ltg, 'b', 'b',label="Shark LTGs $(\\rm B/T)_{\\rm bar}$",linestyle='solid')
    plot_mrelation_fill(mgas_relation_etg, 'r', 'r',label="Shark ETGs $(\\rm B/T)_{\\rm bar}$",linestyle='solid')
    plot_mrelation(mgas_ms_relation_ltg, 'b',label="Shark LTGs $(\\rm B/T)_{\star}$",linestyle='dotted')
    plot_mrelation(mgas_ms_relation_etg, 'r',label="Shark ETGs $(\\rm B/T)_{\star}$",linestyle='dotted')

    # Legend
    common.prepare_legend(ax, ['b','r','b','r','k'],loc=3)

    # Third subplot
    ax = fig.add_subplot(323)
    plt.subplots_adjust(left=0.15)

    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm log_{10}(M_{\\rm HI}/M_{\\star})$"
    xmin, xmax, ymin, ymax = 9, 12, -3, 1
    prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)

    #Predicted relation
    ind = np.where((mh1_gals[0,:] > 0) & (mh1_gals[1,:] != 0) )
    xdata = mh1_gals[0,ind]
    ydata = mh1_gals[1,ind]
    us.density_contour(ax, xdata[0], ydata[0], 30, 30) #, **contour_kwargs)
    plot_mrelation(mh1_relation, 'k')

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    add_observations_to_plot(obs_dir, 'HIGasRatio_NonDetEQUpperLimits.dat', ax, 'v', "xGASS")
    add_observations_to_plot(obs_dir, 'HIGasRatio_NonDetEQZero.dat', ax, '^', "xGASS")

    x, y, yerr_down, yerr_up = common.load_observation(obs_dir, 'Gas/Parkash18.dat', (0, 1, 2, 3))
    ax.errorbar(x,y-x,yerr=[(y-x) - (yerr_down-x),(yerr_up-x) - (y-x)], ls='None', mfc='r', fillstyle='full', ecolor = 'r', mec='r',marker='s',markersize=7, label="Parkash+18")

    m, mrat, merr = common.load_observation(obs_dir, 'Gas/RHI-Mstars_Brown15.dat', [0,1,2])
    errdn = np.log10(mrat) - np.log10(mrat - merr) 
    errup = np.log10(mrat + merr) - np.log10(mrat) 
    ax.errorbar(m,np.log10(mrat),yerr=[errdn,errup], ls='None', mfc='Salmon', fillstyle='full', ecolor = 'Salmon', mec='Salmon',marker='o',markersize=7, label="Brown+15")

    # Legend
    common.prepare_legend(ax, ['k','k','r','Salmon'], loc=1)

    # Fourth subplot
    ax = fig.add_subplot(324)
    xmin, xmax, ymin, ymax = 9, 12, -4.5, 1
    prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)

    plot_mrelation_fill(mh1_relation_ltg, 'b', 'b',linestyle='solid')
    plot_mrelation_fill(mh1_relation_etg, 'r', 'r',linestyle='solid')
    plot_mrelation(mh1_ms_relation_ltg,  'b',linestyle='dotted')
    plot_mrelation(mh1_ms_relation_etg,  'r',linestyle='dotted')

    add_observations_to_plot(obs_dir, 'RHI-Mstars_Callette18-LTGs.dat', ax, 's', "Calette+18 LTGs", color='grey', err_absolute=True)
    add_observations_to_plot(obs_dir, 'RHI-Mstars_Callette18-ETGs.dat', ax, 'o', "Calette+18 ETGs", color='grey', err_absolute=True)

    # Legend
    common.prepare_legend(ax, ['grey','grey','grey'],loc=1)

    # Fifth subplot
    ax = fig.add_subplot(325)
    plt.subplots_adjust(left=0.15)

    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm log_{10}(M_{\\rm H_2}/M_{\\star})$"
    xmin, xmax, ymin, ymax = 9, 12, -3, 1
    prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)

    #Predicted relation
    ind = np.where((mh2_gals[0,:] > 0) & (mh2_gals[1,:] != 0) )
    xdata = mh2_gals[0,ind]
    ydata = mh2_gals[1,ind]
    us.density_contour(ax, xdata[0], ydata[0], 30, 30) #, **contour_kwargs)
    plot_mrelation(mh2_relation, 'k')

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    add_observations_to_plot(obs_dir, 'MolecularGasRatio_NonDetEQUpperLimits.dat', ax, 'v', "xCOLDGASS")
    add_observations_to_plot(obs_dir, 'MolecularGasRatio_NonDetEQZero.dat', ax, '^', "xCOLDGASS")

    common.prepare_legend(ax, ['k','k','k'], loc = 1)

    # Fourth subplot
    ax = fig.add_subplot(326)
    xmin, xmax, ymin, ymax = 9, 12, -4.5, 1
    prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)

    plot_mrelation_fill(mh2_relation_ltg, 'b', 'b',linestyle='solid')
    plot_mrelation_fill(mh2_relation_etg, 'r', 'r',linestyle='solid')
    plot_mrelation(mh2_ms_relation_ltg, 'b',linestyle='dotted')
    plot_mrelation(mh2_ms_relation_etg, 'r',linestyle='dotted')

    add_observations_to_plot(obs_dir, 'RH2-Mstars_Callette18-LTGs.dat', ax, 's', "Calette+18 LTGs",color='grey', err_absolute=True)
    add_observations_to_plot(obs_dir, 'RH2-Mstars_Callette18-ETGs.dat', ax, 'o', "Calette+18 ETGs",color='grey', err_absolute=True)

    # Legend
    common.prepare_legend(ax, ['grey','grey','grey'],loc=1)

    common.savefig(output_dir, fig, "molecular_gas_fraction.pdf")

def plot_h1h2_gas_fraction(plt, output_dir, mhr_relation, mhr_relation_cen, mhr_relation_sat):

    fig = plt.figure(figsize=(5,4.5))

    ax = fig.add_subplot(111)
    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
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


def main(model_dir, output_dir, subvols, obs_dir, snapshot):

    plt = common.load_matplotlib()
    fields = {'galaxies': ('type', 'mstars_disk', 'mstars_bulge',
                           'rstar_disk', 'm_bh', 'matom_disk', 'mmol_disk', 'mgas_disk',
                           'matom_bulge', 'mmol_bulge', 'mgas_bulge')}
    hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)

    (mgas_relation, mgas_relation_cen, mgas_relation_sat,
     mh2_gals, mh1_gals, mgas_gals,
     mh2_relation, mh1_relation, mhr_relation, mhr_relation_cen, mhr_relation_sat,
     mgas_relation_ltg, mh2_relation_ltg, mh1_relation_ltg,
     mgas_relation_etg, mh2_relation_etg, mh1_relation_etg,
     mgas_ms_relation_ltg, mh2_ms_relation_ltg, mh1_ms_relation_ltg,
     mgas_ms_relation_etg, mh2_ms_relation_etg, mh1_ms_relation_etg) = prepare_data(hdf5_data)

    plot_cold_gas_fraction(plt, output_dir, obs_dir, mgas_relation, mgas_relation_cen, mgas_relation_sat)

    plot_molecular_gas_fraction(plt, output_dir, obs_dir, mgas_gals, mgas_relation, mh1_gals, mh1_relation, mh2_gals, mh2_relation, mgas_relation_ltg, 
	mh2_relation_ltg, mh1_relation_ltg, mgas_relation_etg, mh2_relation_etg, mh1_relation_etg, mgas_ms_relation_ltg, mh2_ms_relation_ltg, 
	mh1_ms_relation_ltg, mgas_ms_relation_etg, mh2_ms_relation_etg, mh1_ms_relation_etg)

    plot_h1h2_gas_fraction(plt, output_dir, mhr_relation, mhr_relation_cen, mhr_relation_sat)

if __name__ == '__main__':
    main(*common.parse_args())
