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
from cmath import sqrt
"""Global plots"""

import math

import numpy as np

import common


##################################
# Constants
GyrToYr = 1e9
Zsun = 0.0127
minmass = 1.0
Omegab = 0.0491
G    = 4.299e-9 #Gravity constant in units of (km/s)^2 * Mpc/Msun
rho_crit = 3.0 * pow(100.0,2.0) / 8 / math.pi / G #in units of h^2*Msun/Mpc^3
sbar = rho_crit * Omegab
OmegaM = 0.3121
OmegaL = 0.6879

def prepare_data(hdf5_data, redshifts):

    (h0, volh, _, mHI, mH2, mcold, mcold_metals, mhot, meje, mstar,
     mstar_burst, mBH, sfrdisk, sfrburst, mDM, mcold_halo) = hdf5_data

    sfrall = sfrdisk + sfrburst

    #Add up cold halo component to hot gas.
    mhot = mhot + mcold_halo

    massbar = mcold+mhot+meje+mstar+mBH
    sfr  = sfrall / volh / GyrToYr
    sfrd = sfrdisk  / volh / GyrToYr
    sfrb = sfrburst / volh / GyrToYr

    ind = np.where(mcold <= 0.0)
    mcold[ind] = minmass
    ind = np.where(mstar <= 0.0)
    mstar[ind] = minmass
    ind = np.where(mhot <= 0.0)
    mhot[ind] = minmass
    ind = np.where(meje <= 0.0)
    meje[ind] = minmass
    ind = np.where(mHI <= 0.0)
    mHI[ind] = minmass

    sfre = sfrall/mcold
    sfreH2 = np.zeros(shape = len(mcold))
    ind = np.where(mH2 > 0.0)
    sfreH2[ind] = sfrall[ind]/mH2[ind]

    mhrat = mH2/mHI
    mstarden  = mstar / volh
    mstarbden = mstar_burst / volh
    mH2den    = mH2 / volh


    mHIden   = mHI / volh
    print len(mHIden), len(redshifts)
    h = np.zeros(shape = (len(redshifts)))
    omegaHI = np.zeros(shape = (len(redshifts)))
    for z in range(0,len(redshifts)):
        h[z] = h0 * sqrt(OmegaM*pow(1.0+redshifts[z],3.0) + OmegaL)
        omegaHI[z]  = mHIden[z] / (rho_crit*pow(h0,2.0))

    #Assume a gas-dust mass ratio that scales with metallicity. We use the Remy-Ruyer et al. (2013) G/D ratio.
    mdustden = 0.006 * mcold_metals/Zsun / volh

    mcold_plot = np.zeros(shape = len(mcold))
    mhot_plot = np.zeros(shape = len(mcold))
    meje_plot = np.zeros(shape = len(mcold))
    mstar_plot = np.zeros(shape = len(mcold))
    mHI_plot = np.zeros(shape = len(mcold))
    mH2_plot = np.zeros(shape = len(mcold))

    ind = np.where(massbar > 0.0)
    mcold_plot[ind] = np.log10(mcold[ind]/volh/sbar)
    mhot_plot[ind] = np.log10(mhot[ind]/volh/sbar)
    meje_plot[ind] = np.log10(meje[ind]/volh/sbar)
    mstar_plot[ind] = np.log10(mstar[ind]/volh/sbar)
    mHI_plot[ind] = np.log10(mHI[ind]/volh/sbar)
    mH2_plot[ind] = np.log10(mH2[ind]/volh/sbar)

    mcold_dm_plot = np.zeros(shape = len(mcold))
    mhot_dm_plot = np.zeros(shape = len(mcold))
    meje_dm_plot = np.zeros(shape = len(mcold))
    mstar_dm_plot = np.zeros(shape = len(mcold))
    mbar_dm_plot = np.zeros(shape = len(mcold))
    mHI_dm_plot = np.zeros(shape = len(mcold))
    mH2_dm_plot = np.zeros(shape = len(mcold))

    ind = np.where(mDM > 0.0)
    mcold_dm_plot[ind] = np.log10(mcold[ind]/(mDM[ind]+massbar[ind]))
    mhot_dm_plot[ind] = np.log10(mhot[ind]/(mDM[ind]+massbar[ind]))
    meje_dm_plot[ind] = np.log10(meje[ind]/(mDM[ind]+massbar[ind]))
    mstar_dm_plot[ind] = np.log10(mstar[ind]/(mDM[ind]+massbar[ind]))
    mbar_dm_plot[ind] = np.log10(massbar[ind]/(mDM[ind]+massbar[ind]))
    mHI_dm_plot[ind] = np.log10(mHI[ind]/(mDM[ind]+massbar[ind]))
    mH2_dm_plot[ind] = np.log10(mH2[ind]/(mDM[ind]+massbar[ind]))

    return (mstar_plot, mcold_plot, mhot_plot, meje_plot,
     mstar_dm_plot, mcold_dm_plot, mhot_dm_plot, meje_dm_plot, mbar_dm_plot,
     sfr, sfrd, sfrb, mstarden, mstarbden, sfre, sfreH2, mhrat,
     mHI_plot, mH2_plot, mH2den, mdustden, omegaHI)

def plot_mass_densities(plt, outdir, redshifts, mstar, mcold, mhot, meje):

    fig = plt.figure(figsize=(5,5))

    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(m/m_{\\rm bar,total})$"

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, 0, 10, -6, 0, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    ax.plot(redshifts, mstar,'k', label='stars')
    ax.plot(redshifts, mcold,'b', label='cold gas')
    ax.plot(redshifts, mhot,'r', label='hot gas')
    ax.plot(redshifts, meje,'g', label='ejected gas')

    common.prepare_legend(ax, ['k','b','r','g'])
    common.savefig(outdir, fig, "global.pdf")


def plot_baryon_fractions(plt, outdir, redshifts, mstar_dm, mcold_dm, mhot_dm, meje_dm, mbar_dm):

    fig = plt.figure(figsize=(9.5,9.5))

    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(\\rho/\\rho_{\\rm m})$"

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, 0, 10, -4, 0.5, xtit, ytit, locators=(0.1, 1, 0.1))

    ax.plot(redshifts, mstar_dm, 'k', label='stars')
    ax.plot(redshifts, mcold_dm, 'b', label='cold gas')
    ax.plot(redshifts, mhot_dm, 'r', label='hot gas')
    ax.plot(redshifts, meje_dm, 'g', label='ejected gas')
    ax.plot(redshifts, mbar_dm, 'm', label='total baryons')

    yplot = [0.1866920152, 0.1866920152]
    xplot = [0, 10]

    ax.plot(xplot, np.log10(yplot), 'k', linestyle='dashed', label ='Universal $f_{\\rm baryon}$')

    common.prepare_legend(ax, ['k','b','r','g','m','k'])
    common.savefig(outdir, fig, "baryon_frac.pdf")


def plot_cosmic_sfr(plt, outdir, obsdir, redshifts, h0, sfr, sfrd, sfrb):

    fig = plt.figure(figsize=(5,5))

    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(CSFRD/ M_{\odot}\,yr^{-1}\,cMpc^{-3})$"

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, 0, 10, -6, 0, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    redK11, SFRK11, err_upK11, err_dnK11 = common.load_observation(obsdir, 'SFR/SFRD_Karim11.data', [0,1,2,3])

    hobs = 0.7
    xobs = redK11

    yobs = xobs*0. - 999.
    indx = np.where( SFRK11 > 0)
    yobs[indx] = np.log10(SFRK11[indx] * pow(hobs/h0, 2.0))

    lerr = yobs*0. - 999.
    indx = np.where( (SFRK11-err_dnK11) > 0)
    lerr[indx]  = np.log10(SFRK11[indx] - err_dnK11[indx])

    herr = yobs*0. + 999.
    indx = np.where( (SFRK11+err_upK11) > 0)
    herr[indx]  = np.log10(SFRK11[indx] + err_upK11[indx])

    ax.errorbar(xobs[0:8], yobs[0:8], yerr=[yobs[0:8]-lerr[0:8],herr[0:8]-yobs[0:8]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Karim+11 obs")
    ax.errorbar(xobs[9:17], yobs[9:17], yerr=[yobs[9:17]-lerr[9:17],herr[9:17]-yobs[9:17]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='x',label="Karim+11 extr")

    #Driver (Chabrier IMF), ['Baldry+2012, z<0.06']
    redD17d, redD17u, sfrD17, err1, err2, err3, err4 = common.load_observation(obsdir, 'SFR/Driver17_sfr.dat', [0,1,2,3,4,5,6])

    hobs = 0.7
    xobs = (redD17d+redD17u)/2.0
    yobs = sfrD17 + np.log10(hobs/h0)

    err = yobs*0. - 999.
    err = np.sqrt(pow(err1,2.0)+pow(err2,2.0)+pow(err3,2.0)+pow(err4,2.0))
    ax.errorbar(xobs, yobs, yerr=[err,err], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D',label="Driver+17")

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(sfr > 0)
    ax.plot(redshifts[ind], np.log10(sfr[ind]*pow(h0,2.0)), 'k', label ='total')
    ind = np.where(sfrd > 0)
    ax.plot(redshifts[ind], np.log10(sfrd[ind]*pow(h0,2.0)), 'b', linestyle='dashed',label ='quiescent')
    ind = np.where(sfrb > 0)
    ax.plot(redshifts[ind], np.log10(sfrb[ind]*pow(h0,2.0)),'r', linestyle='dotted', label ='bursts')

    common.prepare_legend(ax, ['k','b','r','grey','grey','grey'])
    common.savefig(outdir, fig, "cosmic_sfr.pdf")


def plot_stellar_mass_cosmic_density(plt, outdir, obsdir, redshifts, h0, mstarden, mstarbden):

    # Plots stellar mass cosmic density
    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(\\rho_{\\rm star}/ M_{\odot}\,cMpc^{-3})$"

    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(111)
    common.prepare_ax(ax, 0, 10, 5, 8.7, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(mstarden > 0)
    ax.plot(redshifts[ind],np.log10(mstarden[ind]*pow(h0,2.0)), 'r', label='SHArk all galaxies')
    ind = np.where(mstarbden > 0)
    ax.plot(redshifts[ind],np.log10(mstarbden[ind]*pow(h0,2.0)), 'b', linestyle='dashed', label='mass formed in SBs')

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    redD17d, redD17u, smdD17, err1, err2, err3, err4 = common.load_observation(obsdir, 'SFR/Driver17_smd.dat', [1,2,3,4,5,6,7])

    hobs = 0.7
    xobs = (redD17d+redD17u)/2.0
    yobs = smdD17 + np.log10(hobs/h0)

    err = yobs*0. - 999.
    err = np.sqrt(pow(err1,2.0)+pow(err2,2.0)+pow(err3,2.0)+pow(err4,2.0))
    ax.errorbar(xobs, yobs, yerr=[err,err], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Driver+17")

    common.prepare_legend(ax, ['r','b','grey'])
    common.savefig(outdir, fig, "cosmic_smd.pdf")


def plot_sft_efficiency(plt, outdir, redshifts, sfre, sfreH2, mhrat):

    fig = plt.figure(figsize=(9.5,11))
    xmin, xmax, ymin, ymax = 0, 10, -6, 0

    # panel 1
    ax = fig.add_subplot(311)
    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(SFR/M_{\\rm cold} Gyr^{-1})$"
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(sfre > 0)
    ax.plot(redshifts[ind],np.log10(sfre[ind]), 'r', label='SHArk')
    common.prepare_legend(ax, ['r'])

    #panel 2
    ax = fig.add_subplot(312)
    ytit="$\\rm log_{10}(SFE_{\\rm H_2}/Gyr^{-1})$"
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(sfreH2 > 0)
    ax.plot(redshifts[ind],np.log10(sfreH2[ind]), 'r', label='SHArk')
    common.prepare_legend(ax, ['r'])

    #panel 3
    ax = fig.add_subplot(313)
    ytit="$\\rm log_{10}(M_{\\rm mol}/M_{\\rm atom})$"
    common.prepare_ax(ax, xmin, xmax, -3, 2, xtit, ytit, locators=(0.1, 1, 0.1))

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(mhrat > 0)
    ax.plot(redshifts[ind],np.log10(mhrat[ind]), 'r', label ='SHArk')
    common.prepare_legend(ax, ['r'])

    common.savefig(outdir, fig, "cosmic_sfe.pdf")


def plot_omega_h2(plt, outdir, obsdir, redshifts, h0, mH2den):

    fig = plt.figure(figsize=(5,5))

    ax = fig.add_subplot(111)
    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(\\rho_{\\rm H_2}/ M_{\odot}\,cMpc^{-3})$"
    common.prepare_ax(ax, 0, 6, 5, 9, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(mH2den > 0)
    ax.plot(redshifts[ind], np.log10(mH2den[ind]*pow(h0,2.0)), 'r', label='SHArk')

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    zD16, zloD16, zupD16, rhoH2D16, rhoH2loD16, rhoH2upD16  = common.load_observation(obsdir, 'SFR/Walter17_H2.dat', [0,1,2,3,4,5])

    hobs = 0.7

    xobs    = zD16
    errxlow = zD16-zloD16
    errxup  = zupD16-zD16
    yobs = np.log10(rhoH2D16) + np.log10(pow(hobs/h0,3.0))
    errylow = np.log10(rhoH2D16) - np.log10(rhoH2loD16)
    erryup  = np.log10(rhoH2upD16) - np.log10(rhoH2D16)

    ax.errorbar(xobs, yobs, xerr=[errxlow,errxup], yerr=[errylow,erryup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='+',label="Decarli+16")
    ax.errorbar(xobs[0:1], yobs[0:1], xerr=[errxlow[0:1],errxup[0:1]], yerr=[errylow[0:1],erryup[0:1]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Bosseli+14")

    # Legend
    common.prepare_legend(ax, ['r','grey','grey','grey'])
    common.savefig(outdir, fig, "omega_H2.pdf")


def plot_mass_cosmic_density(plt, outdir, redshifts, mcold, mHI, mH2):

    fig = plt.figure(figsize=(9.5,9.5))

    ax = fig.add_subplot(111)
    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(\\rho_{\\rm neutral}/ \\rho_{\\rm crit,z=0})$"
    common.prepare_ax(ax, 0, 10, -5, -1, xtit, ytit, locators=(0.1, 1, 0.1))

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ax.plot(redshifts, mcold + np.log10(Omegab), 'k', label='total neutral')
    ax.plot(redshifts, mHI + np.log10(Omegab), 'b', label='atomic')
    ax.plot(redshifts, mH2 + np.log10(Omegab), 'r', label='molecular')

    common.prepare_legend(ax, ['k','b','r'])
    common.savefig(outdir, fig, "omega_neutral.pdf")


def plot_cosmic_dust(plt, outdir, obsdir, redshifts, h0, mdustden):

    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(111)
    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(\\rho_{\\rm dust}/ M_{\odot}\,cMpc^{-3})$"
    common.prepare_ax(ax, 0, 10, 4, 6.3, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(mdustden > 0)
    ax.plot(redshifts[ind],np.log10(mdustden[ind]*pow(h0,2.0)),'r', label ='SHArk')

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    redD17d,redD17u,smdD17,err1,err2,err3,err4 = common.load_observation(obsdir, 'SFR/Driver17_dust.dat', [1,2,3,4,5,6,7])

    hobs = 0.7
    xobs = (redD17d+redD17u)/2.0
    yobs = smdD17 + np.log10(hobs/h0)

    err = yobs*0. - 999.
    err = np.sqrt(pow(err1,2.0)+pow(err2,2.0)+pow(err3,2.0)+pow(err4,2.0))

    ax.errorbar(xobs, yobs, yerr=[err,err], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Driver+17")

    common.prepare_legend(ax, ['r','grey'])
    common.savefig(outdir, fig, "cosmic_dust.pdf")


def plot_omega_HI(plt, outdir, obsdir, redshifts, h0, omegaHI):

    fig = plt.figure(figsize=(5,5))

    ax = fig.add_subplot(111)
    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(\\Omega_{\\rm H_I})$"
    common.prepare_ax(ax, 0, 6, -5, -1, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    # note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(omegaHI > 0)
    ax.plot(redshifts[ind], np.log10(omegaHI[ind]*pow(h0,2.0)), 'r', label='SHArk')

    # Rhee+18 compilation
    redR18,reddR18,reduR18,omegaR18,errdnR18,errupR18 = common.load_observation(obsdir, 'Gas/HI_density_for_Claudia.dat', [1,2,3,7,8,9])

    ax.errorbar(redR18,np.log10(omegaR18*1e-3), xerr=[reddR18,reduR18],yerr=[errdnR18,errupR18], ls='None', mfc='None', ecolor = 'grey', mec='grey', marker='o', label="Rhee+18 (comp)")

    common.prepare_legend(ax, ['r','grey','grey','grey'])
    common.savefig(outdir, fig, "omega_HI.pdf")

def main():

    plt = common.load_matplotlib()
    modeldir, outdir, obsdir, snapshot = common.parse_args()

    fields = {'Global': ('redshifts', 'mHI', 'mH2', 'mcold', 'mcold_metals',
                         'mhot_halo', 'mejected_halo', 'mstars', 'mstars_bursts',
                         'mBH', 'SFR_quiescent', 'SFR_burst', 'mDM', 'mcold_halo')}

    hdf5_data = common.read_data(modeldir, snapshot, fields)
    h0, redshifts = hdf5_data[0], hdf5_data[2]

    (mstar_plot, mcold_plot, mhot_plot, meje_plot,
     mstar_dm_plot, mcold_dm_plot, mhot_dm_plot, meje_dm_plot, mbar_dm_plot,
     sfr, sfrd, sfrb, mstarden, mstarbden, sfre, sfreH2, mhrat,
     mHI_plot, mH2_plot, mH2den, mdustden, omegaHI) = prepare_data(hdf5_data, redshifts)

    plot_mass_densities(plt, outdir, redshifts, mstar_plot, mcold_plot, mhot_plot, meje_plot)
    plot_baryon_fractions(plt, outdir, redshifts, mstar_dm_plot, mcold_dm_plot, mhot_dm_plot, meje_dm_plot, mbar_dm_plot)
    plot_cosmic_sfr(plt, outdir, obsdir, redshifts, h0, sfr, sfrd, sfrb)
    plot_stellar_mass_cosmic_density(plt, outdir, obsdir, redshifts, h0, mstarden, mstarbden)
    plot_sft_efficiency(plt, outdir, redshifts, sfre, sfreH2, mhrat)
    plot_mass_cosmic_density(plt, outdir, redshifts, mcold_plot, mHI_plot, mH2_plot)
    plot_omega_h2(plt, outdir, obsdir, redshifts, h0, mH2den)
    plot_cosmic_dust(plt, outdir, obsdir, redshifts, h0, mdustden)
    plot_omega_HI(plt, outdir, obsdir, redshifts, h0, omegaHI)

if __name__ == '__main__':
    main()