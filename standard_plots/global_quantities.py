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
"""Global plots"""

import math

import numpy as np

import common
import utilities_statistics as us


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
XH = 0.72

def prepare_data(hdf5_data, redshifts):

    (h0, volh, _, mHI, mH2, mcold, mcold_metals, mhot, meje, mlost, mcreated, mstar,
     mstar_burst_mergers, mstar_burst_diskins, mBH, sfrdisk, sfrburst, 
     mDM, mcold_halo, number_major_mergers, number_minor_mergers, 
     number_disk_instabil, max_smbh) = hdf5_data


    history_interactions = np.zeros(shape = (3, len(redshifts)))
    history_interactions[0,:] = number_major_mergers[:]/volh
    history_interactions[1,:] = number_minor_mergers[:]/volh
    history_interactions[2,:] = number_disk_instabil[:]/volh

    sfrall = sfrdisk + sfrburst

    maxden = 6.2863*pow(10.0,9.0)
    ind = 0
    for i,j,k,l in zip(mhot, meje, mstar, mcold):
        tot = i+j+k+l
        if(tot/(volh/pow(h0,2.0)) > maxden):
            print("density exceeding maxden", tot/(volh/pow(h0,2.0))/maxden)
        ind = ind +1

    #Add up cold halo component to hot gas.
    mhot = mhot + mcold_halo

    print(mlost)
    massbar = mcold+mhot+meje+mstar+mBH+mlost
    
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
    mcoldden  = mcold / volh
    mhotden   = mhot / volh
    mejeden   = meje / volh
    mDMden    = mDM / volh
    mstarbden_mergers = mstar_burst_mergers / volh
    mstarbden_diskins = mstar_burst_diskins / volh

    mH2den    = mH2 / volh 
    mHIden    = mHI / volh / XH

    h = np.zeros(shape = (len(redshifts)))
    omegaHI = np.zeros(shape = (len(redshifts)))
    for z in range(0,len(redshifts)):
        h[z] = h0 * math.sqrt(OmegaM*pow(1.0+redshifts[z],3.0) + OmegaL)
        omegaHI[z]  = mHIden[z] / (rho_crit*pow(h0,2.0))

    #Assume a gas-dust mass ratio that scales with metallicity. We use the Remy-Ruyer et al. (2013) G/D ratio.
    mdustden     = 0.006 * mcold_metals/Zsun / volh
    mdustden_mol = 0.006 * mcold_metals/Zsun * mH2 / (mH2 + mHI) / volh

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
    mlost_dm_plot = np.zeros(shape = len(mcold))
    mcreated_dm_plot = np.zeros(shape = len(mcold))

    ind = np.where(mDM > 0.0)
    mcold_dm_plot[ind] = np.log10(mcold[ind]/(mDM[ind]+massbar[ind]))
    mhot_dm_plot[ind] = np.log10(mhot[ind]/(mDM[ind]+massbar[ind]))
    meje_dm_plot[ind] = np.log10(meje[ind]/(mDM[ind]+massbar[ind]))
    mstar_dm_plot[ind] = np.log10(mstar[ind]/(mDM[ind]+massbar[ind]))
    mbar_dm_plot[ind] = np.log10(massbar[ind]/(mDM[ind]+massbar[ind]))
    mHI_dm_plot[ind] = np.log10(mHI[ind]/(mDM[ind]+massbar[ind]))
    mH2_dm_plot[ind] = np.log10(mH2[ind]/(mDM[ind]+massbar[ind]))
    mlost_dm_plot[ind] = np.log10(mlost[ind]/(mDM[ind]+massbar[ind]))
    mcreated_dm_plot[ind] = np.log10(mcreated[ind]/(mDM[ind]+massbar[ind]))
  
    return (mstar_plot, mcold_plot, mhot_plot, meje_plot,
     mstar_dm_plot, mcold_dm_plot, mhot_dm_plot, meje_dm_plot, mbar_dm_plot,
     sfr, sfrd, sfrb, mstarden, mstarbden_mergers, mstarbden_diskins, sfre, sfreH2, mhrat,
     mHI_plot, mH2_plot, mH2den, mdustden, omegaHI, mdustden_mol, mcoldden, mhotden, mejeden,
     history_interactions, mDMden, mlost_dm_plot, mcreated_dm_plot)

def plot_mass_densities(plt, outdir, obsdir, h0, redshifts, mstar, mcold, mhot, meje, mstarden, mcoldden, mhotden, mejeden):

    fig = plt.figure(figsize=(5,15))

    xtit="$\\rm Lookback\, time/Gyr$"
    ytit="$\\rm log_{10}(\\rho/\\rho_{\\rm bar,halos})$"

    ax = fig.add_subplot(411)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, 0, 13.5, -5, 0.1, xtit, ytit, locators=(0.1, 1, 0.1, 1), fontsize=10)
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    new_tick_locations = np.array([0., 2., 4., 6., 8., 10., 12.])

    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(us.redshift(new_tick_locations), fontsize=12)

    ax2.set_xlabel("redshift",fontsize=13)

    ax.plot(us.look_back_time(redshifts), mstar,'k', label='stars')
    ax.plot(us.look_back_time(redshifts), mcold,'b', label='ISM gas')
    ax.plot(us.look_back_time(redshifts), mhot,'r', label='halo gas')
    ax.plot(us.look_back_time(redshifts), meje,'g', label='ejected gas')

    common.prepare_legend(ax, ['k','b','r','g'], fontsize=10)

    xtit="$\\rm Lookback\, time/Gyr$"
    ytit="$\\rm log_{10}(\\rho_{\\rm bar} /\\rm M_{\odot} \\rm Mpc^{-3})$"

    ax = fig.add_subplot(412)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, 0, 13.5, 6, 10, xtit, ytit, locators=(0.1, 1, 0.1, 1), fontsize=10)
    
    msharktot = mstarden + mcoldden + mhotden + mejeden
    ind = np.where(mstarden > 0)
    ax.plot(us.look_back_time(redshifts), np.log10(mstarden[ind]*pow(h0,2.0)),'k', label='Shark')
    ind = np.where(mcoldden > 0)
    ax.plot(us.look_back_time(redshifts), np.log10(mcoldden[ind]*pow(h0,2.0)),'b')
    ind = np.where(mhotden > 0)
    ax.plot(us.look_back_time(redshifts), np.log10(mhotden[ind]*pow(h0,2.0)),'r')
    ind = np.where(mejeden > 0)
    ax.plot(us.look_back_time(redshifts), np.log10(mejeden[ind]*pow(h0,2.0)),'g')
    ind = np.where(msharktot > 0)
    ax.plot(us.look_back_time(redshifts), np.log10(msharktot[ind]*pow(h0,2.0)),'DarkMagenta')

    lbt, eaglesm, eaglesmout, eagleism, eaglehg, eagleejec = common.load_observation(obsdir, 'Models/OtherModels/EAGLE_BaryonGrowthTotal.dat', [0,2,3,4,5,6])
    eaglesm  = np.log10(pow(10.0, eaglesm) + pow(10.0, eaglesmout))
    eagletot = np.log10(pow(10.0, eaglesm) + pow(10.0, eagleism) +  pow(10.0, eaglehg) + pow(10.0, eagleejec))
    ind = np.where(eaglesm > 0)
    ax.plot(lbt[ind], eaglesm[ind] - 6.0, 'k', linestyle='dotted', label ='EAGLE')
    ind = np.where(eagleism > 0)
    ax.plot(lbt[ind], eagleism[ind]- 6.0, 'b', linestyle='dotted')
    ind = np.where(eaglehg > 0)
    ax.plot(lbt[ind], eaglehg[ind] - 6.0, 'r', linestyle='dotted')
    ind = np.where(eagleejec > 0)
    ax.plot(lbt[ind], eagleejec[ind] - 6.0, 'g', linestyle='dotted')
    ind = np.where(eagletot > 0)
    ax.plot(lbt[ind], eagletot[ind] - 6.0, 'DarkMagenta', linestyle='dotted')

    common.prepare_legend(ax, ['k','k','k'], fontsize=10)

    ax = fig.add_subplot(413)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, 0, 13.5, 6, 10, xtit, ytit, locators=(0.1, 1, 0.1, 1), fontsize=10)

    ind = np.where(mstarden > 0)
    ax.plot(us.look_back_time(redshifts), np.log10(mstarden[ind]*pow(h0,2.0)),'k', label='Shark')
    ind = np.where(mcoldden > 0)
    ax.plot(us.look_back_time(redshifts), np.log10(mcoldden[ind]*pow(h0,2.0)),'b')
    ind = np.where(mhotden > 0)
    ax.plot(us.look_back_time(redshifts), np.log10(mhotden[ind]*pow(h0,2.0)),'r')
    ind = np.where(mejeden > 0)
    ax.plot(us.look_back_time(redshifts), np.log10(mejeden[ind]*pow(h0,2.0)),'g')
    ind = np.where(msharktot > 0)
    ax.plot(us.look_back_time(redshifts), np.log10(msharktot[ind]*pow(h0,2.0)),'DarkMagenta')

    lbt, galsm, galism, galhg, galejec = common.load_observation(obsdir, 'Models/OtherModels/global_Mitchell18.dat', [0,2,3,4,5])
    galtot  = np.log10(pow(10.0, galsm) + pow(10.0, galism) + pow(10.0, galhg) + pow(10.0, galejec))
    galsm   = (galsm)
    galism  = (galism)
    galhg   = (galhg)
    galejec = (galejec)

    h = 0.704
    vol = 6.0 #np.log10(1953125.0/pow(h,3.0))
    ind = np.where(galsm > 0)
    ax.plot(lbt[ind], galsm[ind] - vol, 'k', linestyle='dashed', label ='GALFORM M18', linewidth=1)
    ind = np.where(galism > 0)
    ax.plot(lbt[ind], galism[ind]- vol, 'b', linestyle='dashed', linewidth=1)
    ind = np.where(galhg > 0)
    ax.plot(lbt[ind], galhg[ind] - vol, 'r', linestyle='dashed', linewidth=1)
    ind = np.where(galejec > 0)
    ax.plot(lbt[ind], galejec[ind] - vol, 'g', linestyle='dashed', linewidth=1)
    ind = np.where(galtot > 0)
    ax.plot(lbt[ind], galtot[ind] - vol, 'DarkMagenta', linestyle='dashed', linewidth=1)


    common.prepare_legend(ax, ['k','k','k'], fontsize=10)

    ax = fig.add_subplot(414)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, 0, 13.5, 6, 10, xtit, ytit, locators=(0.1, 1, 0.1, 1), fontsize=10)

    ind = np.where(mstarden > 0)
    ax.plot(us.look_back_time(redshifts), np.log10(mstarden[ind]*pow(h0,2.0)),'k', label='Shark')
    ind = np.where(mcoldden > 0)
    ax.plot(us.look_back_time(redshifts), np.log10(mcoldden[ind]*pow(h0,2.0)),'b')
    ind = np.where(mhotden > 0)
    ax.plot(us.look_back_time(redshifts), np.log10(mhotden[ind]*pow(h0,2.0)),'r')
    ind = np.where(mejeden > 0)
    ax.plot(us.look_back_time(redshifts), np.log10(mejeden[ind]*pow(h0,2.0)),'g')
    ind = np.where(msharktot > 0)
    ax.plot(us.look_back_time(redshifts), np.log10(msharktot[ind]*pow(h0,2.0)),'DarkMagenta')

    zbt, galism, galsm, galhg, galejec = common.load_observation(obsdir, 'Models/OtherModels/BaryonBudgetLgalaxies.dat', [1,2,3,4,5])
    lbt = us.look_back_time(zbt)
    galtot  = np.log10(galsm + galism + galhg + galejec) + 10.0
    galsm   = np.log10(galsm)+10.0
    galism  = np.log10(galism)+10.0
    galhg   = np.log10(galhg)+10.0
    galejec = np.log10(galejec)+10.0

    h = 0.673
    vol = np.log10(125000000.0/pow(h,3.0))
    ind = np.where(galsm > 0)
    ax.plot(lbt[ind], galsm[ind] - vol - np.log10(h), 'k', linestyle='dashdot', label ='L-galaxies H15', linewidth=1)
    ind = np.where(galism > 0)
    ax.plot(lbt[ind], galism[ind]- vol - np.log10(h), 'b', linestyle='dashdot', linewidth=1)
    ind = np.where(galhg > 0)
    ax.plot(lbt[ind], galhg[ind] - vol - np.log10(h), 'r', linestyle='dashdot', linewidth=1)
    ind = np.where(galejec > 0)
    ax.plot(lbt[ind], galejec[ind] - vol - np.log10(h), 'g', linestyle='dashdot', linewidth=1)
    ind = np.where(galtot > 0)
    ax.plot(lbt[ind], galtot[ind] - vol - np.log10(h), 'DarkMagenta', linestyle='dashdot', linewidth=1)

    common.prepare_legend(ax, ['k','k','k'], fontsize=10)


    common.savefig(outdir, fig, "global.pdf")

    fig = plt.figure(figsize=(5,4))

    xtit="$\\rm Lookback\, time/Gyr$"
    ytit="$\\rm log_{10}(\\rho_{\\rm gas, halo} /\\rm M_{\odot} \\rm Mpc^{-3})$"

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)
    common.prepare_ax(ax, 0, 13.5, 5, 10, xtit, ytit, locators=(0.1, 1, 0.1, 1), fontsize=10)

    ind = np.where(mhotden > 0)
    ax.plot(us.look_back_time(redshifts), np.log10(mhotden[ind]*pow(h0,2.0)),'r', label='Shark')

    lbt, eaglesm, eaglesmout, eagleism, eaglehg, eagleejec = common.load_observation(obsdir, 'Models/OtherModels/EAGLE_BaryonGrowthTotal.dat', [0,2,3,4,5,6])
    eaglesm  = np.log10(pow(10.0, eaglesm) + pow(10.0, eaglesmout))
    eagletot = np.log10(pow(10.0, eaglesm) + pow(10.0, eagleism) +  pow(10.0, eaglehg) + pow(10.0, eagleejec))
    ind = np.where(eaglehg > 0)
    ax.plot(lbt[ind], eaglehg[ind] - 6.0, 'r', linestyle='dotted', label ='EAGLE')


    lbt, galsm, galism, galhg, galejec = common.load_observation(obsdir, 'Models/OtherModels/global_Mitchell18.dat', [0,2,3,4,5])
    galtot  = np.log10(pow(10.0, galsm) + pow(10.0, galism) + pow(10.0, galhg) + pow(10.0, galejec))
    galhg   = (galhg)

    h = 0.704
    vol = 6.0 #np.log10(1953125.0/pow(h,3.0))
    ind = np.where(galhg > 0)
    ax.plot(lbt[ind], galhg[ind] - vol, 'r', linestyle='dashed', linewidth=1, label ='GALFORM M18')
    ind = np.where(galejec > 0)

    zbt, galism, galsm, galhg, galejec = common.load_observation(obsdir, 'Models/OtherModels/BaryonBudgetLgalaxies.dat', [1,2,3,4,5])
    lbt = us.look_back_time(zbt)
    galhg   = np.log10(galhg)+10.0

    h = 0.673
    vol = np.log10(125000000.0/pow(h,3.0))
    ind = np.where(galhg > 0)
    ax.plot(lbt[ind], galhg[ind] - vol - np.log10(h), 'r', linestyle='dashdot', linewidth=1, label ='L-galaxies H15')

    common.prepare_legend(ax, ['k','k','k','k'], fontsize=10)


    common.savefig(outdir, fig, "global_hotgas.pdf")


def plot_baryon_fractions(plt, outdir, redshifts, mstar_dm, mcold_dm, mhot_dm, meje_dm, mbar_dm, mlost_dm, mcreated_dm):

    fig = plt.figure(figsize=(9.5,9.5))

    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(\\rho/\\rho_{\\rm m})$"

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, 0, 10, -4, 0.5, xtit, ytit, locators=(0.1, 1, 0.1))

    ax.plot(redshifts, mstar_dm, 'k', label='stars')
    ax.plot(redshifts, mcold_dm, 'b', label='ISM gas')
    ax.plot(redshifts, mhot_dm, 'r', label='halo gas')
    ax.plot(redshifts, meje_dm, 'g', label='ejected gas')
    ax.plot(redshifts, mlost_dm, 'orange', label='lost')
    ax.plot(redshifts, mcreated_dm, 'orange', linestyle='dashed', label='created')

    ax.plot(redshifts, mbar_dm, 'm', label='total baryons')

    yplot = [0.1866920152, 0.1866920152]
    xplot = [0, 10]

    ax.plot(xplot, np.log10(yplot), 'k', linestyle='dashed', label ='Universal $f_{\\rm baryon}$')

    common.prepare_legend(ax, ['k','b','r','g','orange','orange','m','k'])
    common.savefig(outdir, fig, "baryon_frac.pdf")


def plot_cosmic_sfr(plt, outdir, obsdir, redshifts, h0, sfr, sfrd, sfrb, history_interactions, mDMden):

    fig = plt.figure(figsize=(5,9))

    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(CSFRD/ M_{\odot}\,yr^{-1}\,cMpc^{-3})$"

    ax = fig.add_subplot(211)
    plt.subplots_adjust(left=0.15)

    common.prepare_ax(ax, 0, 10, -3, -0.5, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #Madau & Dickinson 2014
    reddnM14, redupM14, sfrM14, sfrM14errup, sfrM14errdn = common.load_observation(obsdir, 'Global/SFRD_Madau14.dat', [0,1,2,3,4])
    #authors assume a Salpeter IMF, so a correction of np.log10(0.63) is necessary.
    sfrM14errdn = abs(sfrM14errdn)
    hobs = 0.7
    sfrM14 = sfrM14 + np.log10(pow(hobs/h0, 2.0)) + np.log10(0.63)
    ax.errorbar((reddnM14 + redupM14) / 2.0, sfrM14, xerr=abs(redupM14-reddnM14)/2.0, yerr=[sfrM14errdn, sfrM14errup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D', markersize=1.5)

    #D'Silva+23 (Chabrier IMF)
    sfrD23, sfrD23errup, sfrD23errdn, zD23, zD23errup, zD23errdn = common.load_observation(obsdir, 'Global/DSilva23_sfr.dat', [0,1,2,3,4,5])


    #Adams23 (Chabrier IMF)
    zA23, sfrA23, sfrA23errdn, sfrA23errup = common.load_observation(obsdir, 'Global/Adams23_CSFRDCompilation.dat', [0,1,2,3])
    sfrA23errdn = sfrA23 - sfrA23errdn #make them relative errors
    sfrA23errup = sfrA23errup - sfrA23
    hobs = 0.7
    yobsA23 = sfrA23 + np.log10(hobs/h0)
    ax.errorbar(zA23, yobsA23, yerr=[sfrA23errdn, sfrA23errup], ls='None', mfc='None', ecolor = 'darkgreen', mec='darkgreen',marker='s')



    hobs = 0.7
    yobsD23 = sfrD23 + np.log10(hobs/h0)
    ax.errorbar(zD23, yobsD23, xerr=[zD23errdn, zD23errup], yerr=[sfrD23errdn, sfrD23errup], ls='None', mfc='None', ecolor = 'navy', mec='navy',marker='o')

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(sfr > 0)
    ax.plot(redshifts[ind], np.log10(sfr[ind]*pow(h0,2.0)), 'k', linewidth=1, label ='total')
    #print("SFR density of the Universe")
    #for a,b in zip(redshifts[ind],  np.log10(sfr[ind]*pow(h0,2.0))):
    #    print(a,b)

    ind = np.where(sfrd > 0)
    ax.plot(redshifts[ind], np.log10(sfrd[ind]*pow(h0,2.0)), 'b', linestyle='dashed', linewidth=1, label ='quiescent')
    ind = np.where(sfrb > 0)
    ax.plot(redshifts[ind], np.log10(sfrb[ind]*pow(h0,2.0)),'r', linestyle='dotted',  linewidth=1, label ='bursts')

    z, sfr_modelvar = common.load_observation(obsdir, 'Models/SharkVariations/Global_OtherModels.dat', [0, 3])
    sfr_modelvar_burst3 = sfr_modelvar[0:179]
    sfr_modelvar_nu0p5  = sfr_modelvar[179:359]
    sfr_modelvar_burst20= sfr_modelvar[360:539]

    ind = np.where(sfr_modelvar_burst20 > -10)
    ax.plot(z[ind], sfr_modelvar_burst20[ind], 'Sienna', linestyle='dotted', label ='$\\eta_{\\rm burst}=20$')
    ind = np.where(sfr_modelvar_burst3 > -10)
    ax.plot(z[ind], sfr_modelvar_burst3[ind], 'DarkSlateGray', linestyle='dashdot', label ='$\\eta_{\\rm burst}=3$')
    ind = np.where(sfr_modelvar_nu0p5 > -10)
    ax.plot(z[ind], sfr_modelvar_nu0p5[ind], 'SlateGray', linestyle='dotted', label ='$\\nu_{\\rm SF}=0.5 \\rm Gyr^{-1}$')

    common.prepare_legend(ax, ['k','b','r','Sienna','DarkSlateGray','SlateGray','grey','grey','grey'], bbox_to_anchor=(0.52, 0.47))

    xtit="$\\rm Lookback\, time/Gyr$"
    ax = fig.add_subplot(212)
    plt.subplots_adjust(left=0.15)

    common.prepare_ax(ax, 0, 13.5, -3, -0.5, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    lbtdown = us.look_back_time(reddnM14)
    lbtup = us.look_back_time(redupM14)
    ax.errorbar(us.look_back_time((reddnM14 + redupM14) / 2.0), sfrM14, xerr=abs(lbtdown-lbtup)/2.0, yerr=[sfrM14errdn, sfrM14errup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D', markersize=1.5, label='Madau+14')
    ax.errorbar(us.look_back_time(zD23), yobsD23, yerr=[sfrD23errdn, sfrD23errup], ls='None', mfc='None', ecolor = 'navy', mec='navy',marker='o', label="D'Silva+23")

    ind = np.where(sfr > 0)
    ax.plot(us.look_back_time(redshifts[ind]), np.log10(sfr[ind]*pow(h0,2.0)), 'k',  linewidth=1)

    
    #print("cosmic SFR")
    #for a,b,c,d,e in zip(redshifts[ind], us.look_back_time(redshifts[ind]), np.log10(sfr[ind]*pow(h0,2.0)), np.log10(sfrd[ind]*pow(h0,2.0)), np.log10(sfrb[ind]*pow(h0,2.0))):
    #    print(a,b,c,d,e)
    #print("finished printing cosmic sfr")
    ind = np.where(sfrd > 0)
    ax.plot(us.look_back_time(redshifts[ind]), np.log10(sfrd[ind]*pow(h0,2.0)), 'b', linestyle='dashed', linewidth=1)
    ind = np.where(sfrb > 0)
    ax.plot(us.look_back_time(redshifts[ind]), np.log10(sfrb[ind]*pow(h0,2.0)),'r', linestyle='dotted',  linewidth=1)

    ind = np.where(sfr_modelvar_burst20 > -10)
    ax.plot(us.look_back_time(z[ind]), sfr_modelvar_burst20[ind], 'Sienna', linestyle='dotted')
    ind = np.where(sfr_modelvar_burst3 > -10)
    ax.plot(us.look_back_time(z[ind]), sfr_modelvar_burst3[ind], 'DarkSlateGray', linestyle='dashdot')
    ind = np.where(sfr_modelvar_nu0p5 > -10)
    ax.plot(us.look_back_time(z[ind]), sfr_modelvar_nu0p5[ind], 'SlateGray', linestyle='dotted')

    common.prepare_legend(ax, ['grey','navy','grey'], loc=2)
    common.savefig(outdir, fig, "cosmic_sfr.pdf")


    #create plot with interaction history
    fig = plt.figure(figsize=(5,4))

    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(density\,rate/Mpc^{-3} h^{-3} Gyr^{-1})$"

    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.15)

    common.prepare_ax(ax, 0, 10, -4, 0, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    delta_time = np.zeros(shape = (len(redshifts))) 
    for i in range (0,len(redshifts)):
        if(i == 0):
               delta_time[i] = 13.7969 - us.look_back_time(redshifts[i])
        if(i < len(redshifts)):
               delta_time[i] = us.look_back_time(redshifts[i-1])-us.look_back_time(redshifts[i])

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(history_interactions[0,:]+history_interactions[1,:] > 0)
    yplot = np.log10((history_interactions[0,ind]+history_interactions[1,ind])/delta_time[ind])
    ax.plot(redshifts[ind], yplot[0], 'k', linewidth=1, label ='mergers')
    ind = np.where(history_interactions[2,:] > 0)
    yplot = np.log10(history_interactions[2,ind]/delta_time[ind])
    ax.plot(redshifts[ind], yplot[0], 'b', linewidth=1, label ='disk instabilities')

    ind = np.where(mDMden[:] > 0)
    yplot = np.log10(mDMden[ind])-11.0
    ax.plot(redshifts[ind], yplot, 'r', linewidth=1, label ='DM mass(-11dex)')

    common.prepare_legend(ax, ['k','b','r'], loc=3)

    common.savefig(outdir, fig, "interaction_history.pdf")

    #create plot with interaction history
    fig = plt.figure(figsize=(5,4))

    xtit="$\\rm redshift$"
    ytit="$\\rm delta\, time/Gyr$"

    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.15)

    common.prepare_ax(ax, 0, 10, 0, 0.3, xtit, ytit, locators=(0.1, 1, 0.1, 0.1))

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ax.plot(redshifts, delta_time, 'k', linewidth=1)

    common.savefig(outdir, fig, "delta_time_history.pdf")

    fig = plt.figure(figsize=(6,5.5))

    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(CSFRD/ M_{\odot}\,yr^{-1}\,cMpc^{-3})$"

    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.15)

    common.prepare_ax(ax, 0, 15, -6, -0.3, xtit, ytit, locators=(1, 1, 1, 1))
    #plt.xscale('log')

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    reddnM14, redupM14, sfrM14, sfrM14errup, sfrM14errdn = common.load_observation(obsdir, 'Global/SFRD_Madau14.dat', [0,1,2,3,4])
    #authors assume a Salpeter IMF, so a correction of np.log10(0.63) is necessary.
    sfrM14errdn = abs(sfrM14errdn)
    hobs = 0.7
    sfrM14 = sfrM14 + np.log10(pow(hobs/h0, 2.0)) + np.log10(0.63)
    #ax.errorbar((reddnM14 + redupM14) / 2.0, sfrM14, xerr=abs(redupM14-reddnM14)/2.0, yerr=[sfrM14errdn, sfrM14errup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D', markersize=1.5, label='Madau+14')

    ax.errorbar(zD23, yobsD23, xerr=[zD23errdn, zD23errup], yerr=[sfrD23errdn, sfrD23errup], ls='None', mfc='None', ecolor = 'navy', mec='navy',marker='o', label ='D\'Silva+23')
    ax.errorbar(zA23, yobsA23, yerr=[sfrA23errdn, sfrA23errup], ls='None', mfc='None', ecolor = 'darkgreen', mec='darkgreen',marker='s', label = 'Adams+23')
    #Driver (Chabrier IMF)
    redD17d, redD17u, sfrD17, err1, err2, err3, err4 = common.load_observation(obsdir, 'Global/Driver18_sfr.dat', [0,1,2,3,4,5,6])
    hobs = 0.7
    xobsD17 = (redD17d+redD17u)/2.0
    yobsD17 = sfrD17 + np.log10(hobs/h0)
    errD17 = yobsD17*0. - 999.
    errD17 = np.sqrt(pow(err1,2.0)+pow(err2,2.0)+pow(err3,2.0)+pow(err4,2.0))
    ax.errorbar(xobsD17, yobsD17, yerr=[errD17,errD17], ls='None', mfc='None', ecolor = 'darkorange', mec='darkorange',marker='o', label = 'Driver+18')


    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(sfr > 0)
    ax.plot(redshifts[ind], np.log10(sfr[ind]*pow(h0,2.0)), 'k', linewidth=1, label ='total (v2.0)')

    ind = np.where(sfrd > 0)
    ax.plot(redshifts[ind], np.log10(sfrd[ind]*pow(h0,2.0)), 'b', linestyle='solid', linewidth=1, label ='disks')
    ind = np.where(sfrb > 0)
    ax.plot(redshifts[ind], np.log10(sfrb[ind]*pow(h0,2.0)),'r', linestyle='solid',  linewidth=1, label ='bursts')
    #print("#will print SFR density")
    #print("#redshift LBG/Gyr SFRD_tot[Msun/yr/Mpc^3] SFRD_disk[Msun/yr/Mpc^3] SFRD_bulges[Msun/yr/Mpc^3]")
    #for a,b,c,d,e in zip(redshifts, us.look_back_time(redshifts), np.log10(sfr*pow(h0,2.0)), np.log10(sfrd*pow(h0,2.0)), np.log10(sfrb*pow(h0,2.0))):
    #    print(a,b,c,d,e)


    zin, sfr_l18, sfr_l18d, sfr_l18b = common.load_observation(obsdir, 'Models/SharkVariations/Global_Lagos18.dat', [0, 2, 3, 4])
    #zin, sfr_l18, sfr_l18d, sfr_l18b = common.load_observation(obsdir, 'Models/SharkVariations/Global_SFR_Lagos23_OldTrees.dat', [0, 2, 3, 4])

    ax.plot(zin, sfr_l18, 'k', linewidth=1, linestyle='dashed', label ='v1.1 (L18)')
    ax.plot(zin, sfr_l18d, 'b', linewidth=1,linestyle='dashed') #,  label ='disks (L18)')
    ax.plot(zin, sfr_l18b, 'r', linewidth=1,linestyle='dashed') #,  label ='bursts (L18)')

    common.prepare_legend(ax, ['k','b','r','k','navy', 'darkgreen','DarkOrange'], loc=3) #bbox_to_anchor=(0.52, 0.47))

    common.savefig(outdir, fig, "cosmic_sfr_compL18.pdf")


def plot_stellar_mass_cosmic_density(plt, outdir, obsdir, redshifts, h0, mstarden, mstarbden_mergers, mstarbden_diskins):

    # Plots stellar mass cosmic density
    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(\\rho_{\\rm star}/ M_{\odot}\,cMpc^{-3})$"

    fig = plt.figure(figsize=(5,9))
    ax = fig.add_subplot(211)
    plt.subplots_adjust(left=0.15)

    common.prepare_ax(ax, 0, 10, 5, 8.7, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(mstarden > 0)
    ax.plot(redshifts[ind],np.log10(mstarden[ind]*pow(h0,2.0)), 'k')

    mstardisk = mstarden - (mstarbden_mergers+mstarbden_diskins)

    ind = np.where(mstarbden_mergers + mstarbden_diskins > 0)
    ax.plot(redshifts[ind],np.log10((mstarbden_mergers[ind]+mstarbden_diskins[ind])*pow(h0,2.0)), 'r', linestyle='dashed')
    ind = np.where(mstardisk > 0)
    ax.plot(redshifts[ind],np.log10(mstardisk[ind]*pow(h0,2.0)), 'b', linestyle='dotted')

    print("#will print stellar density")
    print("#redshift LBG/Gyr SMD_tot[Msun/Mpc^3] SMD_disk[Msun/Mpc^3] SMD_bulge_mergers[Msun/Mpc^3] SMD_bulge_diskins[Msun/Mpc^3]")
    for a,b,c,d,e,f in zip(redshifts, us.look_back_time(redshifts), np.log10(mstarden*pow(h0,2.0)), np.log10(mstardisk*pow(h0,2.0)), np.log10(mstarbden_mergers*pow(h0,2.0)), np.log10(mstarbden_diskins*pow(h0,2.0))):
        print(a,b,c,d,e,f)

    z, sm_modelvar = common.load_observation(obsdir, 'Models/SharkVariations/Global_OtherModels.dat', [0, 4])
    sm_modelvar_burst3  = sm_modelvar[0:179]
    sm_modelvar_nu0p5   = sm_modelvar[181:360]
    sm_modelvar_burst20 = sm_modelvar[360:539]

    ind = np.where(sm_modelvar_burst20 > -10)
    ax.plot(z[ind], sm_modelvar_burst20[ind], 'Sienna', linestyle='dotted')
    ind = np.where(sm_modelvar_burst3 > -10)
    ax.plot(z[ind], sm_modelvar_burst3[ind], 'DarkSlateGray', linestyle='dashdot')
    ind = np.where(sm_modelvar_nu0p5 > -10)
    ax.plot(z[ind], sm_modelvar_nu0p5[ind], 'SlateGray', linestyle='dotted')

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    redD17d, redD17u, smdD17, err1, err2, err3, err4 = common.load_observation(obsdir, 'Global/Driver18_smd.dat', [1,2,3,4,5,6,7])

    hobs = 0.7
    xobs = (redD17d+redD17u)/2.0
    yobs = smdD17 + np.log10(hobs/h0)

    err = yobs*0. - 999.
    err = np.sqrt(pow(err1,2.0)+pow(err2,2.0)+pow(err3,2.0)+pow(err4,2.0))
    ax.errorbar(xobs, yobs, yerr=[err,err], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o', label="Driver+18")

    common.prepare_legend(ax, ['grey'], loc=3)


    xtit="$\\rm Lookback\, time/Gyr$"
    ax = fig.add_subplot(212)
    plt.subplots_adjust(left=0.15)
    common.prepare_ax(ax, 0, 13.5, 5, 8.7, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(mstarden > 0)
    ax.plot(us.look_back_time(redshifts[ind]),np.log10(mstarden[ind]*pow(h0,2.0)), 'k', label='Shark all galaxies')
    mstardisk = mstarden - (mstarbden_mergers+mstarbden_diskins)

    ind = np.where(mstarbden_mergers + mstarbden_diskins> 0)
    ax.plot(us.look_back_time(redshifts[ind]),np.log10((mstarbden_mergers[ind] + mstarbden_diskins[ind])*pow(h0,2.0)), 'r', linestyle='dashed', label='formed in bursts')
    ind = np.where(mstardisk > 0)
    ax.plot(us.look_back_time(redshifts[ind]),np.log10(mstardisk[ind]*pow(h0,2.0)), 'b', linestyle='dotted', label='formed in disks')

    ind = np.where(sm_modelvar_burst20 > -10)
    ax.plot(us.look_back_time(z[ind]), sm_modelvar_burst20[ind], 'Sienna', linestyle='dotted',  label ='$\\eta_{\\rm burst}=20$')
    ind = np.where(sm_modelvar_burst3 > -10)
    ax.plot(us.look_back_time(z[ind]), sm_modelvar_burst3[ind], 'DarkSlateGray', linestyle='dashdot', label ='$\\eta_{\\rm burst}=3$')
    ind = np.where(sm_modelvar_nu0p5 > -10)
    ax.plot(us.look_back_time(z[ind]), sm_modelvar_nu0p5[ind], 'SlateGray', linestyle='dotted', label ='$\\nu_{\\rm SF}=0.5 \\rm Gyr^{-1}$')

    ax.errorbar(us.look_back_time(xobs), yobs, yerr=[err,err], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

    common.prepare_legend(ax, ['k','r','b', 'Sienna','DarkSlateGray','SlateGray'], loc=3)

    common.savefig(outdir, fig, "cosmic_smd.pdf")

    # Plots stellar mass cosmic density
    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(\\rho_{\\rm star}/ M_{\odot}\,cMpc^{-3})$"

    fig = plt.figure(figsize=(6,5.5))
    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.15)

    common.prepare_ax(ax, 0, 15, 1, 9, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(mstarden > 0)
    ax.plot(redshifts[ind],np.log10(mstarden[ind]*pow(h0,2.0)), 'k', label='total (v2.0)')
    mstardisk = mstarden - (mstarbden_mergers+mstarbden_diskins)

    #for a,b,c,d in zip(us.look_back_time(redshifts), mstardisk/mstarden, mstarbden_mergers/mstarden, mstarbden_diskins/mstarden):
    #    print(a,b,c,d)


    ind = np.where(mstardisk > 0)
    ax.plot(redshifts[ind],np.log10(mstardisk[ind]*pow(h0,2.0)), 'b', linestyle='solid', label='disks')

    ind = np.where(mstarbden_mergers + mstarbden_diskins > 0)
    ax.plot(redshifts[ind],np.log10((mstarbden_mergers[ind]+mstarbden_diskins[ind])*pow(h0,2.0)), 'r', linestyle='solid', label='bulges')


    zin, sd_l18, sd_l18d, sd_l18b = common.load_observation(obsdir, 'Models/SharkVariations/Global_SMD_Lagos18.dat', [0, 2, 3, 4])
    ax.plot(zin, sd_l18, 'k', linewidth=1, linestyle='dashed', label ='v1.1 (L18)')
    ax.plot(zin, sd_l18d, 'b', linewidth=1,linestyle='dashed') #,  label ='disks (L18)')
    ax.plot(zin, sd_l18b, 'r', linewidth=1,linestyle='dashed') #,  label ='bursts (L18)')


    reddnM14, redupM14, sfrM14, sfrM14errup, sfrM14errdn = common.load_observation(obsdir, 'Global/SMD_Madau14.dat', [0,1,2,3,4])
    #authors assume a Salpeter IMF, so a correction of np.log10(0.63) is necessary.
    sfrM14errdn = abs(sfrM14errdn)
    hobs = 0.7
    sfrM14 = sfrM14 + np.log10(pow(hobs/h0, 2.0)) + np.log10(0.63)
    ax.errorbar((reddnM14 + redupM14) / 2.0, sfrM14, xerr=abs(redupM14-reddnM14)/2.0, yerr=[sfrM14errdn, sfrM14errup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D', markersize=1.5, label='Madau+14')

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    redD17d, redD17u, smdD17, err1, err2, err3, err4 = common.load_observation(obsdir, 'Global/Driver18_smd.dat', [1,2,3,4,5,6,7])

    hobs = 0.7
    xobs = (redD17d+redD17u)/2.0
    yobs = smdD17 + np.log10(hobs/h0)

    err = yobs*0. - 999.
    err = np.sqrt(pow(err1,2.0)+pow(err2,2.0)+pow(err3,2.0)+pow(err4,2.0))
    ax.errorbar(xobs, yobs, yerr=[err,err], ls='None', mfc='None', ecolor = 'darkorange', mec='darkorange',marker='o', label="Driver+18")

    zdnW22, zupW22, rhoW22, errrhoW22u, errrhoW22dn = common.load_observation(obsdir, 'Global/Weaver22_SMD.dat', [0,1,2,3,4])
    rhoW22_cosmocorr = rhoW22 * hobs/h0 * 1e7
    zW22 = (zupW22 + zdnW22)/2.0
    errrhoW22u = np.log10(rhoW22 + errrhoW22u) - np.log10(rhoW22)
    errrhoW22dn = np.log10(rhoW22) - np.log10(rhoW22 - errrhoW22dn)
    ax.errorbar(zW22, np.log10(rhoW22_cosmocorr), xerr=[zW22 - zdnW22, zupW22 - zW22], yerr=[errrhoW22dn, errrhoW22u], ls='None', ecolor = 'navy',  mec='navy', marker='*', label='Weaver+22')

    zdnS23, zmidS23, zupS23, rhoS23, rhoS23_dn_s, rhoS23_up_s, rhoS23_dn_l, rhoS23_up_l = common.load_observation(obsdir, 'Global/Santini23_JWST_SMD.dat', [0, 1, 2, 3, 4, 5, 6, 7])
    rhoS23_cosmocorr = rhoS23 * hobs/h0
    rhoS23_errdn = np.log10(rhoS23) - np.log10(rhoS23_dn_s)
    rhoS23_errup = np.log10(rhoS23_up_s) - np.log10(rhoS23)
    rhoS23_errdn2 = np.log10(rhoS23) - np.log10(rhoS23_dn_l)
    rhoS23_errup2 = np.log10(rhoS23_up_l) - np.log10(rhoS23)

    ax.errorbar(zmidS23, np.log10(rhoS23_cosmocorr), xerr=[zmidS23 - zdnS23, zupS23 - zmidS23], yerr=[rhoS23_errdn2, rhoS23_errup2],  ls='None', ecolor = 'darkgreen',  mec='darkgreen', marker='s', label='Santini+23')

    common.prepare_legend(ax, ['k','b','r','k','grey', 'darkorange', 'navy', 'darkgreen'], loc=3)

    common.savefig(outdir, fig, "cosmic_smd_compL18.pdf")


def plot_sft_efficiency(plt, outdir, redshifts, sfre, sfreH2, mhrat):

    fig = plt.figure(figsize=(9.5,11))
    xmin, xmax, ymin, ymax = 0, 10, -6, 0

    # panel 1
    ax = fig.add_subplot(311)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(SFR/M_{\\rm cold} Gyr^{-1})$"
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(sfre > 0)
    ax.plot(redshifts[ind],np.log10(sfre[ind]), 'r', label='Shark')
    common.prepare_legend(ax, ['r'])

    #panel 2
    ax = fig.add_subplot(312)
    ytit="$\\rm log_{10}(SFE_{\\rm H_2}/Gyr^{-1})$"
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(sfreH2 > 0)
    ax.plot(redshifts[ind],np.log10(sfreH2[ind]), 'r', label='Shark')
    common.prepare_legend(ax, ['r'])

    #panel 3
    ax = fig.add_subplot(313)
    ytit="$\\rm log_{10}(M_{\\rm mol}/M_{\\rm atom})$"
    common.prepare_ax(ax, xmin, xmax, -3, 2, xtit, ytit, locators=(0.1, 1, 0.1))

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(mhrat > 0)
    ax.plot(redshifts[ind],np.log10(mhrat[ind]), 'r', label ='Shark')
    common.prepare_legend(ax, ['r'])

    common.savefig(outdir, fig, "cosmic_sfe.pdf")

    fig = plt.figure(figsize=(5,5))
    xmin, xmax, ymin, ymax = 0, 6, 0, 5

    # panel 1
    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    xtit="$\\rm redshift$"
    ytit="$\\rm \\Omega_{\\rm mol}/\\Omega_{\\rm atom}$"
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(mhrat > 0)
    ax.plot(redshifts[ind],mhrat[ind], 'r', label='Shark')
    common.prepare_legend(ax, ['r'],loc='upper left')
    common.savefig(outdir, fig, "cosmic_H2HIrat.pdf")


def plot_omega_h2(plt, outdir, obsdir, redshifts, h0, mH2den):


    def load_observations_h2(ax, obsdir, h0, caption=False):
        #Walter ASPECS ALMA program
        zloD16, zupD16, rhoH2loD16, rhoH2upD16  = common.load_observation(obsdir, 'Global/Decarli19_H2.dat', [0,1,2,3])
        zD16 =(zupD16 + zloD16)/2.0
        rhoH2D16 = (rhoH2loD16 + rhoH2upD16)/2.0
        hobs = 0.7
        xobs    = zD16
        errxlow = zD16-zloD16
        errxup  = zupD16-zD16
        yobs = rhoH2D16 + np.log10(pow(hobs/h0,3.0))
        errylow = rhoH2D16 - rhoH2loD16
        erryup  = rhoH2upD16 - rhoH2D16
        ax.errorbar(xobs, yobs, xerr=[errxlow,errxup], yerr=[errylow,erryup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='d',label="Decarli+19" if caption == True else None)

        #COLDz
        zloD16, zupD16, rhoH2loD16, rhoH2D16, rhoH2upD16  = common.load_observation(obsdir, 'Global/Riechers19_H2.dat', [0,1,2,3,4])
        zD16 =(zupD16 + zloD16)/2.0
        hobs = 0.7
        xobs    = zD16
        errxlow = zD16-zloD16
        errxup  = zupD16-zD16
        yobs = np.log10(rhoH2D16) + np.log10(pow(hobs/h0,3.0))
        errylow = np.log10(rhoH2D16) - np.log10(rhoH2loD16)
        erryup  = np.log10(rhoH2upD16) - np.log10(rhoH2D16)
        ax.errorbar(xobs, yobs, xerr=[errxlow,errxup], yerr=[errylow,erryup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s',label="Riechers+19" if caption == True else None)

        #ALMACAL-CO
        zD16, zloD16, zupD16, rhoH2D16, rhoH2loD16, rhoH2upD16  = common.load_observation(obsdir, 'Global/ASPECTS_Cont_H2.dat', [0,1,2,3,4,5])
        hobs = 0.7
        xobs    = zD16
        errxlow = zD16-zloD16
        errxup  = zupD16-zD16
        yobs = rhoH2D16 + np.log10(pow(hobs/h0,3.0))
        errylow = rhoH2D16 - rhoH2loD16
        erryup  = rhoH2upD16 - rhoH2D16
        ax.errorbar(xobs, yobs, xerr=[errxlow,errxup], yerr=[errylow,erryup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D',label="Decarli+20" if caption == True else None)


        #ALMACAL-CO
        zloD16, zupD16, rhoH2loD16, rhoH2upD16  = common.load_observation(obsdir, 'Global/Hamanowicz22_H2.dat', [0,1,2,3])
        zD16 =(zupD16 + zloD16)/2.0
        medrhoH2 = (rhoH2loD16 + rhoH2upD16)/2.0
        rhoH2D16 = np.log10(medrhoH2)
        hobs = 0.7
        xobs    = zD16
        errxlow = zD16-zloD16
        errxup  = zupD16-zD16
        yobs = rhoH2D16 + np.log10(pow(hobs/h0,3.0))
        errylow = np.log10(medrhoH2) - np.log10(rhoH2loD16)
        erryup  = np.log10(rhoH2upD16) - np.log10(medrhoH2)
        ax.errorbar(xobs, yobs, xerr=[errxlow,errxup], yerr=[errylow,erryup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='^',label="Hamanowicz+23" if caption == True else None)

        #z0 data
        zD16, zloD16, zupD16, rhoH2D16, rhoH2loD16, rhoH2upD16  = common.load_observation(obsdir, 'Global/H2_z0.dat', [0,1,2,3,4,5])
        xobs    = zD16
        errxlow = zD16-zloD16
        errxup  = zupD16-zD16
        yobs = np.log10(rhoH2D16) + np.log10(pow(hobs/h0,3.0))
        errylow = np.log10(rhoH2D16) - np.log10(rhoH2loD16)
        erryup  = np.log10(rhoH2upD16) - np.log10(rhoH2D16)
        ax.errorbar(xobs[0:1], yobs[0:1], xerr=[errxlow[0:1],errxup[0:1]], yerr=[errylow[0:1],erryup[0:1]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Boselli+14" if caption == True else None)
        ax.errorbar(xobs[1:2], yobs[1:2], xerr=[errxlow[1:2],errxup[1:2]], yerr=[errylow[1:2],erryup[1:2]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*',label="Fletcher+21" if caption == True else None)


    fig = plt.figure(figsize=(5,4.5))

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    xtit="$\\rm Lookback\, time/Gyr$"
    ytit="$\\rm log_{10}(\\rho_{\\rm H_2}/ M_{\odot}\,cMpc^{-3})$"
    common.prepare_ax(ax, 0, 13.5, 6.2, 8.4, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    new_tick_locations = np.array([0., 2., 4., 6., 8., 10., 12.])

    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(us.redshift(new_tick_locations), fontsize=12)

    ax2.set_xlabel("redshift",fontsize=13)

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(mH2den > 0)
    ax.plot(us.look_back_time(redshifts[ind]), np.log10(mH2den[ind]*pow(h0,2.0)) + np.log10(XH), 'r')
    print("will print Omega H2")
    for a,b,c in zip(redshifts, us.look_back_time(redshifts), np.log10(mH2den*pow(h0,2.0)) + np.log10(XH)):
        print(a,b,c)

    #print("H2 density of the Universe")
    #for a,b in zip(redshifts[ind],  np.log10(mH2den[ind]*pow(h0,2.0))):
    #    print(a,b)

    z, h2_modelvar = common.load_observation(obsdir, 'Models/SharkVariations/Global_OtherModels.dat', [0, 2])
    h2_modelvar_burst3 = h2_modelvar[0:179]
    h2_modelvar_nu0p5  = h2_modelvar[181:360]
    h2_modelvar_burst20= h2_modelvar[360:539]

    ind = np.where(h2_modelvar_burst20 > -10)
    ax.plot(us.look_back_time(z[ind]), h2_modelvar_burst20[ind], 'Sienna', linestyle='dotted')
    ind = np.where(h2_modelvar_burst3 > -10)
    ax.plot(us.look_back_time(z[ind]), h2_modelvar_burst3[ind], 'Crimson', linestyle='dashdot')
    ind = np.where(h2_modelvar_nu0p5 > -10)
    ax.plot(us.look_back_time(z[ind]), h2_modelvar_nu0p5[ind], 'Salmon', linestyle='dotted')

    load_observations_h2(ax, obsdir, h0, caption=True)

    # Legend
    common.prepare_legend(ax, ['grey','grey','grey'], loc=0)

    common.savefig(outdir, fig, "omega_H2.pdf")


    fig = plt.figure(figsize=(6,5.5))
    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(\\rho_{\\rm H_2}/ M_{\odot}\,cMpc^{-3})$"
    common.prepare_ax(ax, 0, 6, 5.2, 8.4, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(mH2den > 0)
    ax.plot(redshifts[ind], np.log10(mH2den[ind]*pow(h0,2.0)) + np.log10(XH), 'red', label='H$_2$ in gals (v2.0)')

    zin, omegaH2l18  = common.load_observation(obsdir, 'Models/SharkVariations/Global_OmegaGas_Lagos18.dat', [0,2])
    ax.plot(zin, omegaH2l18, 'k', linestyle='dashed', label = 'H$_2$ in gals (L18)')

    load_observations_h2(ax, obsdir, h0, caption=True)

    # Legend
    common.prepare_legend(ax, ['red','k','grey','grey','grey','grey','grey','grey'], loc=0)
    common.savefig(outdir, fig, "omega_H2_compL18.pdf")


def plot_mass_cosmic_density(plt, outdir, redshifts, mcold, mHI, mH2):

    fig = plt.figure(figsize=(5,4.5))

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    xtit="$\\rm Lookback\,time/Gyr$"
    ytit="$\\rm log_{10}(\\Omega_{\\rm gas})$"
    common.prepare_ax(ax, 0, 13.5, -4, -2.7, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    new_tick_locations = np.array([0., 2., 4., 6., 8., 10., 12.])

    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(us.redshift(new_tick_locations), fontsize=12)

    ax2.set_xlabel("redshift",fontsize=13)

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ax.plot(us.look_back_time(redshifts), mcold + np.log10(Omegab) - np.log10(XH), 'k', label='total neutral ISM')
    ax.plot(us.look_back_time(redshifts), mHI + np.log10(Omegab) - np.log10(XH), 'b', linestyle = 'dotted', label='atomic')
    ax.plot(us.look_back_time(redshifts), mH2 + np.log10(Omegab) - np.log10(XH), 'r', linestyle = 'dashed',label='molecular')

    common.prepare_legend(ax, ['k','b','r'], loc=1)
    common.savefig(outdir, fig, "omega_neutral.pdf")


def plot_mass_cosmic_density(plt, outdir, redshifts, mcold, mHI, mH2):

    fig = plt.figure(figsize=(5,4.5))

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    xtit="$\\rm Lookback\,time/Gyr$"
    ytit="$\\rm log_{10}(\\Omega_{\\rm gas})$"
    common.prepare_ax(ax, 0, 13.5, -4, -2.7, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    new_tick_locations = np.array([0., 2., 4., 6., 8., 10., 12.])

    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(us.redshift(new_tick_locations), fontsize=12)

    ax2.set_xlabel("redshift",fontsize=13)

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ax.plot(us.look_back_time(redshifts), mcold + np.log10(Omegab) - np.log10(XH), 'k', label='total neutral ISM')
    ax.plot(us.look_back_time(redshifts), mHI + np.log10(Omegab) - np.log10(XH), 'b', linestyle = 'dotted', label='atomic')
    ax.plot(us.look_back_time(redshifts), mH2 + np.log10(Omegab) - np.log10(XH), 'r', linestyle = 'dashed',label='molecular')

    common.prepare_legend(ax, ['k','b','r'], loc=1)
    common.savefig(outdir, fig, "omega_neutral.pdf")


def plot_cosmic_dust(plt, outdir, obsdir, redshifts, h0, mdustden, mdustden_mol):

    fig = plt.figure(figsize=(5,4.5))
    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    xtit="$\\rm Lookback\,time/Gyr$"
    ytit="$\\rm log_{10}(\\rho_{\\rm dust}/ M_{\odot}\,cMpc^{-3})$"
    common.prepare_ax(ax, 0, 13.5, 4, 6.3, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    new_tick_locations = np.array([0., 2., 4., 6., 8., 10., 12.])

    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(us.redshift(new_tick_locations), fontsize=12)

    ax2.set_xlabel("redshift",fontsize=13)

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(mdustden > 0)
    ax.plot(us.look_back_time(redshifts[ind]),np.log10(mdustden[ind]*pow(h0,2.0)),'r', label ='Shark all metals')
    
    ind = np.where(mdustden_mol > 0)
    ax.plot(us.look_back_time(redshifts[ind]),np.log10(mdustden_mol[ind]*pow(h0,2.0)),'r', linestyle = 'dashed', label ='Shark metals in molecular gas')
    
    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    redD17d,redD17u,smdD17,err1,err2,err3,err4 = common.load_observation(obsdir, 'Global/Driver18_dust.dat', [1,2,3,4,5,6,7])

    hobs = 0.7
    xobs = (redD17d+redD17u)/2.0
    yobs = smdD17 + np.log10(hobs/h0)

    err = yobs*0. - 999.
    err = np.sqrt(pow(err1,2.0)+pow(err2,2.0)+pow(err3,2.0)+pow(err4,2.0))

    ax.errorbar(us.look_back_time(xobs), yobs, yerr=[err,err], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Driver+18")

    common.prepare_legend(ax, ['r','r','grey'])
    common.savefig(outdir, fig, "cosmic_dust.pdf")


def plot_omega_HI(plt, outdir, obsdir, redshifts, h0, omegaHI, mcold):

    fig = plt.figure(figsize=(5,4.5))

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    xtit="$\\rm Lookback\, time/Gyr$"
    ytit="$\\rm log_{10}(\\Omega_{\\rm H_I})$"
    common.prepare_ax(ax, 0, 13.5, -5, -2, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    new_tick_locations = np.array([0., 2., 4., 6., 8., 10., 12.])

    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(us.redshift(new_tick_locations), fontsize=12)

    ax2.set_xlabel("redshift",fontsize=13)

    # note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(omegaHI > 0)
    ax.plot(us.look_back_time(redshifts[ind]), np.log10(omegaHI[ind]*pow(h0,2.0)) + np.log10(XH), 'r', label='Shark')
    print("will print Omega HI")
    for a,b,c in zip(redshifts, us.look_back_time(redshifts), np.log10(omegaHI*pow(h0,2.0)) + np.log10(XH)):
        print(a,b,c)

    z, hi_modelvar = common.load_observation(obsdir, 'Models/SharkVariations/Global_OtherModels.dat', [0, 1])
    hi_modelvar_burst3 = hi_modelvar[0:179]
    hi_modelvar_nu0p5  = hi_modelvar[181:360]
    hi_modelvar_burst20= hi_modelvar[360:539]

    ind = np.where(hi_modelvar_burst20 > -10)
    ax.plot(us.look_back_time(z[ind]), hi_modelvar_burst20[ind], 'Sienna', linestyle='dotted', label ='$\\eta_{\\rm burst}=20$')
    ind = np.where(hi_modelvar_burst3 > -10)
    ax.plot(us.look_back_time(z[ind]), hi_modelvar_burst3[ind], 'Crimson', linestyle='dashdot', label ='$\\eta_{\\rm burst}=3$')
    ind = np.where(hi_modelvar_nu0p5 > -10)
    ax.plot(us.look_back_time(z[ind]), hi_modelvar_nu0p5[ind], 'Salmon', linestyle='dotted', label ='$\\nu_{\\rm SF}=0.5 \\rm Gyr^{-1}$')

    xcgm = np.zeros(shape = 2)
    ycgm = np.zeros(shape = 2)

    xcgm[:] = 2.0
    ycgm[0] = -5.0
    ycgm[1] = -2.0
    ax.plot(us.look_back_time(xcgm),ycgm, 'k', linestyle='dotted', linewidth=0.85)
    ax.arrow(us.look_back_time(2.0), -2.5, 0.75, 0, head_width=0.05, head_length=0.1, fc='k', ec='k')
    ax.text(10.55, -2.4, 'CGM?', fontsize=12)

    # Rhee+18 compilation
    redR18,reddR18,reduR18,omegaR18,errdnR18,errupR18 = common.load_observation(obsdir, 'Global/HI_density_Rhee18.dat', [1,2,3,7,8,9])

    ax.errorbar(us.look_back_time(redR18),np.log10(omegaR18*1e-3), xerr=[reddR18,reduR18],yerr=[errdnR18,errupR18], ls='None', mfc='None', ecolor = 'grey', mec='grey', marker='o', label="Rhee+18 (comp)")


    common.prepare_legend(ax, ['r','Sienna','Crimson','Salmon','grey'])
    common.savefig(outdir, fig, "omega_HI.pdf")





    fig = plt.figure(figsize=(6,5.5))
    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(\\Omega_{\\rm H_I})$"
    common.prepare_ax(ax, 0, 6, -5, -2, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    # note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(omegaHI > 0)
    ax.plot(redshifts[ind], np.log10(omegaHI[ind]*pow(h0,2.0)) + np.log10(XH), 'red', label='HI in gals (v2.0)')
    ax.plot(redshifts[ind], mcold + np.log10(Omegab) - np.log10(XH), 'red', linestyle = 'dotted', label='HI+H$_{2}$ in gals (v2.0)')

    zin, omegaHIl18  = common.load_observation(obsdir, 'Models/SharkVariations/Global_OmegaGas_Lagos18.dat', [0,3])
    ax.plot(zin, omegaHIl18, 'k', linestyle='dashed', label = 'HI in gals (L18)')

    xcgm = np.zeros(shape = 2)
    ycgm = np.zeros(shape = 2)
    xcgm[:] = 2.0
    ycgm[0] = -5.0
    ycgm[1] = -2.0
    ax.plot(xcgm,ycgm, 'k', linestyle='dotted', linewidth=0.85)
    ax.arrow(2.0, -2.5, 0.75, 0, head_width=0.05, head_length=0.1, fc='k', ec='k')
    ax.text(2.1, -2.4, 'CGM?', fontsize=12)


    # Rhee+18 compilation
    redR18,reddR18,reduR18,omegaR18,errdnR18,errupR18 = common.load_observation(obsdir, 'Global/HI_density_Rhee18.dat', [1,2,3,7,8,9])
    ax.errorbar(redR18,np.log10(omegaR18*1e-3), xerr=[reddR18,reduR18],yerr=[errdnR18,errupR18], ls='None', mfc='None', ecolor = 'grey', mec='grey', marker='o', label="Rhee+18 (comp)")

    # Heintz+23 
    redK23, rederrK23, rhoHIK23, rhoHIupK23, rhoHIdnK23 = common.load_observation(obsdir, 'Global/HI_density_Heintz23.dat', [0,1,2,3, 4])
    hK13 = h0 * np.sqrt(OmegaM*(1.0+redK23)**3.0 + OmegaL)
    ax.errorbar(redK23,np.log10(rhoHIK23/((rho_crit*pow(h0,2.0)))), xerr=rederrK23, yerr=[np.log10(rhoHIK23)-np.log10(rhoHIdnK23), np.log10(rhoHIupK23)-np.log10(rhoHIK23)], ls='None', mfc='None', ecolor = 'darkgreen', mec='darkgreen', marker='s', label="Heintz+22 (HI in gals)")

    common.prepare_legend(ax, ['red','red', 'k','grey','darkgreen'])
    common.savefig(outdir, fig, "omega_HI_compL18.pdf")


def main(modeldir, outdir, redshift_table, subvols, obsdir):

    plt = common.load_matplotlib()
    fields = {'global': ('redshifts', 'm_hi', 'm_h2', 'mcold', 'mcold_metals',
                         'mhot_halo', 'mejected_halo', 'mbar_lost', 'mbar_created', 'mstars', 'mstars_bursts_mergers', 'mstars_bursts_diskinstabilities',
                         'm_bh', 'sfr_quiescent', 'sfr_burst', 'm_dm', 'mcold_halo', 'number_major_mergers', 
                         'number_minor_mergers', 'number_disk_instabilities', 'smbh_maximum')}

    # Read data from each subvolume at a time and add it up
    # rather than appending it all together
    for idx, subvol in enumerate(subvols):
        subvol_data = common.read_data(modeldir, redshift_table[0], fields, [subvol])
        max_bhs_subvol = subvol_data[20].copy()
        if idx == 0:
            hdf5_data        = subvol_data
            max_smbh         = max_bhs_subvol
        else:
            max_smbh = np.maximum(max_smbh, max_bhs_subvol)
            for subvol_datum, hdf5_datum in zip(subvol_data[3:], hdf5_data[3:]):
                hdf5_datum += subvol_datum
                #select the most massive black hole from the last list item

    # Also make sure that the total volume takes into account the number of subvolumes read
    hdf5_data[1] = hdf5_data[1] * len(subvols)

    h0, redshifts = hdf5_data[0], hdf5_data[2]

    (mstar_plot, mcold_plot, mhot_plot, meje_plot,
     mstar_dm_plot, mcold_dm_plot, mhot_dm_plot, meje_dm_plot, mbar_dm_plot,
     sfr, sfrd, sfrb, mstarden, mstarbden_mergers, mstarbden_diskins, sfre, sfreH2, mhrat,
     mHI_plot, mH2_plot, mH2den, mdustden, omegaHI, mdustden_mol, mcoldden, mhotden, 
     mejeden, history_interactions, mDMden, mlost_dm_plot, mcreated_dm_plot) = prepare_data(hdf5_data, redshifts)

    plot_mass_densities(plt, outdir, obsdir, h0, redshifts, mstar_plot, mcold_plot, mhot_plot, meje_plot, mstarden, mcoldden, mhotden, mejeden)
    plot_baryon_fractions(plt, outdir, redshifts, mstar_dm_plot, mcold_dm_plot, mhot_dm_plot, meje_dm_plot, mbar_dm_plot, mlost_dm_plot, mcreated_dm_plot)
    plot_cosmic_sfr(plt, outdir, obsdir, redshifts, h0, sfr, sfrd, sfrb, history_interactions, mDMden)
    plot_stellar_mass_cosmic_density(plt, outdir, obsdir, redshifts, h0, mstarden, mstarbden_mergers, mstarbden_diskins)
    plot_sft_efficiency(plt, outdir, redshifts, sfre, sfreH2, mhrat)
    plot_mass_cosmic_density(plt, outdir, redshifts, mcold_plot, mHI_plot, mH2_plot)
    plot_omega_h2(plt, outdir, obsdir, redshifts, h0, mH2den)
    plot_cosmic_dust(plt, outdir, obsdir, redshifts, h0, mdustden, mdustden_mol)
    plot_omega_HI(plt, outdir, obsdir, redshifts, h0, omegaHI, mcold_plot)

if __name__ == '__main__':
    main(*common.parse_args())
