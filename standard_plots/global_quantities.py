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

    (h0, volh, _, mHI, mH2, mcold, mcold_metals, mhot, meje, mstar,
     mstar_burst_mergers, mstar_burst_diskins, mBH, sfrdisk, sfrburst, mDM, mcold_halo) = hdf5_data

    sfrall = sfrdisk + sfrburst

    maxden = 6.2863*pow(10.0,9.0)
    ind = 0
    for i,j,k,l in zip(mhot, meje, mstar, mcold):
        tot = i+j+k+l
        if(tot/(volh/pow(h0,2.0)) > maxden):
            print("density exceeding maxden", tot/(volh/pow(h0,2.0))/maxden)
            #deltam = tot - maxden * (volh/pow(h0,2.0))
            #print "change",i, i - deltam
            #mhot[ind] = i - deltam
        ind = ind +1


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
    mcoldden  = mcold / volh
    mhotden   = mhot / volh
    mejeden   = meje / volh

    mstarbden_mergers = mstar_burst_mergers / volh
    mstarbden_diskins = mstar_burst_diskins / volh

    mH2den    = mH2 / volh


    mHIden   = mHI / volh
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
     sfr, sfrd, sfrb, mstarden, mstarbden_mergers, mstarbden_diskins, sfre, sfreH2, mhrat,
     mHI_plot, mH2_plot, mH2den, mdustden, omegaHI, mdustden_mol, mcoldden, mhotden, mejeden)

def plot_mass_densities(plt, outdir, obsdir, h0, redshifts, mstar, mcold, mhot, meje, mstarden, mcoldden, mhotden, mejeden):

    fig = plt.figure(figsize=(5,15))

    xtit="$\\rm Lookback\, time/Gyr$"
    ytit="$\\rm log_{10}(m/m_{\\rm bar,total})$"

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

    ind = np.where(mstarden > 0)
    ax.plot(us.look_back_time(redshifts), np.log10(mstarden[ind]*pow(h0,2.0)),'k', label='Shark')
    ind = np.where(mcoldden > 0)
    ax.plot(us.look_back_time(redshifts), np.log10(mcoldden[ind]*pow(h0,2.0)),'b')
    ind = np.where(mhotden > 0)
    ax.plot(us.look_back_time(redshifts), np.log10(mhotden[ind]*pow(h0,2.0)),'r')
    ind = np.where(mejeden > 0)
    ax.plot(us.look_back_time(redshifts), np.log10(mejeden[ind]*pow(h0,2.0)),'g')

    lbt, eaglesm, eaglesmout, eagleism, eaglehg, eagleejec = common.load_observation(obsdir + 'Models/OtherModels/', 'EAGLE_BaryonGrowthTotal.txt', [0,2,3,4,5,6])
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

#
#    lbt, galsm, galism, galhg = common.load_observation(obsdir, 'global_Lacey16', [1,5,4,3])
#    lbt = max(lbt) - lbt
#    ind = np.where(galsm > 0)
#    ax.plot(lbt[ind], np.log10(galsm[ind]) + 10, 'k', linestyle='dashed', label ='GALFORM')
#    ind = np.where(galism > 0)
#    ax.plot(lbt[ind], np.log10(galism[ind]) + 10, 'b', linestyle='dashed')
#    ind = np.where(galhg > 0)
#    ax.plot(lbt[ind], np.log10(galhg[ind]) + 10, 'r', linestyle='dashed')

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

    lbt, galsm, galism, galhg, galejec = common.load_observation(obsdir + 'Models/OtherModels/', 'global_Mitchell18', [0,2,3,4,5])
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

    #lbt, eaglesm, eagleism, eaglehg, eagleejec = common.load_observation(obsdir, 'for_claudia_galform.txt', [0,2,3,4,5])
    #eagletot = np.log10(pow(10.0, eaglesm) + pow(10.0, eagleism) +  pow(10.0, eaglehg) + pow(10.0, eagleejec))
    #ind = np.where(eaglesm > 0)
    #ax.plot(lbt[ind], eaglesm[ind] - 6.0 - np.log10(h0), 'k', linestyle='dotted', label ='GALFORM Peter')
    #ind = np.where(eagleism > 0)
    #ax.plot(lbt[ind], eagleism[ind]- 6.0 - np.log10(h0), 'b', linestyle='dotted')
    #ind = np.where(eaglehg > 0)
    #ax.plot(lbt[ind], eaglehg[ind] - 6.0 - np.log10(h0), 'r', linestyle='dotted')
    #ind = np.where(eagleejec > 0)
    #ax.plot(lbt[ind], eagleejec[ind] - 6.0 - np.log10(h0), 'g', linestyle='dotted')

#    lbt, galsm, galism, galhg = common.load_observation(obsdir, 'global_Lacey16', [1,5,4,3])
#    lbt = max(lbt) - lbt
#    ind = np.where(galsm > 0)
#    ax.plot(lbt[ind], np.log10(galsm[ind]) + 10, 'k', linestyle='dashed', label ='GALFORM')
#    ind = np.where(galism > 0)
#    ax.plot(lbt[ind], np.log10(galism[ind]) + 10, 'b', linestyle='dashed')
#    ind = np.where(galhg > 0)
#    ax.plot(lbt[ind], np.log10(galhg[ind]) + 10, 'r', linestyle='dashed')

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

    zbt, galism, galsm, galhg, galejec = common.load_observation(obsdir + 'Models/OtherModels/', 'BaryonBudgetLgalaxies.dat', [1,2,3,4,5])
    lbt = us.look_back_time(zbt)
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

    #lbt, eaglesm, eagleism, eaglehg, eagleejec = common.load_observation('/group/pawsey0119/clagos/Data/', 'for_claudia_galform.txt', [0,2,3,4,5])
    #eagletot = np.log10(pow(10.0, eaglesm) + pow(10.0, eagleism) +  pow(10.0, eaglehg) + pow(10.0, eagleejec))
    #ind = np.where(eaglesm > 0)
    #ax.plot(lbt[ind], eaglesm[ind] - 6.0 - np.log10(h0), 'k', linestyle='dotted', label ='GALFORM Peter')
    #ind = np.where(eagleism > 0)
    #ax.plot(lbt[ind], eagleism[ind]- 6.0 - np.log10(h0), 'b', linestyle='dotted')
    #ind = np.where(eaglehg > 0)
    #ax.plot(lbt[ind], eaglehg[ind] - 6.0 - np.log10(h0), 'r', linestyle='dotted')
    #ind = np.where(eagleejec > 0)
    #ax.plot(lbt[ind], eagleejec[ind] - 6.0 - np.log10(h0), 'g', linestyle='dotted')

#    lbt, galsm, galism, galhg = common.load_observation('/group/pawsey0119/clagos/Data/', 'global_Lacey16', [1,5,4,3])
#    lbt = max(lbt) - lbt
#    ind = np.where(galsm > 0)
#    ax.plot(lbt[ind], np.log10(galsm[ind]) + 10, 'k', linestyle='dashed', label ='GALFORM')
#    ind = np.where(galism > 0)
#    ax.plot(lbt[ind], np.log10(galism[ind]) + 10, 'b', linestyle='dashed')
#    ind = np.where(galhg > 0)
#    ax.plot(lbt[ind], np.log10(galhg[ind]) + 10, 'r', linestyle='dashed')

    common.prepare_legend(ax, ['k','k','k'], fontsize=10)


    common.savefig(outdir, fig, "global.pdf")


def plot_baryon_fractions(plt, outdir, redshifts, mstar_dm, mcold_dm, mhot_dm, meje_dm, mbar_dm):

    fig = plt.figure(figsize=(9.5,9.5))

    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(\\rho/\\rho_{\\rm m})$"

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, 0, 10, -4, 0.5, xtit, ytit, locators=(0.1, 1, 0.1))

    ax.plot(redshifts, mstar_dm, 'k', label='stars')
    ax.plot(redshifts, mcold_dm, 'b', label='ISM gas')
    ax.plot(redshifts, mhot_dm, 'r', label='halo gas')
    ax.plot(redshifts, meje_dm, 'g', label='ejected gas')
    ax.plot(redshifts, mbar_dm, 'm', label='total baryons')

    yplot = [0.1866920152, 0.1866920152]
    xplot = [0, 10]

    ax.plot(xplot, np.log10(yplot), 'k', linestyle='dashed', label ='Universal $f_{\\rm baryon}$')

    common.prepare_legend(ax, ['k','b','r','g','m','k'])
    common.savefig(outdir, fig, "baryon_frac.pdf")


def plot_cosmic_sfr(plt, outdir, obsdir, redshifts, h0, sfr, sfrd, sfrb):

    fig = plt.figure(figsize=(5,9))

    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(CSFRD/ M_{\odot}\,yr^{-1}\,cMpc^{-3})$"

    ax = fig.add_subplot(211)
    plt.subplots_adjust(left=0.15)

    common.prepare_ax(ax, 0, 10, -3, -0.5, xtit, ytit, locators=(0.1, 1, 0.1, 1))

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

    ax.errorbar(xobs[0:8], yobs[0:8], yerr=[yobs[0:8]-lerr[0:8],herr[0:8]-yobs[0:8]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D')
    ax.errorbar(xobs[9:17], yobs[9:17], yerr=[yobs[9:17]-lerr[9:17],herr[9:17]-yobs[9:17]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='x')

    #Driver (Chabrier IMF), ['Baldry+2012, z<0.06']
    redD17d, redD17u, sfrD17, err1, err2, err3, err4 = common.load_observation(obsdir, 'SFR/Driver17_sfr.dat', [0,1,2,3,4,5,6])

    hobs = 0.7
    xobsD17 = (redD17d+redD17u)/2.0
    yobsD17 = sfrD17 + np.log10(hobs/h0)

    errD17 = yobs*0. - 999.
    errD17 = np.sqrt(pow(err1,2.0)+pow(err2,2.0)+pow(err3,2.0)+pow(err4,2.0))
    ax.errorbar(xobsD17, yobsD17, yerr=[errD17,errD17], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

    #note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
    ind = np.where(sfr > 0)
    ax.plot(redshifts[ind], np.log10(sfr[ind]*pow(h0,2.0)), 'k', linewidth=1, label ='total')
    ind = np.where(sfrd > 0)
    ax.plot(redshifts[ind], np.log10(sfrd[ind]*pow(h0,2.0)), 'b', linestyle='dashed', linewidth=1, label ='quiescent')
    ind = np.where(sfrb > 0)
    ax.plot(redshifts[ind], np.log10(sfrb[ind]*pow(h0,2.0)),'r', linestyle='dotted',  linewidth=1, label ='bursts')

    sfr_modelvar = common.load_observation(obsdir + 'Models/SharkVariations', 'Global_OtherModels.dat', [2])
    sfr_modelvar_burst1 = sfr_modelvar[0:179]
    sfr_modelvar_nu0p5  = sfr_modelvar[181:360]

    ind = np.where(sfr_modelvar_burst1 > -10)
    ax.plot(redshifts[ind], sfr_modelvar_burst1[ind], 'DarkSlateGray', linestyle='dashdot', label ='$\\eta_{\\rm burst}=1$')
    ind = np.where(sfr_modelvar_nu0p5 > -10)
    ax.plot(redshifts[ind], sfr_modelvar_nu0p5[ind], 'SlateGray', linestyle='dotted', label ='$\\nu_{\\rm SF}=0.5 \\rm Gyr^{-1}$')

    common.prepare_legend(ax, ['k','b','r','DarkSlateGray','SlateGray','grey','grey','grey'], loc=1)

    xtit="$\\rm Lookback\, time/Gyr$"
    ax = fig.add_subplot(212)
    plt.subplots_adjust(left=0.15)

    common.prepare_ax(ax, 0, 13.5, -3, -0.5, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.errorbar(us.look_back_time(xobs[0:8]), yobs[0:8], yerr=[yobs[0:8]-lerr[0:8],herr[0:8]-yobs[0:8]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D',label="Karim+11 obs")
    ax.errorbar(us.look_back_time(xobs[9:17]), yobs[9:17], yerr=[yobs[9:17]-lerr[9:17],herr[9:17]-yobs[9:17]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='x', label="Karim+11 extr")
    ax.errorbar(us.look_back_time(xobsD17), yobsD17, yerr=[errD17,errD17], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o', label="Driver+18")

    ind = np.where(sfr > 0)
    ax.plot(us.look_back_time(redshifts[ind]), np.log10(sfr[ind]*pow(h0,2.0)), 'k',  linewidth=1)
    ind = np.where(sfrd > 0)
    ax.plot(us.look_back_time(redshifts[ind]), np.log10(sfrd[ind]*pow(h0,2.0)), 'b', linestyle='dashed', linewidth=1)
    ind = np.where(sfrb > 0)
    ax.plot(us.look_back_time(redshifts[ind]), np.log10(sfrb[ind]*pow(h0,2.0)),'r', linestyle='dotted',  linewidth=1)

    ind = np.where(sfr_modelvar_burst1 > -10)
    ax.plot(us.look_back_time(redshifts[ind]), sfr_modelvar_burst1[ind], 'DarkSlateGray', linestyle='dashdot')
    ind = np.where(sfr_modelvar_nu0p5 > -10)
    ax.plot(us.look_back_time(redshifts[ind]), sfr_modelvar_nu0p5[ind], 'SlateGray', linestyle='dotted')

    common.prepare_legend(ax, ['grey','grey','grey'], loc=2)


    common.savefig(outdir, fig, "cosmic_sfr.pdf")


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
    ind = np.where(mstarbden_mergers > 0)
    ax.plot(redshifts[ind],np.log10(mstarbden_mergers[ind]*pow(h0,2.0)), 'r', linestyle='dashed')
    ind = np.where(mstarbden_diskins > 0)
    ax.plot(redshifts[ind],np.log10(mstarbden_diskins[ind]*pow(h0,2.0)), 'b', linestyle='dotted')


    sm_modelvar = common.load_observation(obsdir + 'Models/SharkVariations', 'Global_OtherModels.dat', [3])
    sm_modelvar_burst1 = sm_modelvar[0:179]
    sm_modelvar_nu0p5  = sm_modelvar[181:360]

    ind = np.where(sm_modelvar_burst1 > -10)
    ax.plot(redshifts[ind], sm_modelvar_burst1[ind], 'DarkSlateGray', linestyle='dashdot')
    ind = np.where(sm_modelvar_nu0p5 > -10)
    ax.plot(redshifts[ind], sm_modelvar_nu0p5[ind], 'SlateGray', linestyle='dotted')

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    redD17d, redD17u, smdD17, err1, err2, err3, err4 = common.load_observation(obsdir, 'SFR/Driver17_smd.dat', [1,2,3,4,5,6,7])

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
    ind = np.where(mstarbden_mergers > 0)
    ax.plot(us.look_back_time(redshifts[ind]),np.log10(mstarbden_mergers[ind]*pow(h0,2.0)), 'r', linestyle='dashed', label='formed in galaxy mergers')
    ind = np.where(mstarbden_diskins > 0)
    ax.plot(us.look_back_time(redshifts[ind]),np.log10(mstarbden_diskins[ind]*pow(h0,2.0)), 'b', linestyle='dotted', label='formed in disk instabilities')

    ind = np.where(sm_modelvar_burst1 > -10)
    ax.plot(us.look_back_time(redshifts[ind]), sm_modelvar_burst1[ind], 'DarkSlateGray', linestyle='dashdot', label ='$\\eta_{\\rm burst}=1$')
    ind = np.where(sm_modelvar_nu0p5 > -10)
    ax.plot(us.look_back_time(redshifts[ind]), sm_modelvar_nu0p5[ind], 'SlateGray', linestyle='dotted', label ='$\\nu_{\\rm SF}=0.5 \\rm Gyr^{-1}$')

    ax.errorbar(us.look_back_time(xobs), yobs, yerr=[err,err], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

    common.prepare_legend(ax, ['k','r','b', 'DarkSlateGray','SlateGray'], loc=3)

    common.savefig(outdir, fig, "cosmic_smd.pdf")


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


def plot_omega_h2(plt, outdir, obsdir, redshifts, h0, mH2den):

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
    ax.plot(us.look_back_time(redshifts[ind]), np.log10(mH2den[ind]*pow(h0,2.0)) + np.log10(XH), 'r', label='Shark')

    h2_modelvar = common.load_observation(obsdir + 'Models/SharkVariations', 'Global_OtherModels.dat', [1])
    h2_modelvar_burst1 = h2_modelvar[0:179]
    h2_modelvar_nu0p5  = h2_modelvar[181:360]

    ind = np.where(h2_modelvar_burst1 > -10)
    ax.plot(us.look_back_time(redshifts[ind]), h2_modelvar_burst1[ind], 'Crimson', linestyle='dashdot', label ='$\\eta_{\\rm burst}=1$')
    ind = np.where(h2_modelvar_nu0p5 > -10)
    ax.plot(us.look_back_time(redshifts[ind]), h2_modelvar_nu0p5[ind], 'Salmon', linestyle='dotted', label ='$\\nu_{\\rm SF}=0.5 \\rm Gyr^{-1}$')

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    zD16, zloD16, zupD16, rhoH2D16, rhoH2loD16, rhoH2upD16  = common.load_observation(obsdir, 'SFR/Walter17_H2.dat', [0,1,2,3,4,5])

    hobs = 0.7

    xobs    = zD16
    errxlow = zD16-zloD16
    errxup  = zupD16-zD16
    yobs = np.log10(rhoH2D16) + np.log10(pow(hobs/h0,3.0))
    errylow = np.log10(rhoH2D16) - np.log10(rhoH2loD16)
    erryup  = np.log10(rhoH2upD16) - np.log10(rhoH2D16)

    ax.errorbar(us.look_back_time(xobs), yobs, xerr=[errxlow,errxup], yerr=[errylow,erryup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='+',label="Decarli+16")
    ax.errorbar(us.look_back_time(xobs[0:1]), yobs[0:1], xerr=[errxlow[0:1],errxup[0:1]], yerr=[errylow[0:1],erryup[0:1]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Boselli+14")

    # Legend
    common.prepare_legend(ax, ['r','Crimson','Salmon','grey','grey','grey'], loc=0)

    common.savefig(outdir, fig, "omega_H2.pdf")


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
    ax.plot(us.look_back_time(redshifts), mcold + np.log10(Omegab), 'k', label='total neutral ISM')
    ax.plot(us.look_back_time(redshifts), mHI + np.log10(Omegab), 'b', linestyle = 'dotted', label='atomic')
    ax.plot(us.look_back_time(redshifts), mH2 + np.log10(Omegab), 'r', linestyle = 'dashed',label='molecular')

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
    redD17d,redD17u,smdD17,err1,err2,err3,err4 = common.load_observation(obsdir, 'SFR/Driver17_dust.dat', [1,2,3,4,5,6,7])

    hobs = 0.7
    xobs = (redD17d+redD17u)/2.0
    yobs = smdD17 + np.log10(hobs/h0)

    err = yobs*0. - 999.
    err = np.sqrt(pow(err1,2.0)+pow(err2,2.0)+pow(err3,2.0)+pow(err4,2.0))

    ax.errorbar(us.look_back_time(xobs), yobs, yerr=[err,err], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Driver+18")

    common.prepare_legend(ax, ['r','r','grey'])
    common.savefig(outdir, fig, "cosmic_dust.pdf")


def plot_omega_HI(plt, outdir, obsdir, redshifts, h0, omegaHI):

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

    hi_modelvar = common.load_observation(obsdir + 'Models/SharkVariations', 'Global_OtherModels.dat', [0])
    hi_modelvar_burst1 = hi_modelvar[0:179]
    hi_modelvar_nu0p5  = hi_modelvar[181:360]

    ind = np.where(hi_modelvar_burst1 > -10)
    ax.plot(us.look_back_time(redshifts[ind]), hi_modelvar_burst1[ind], 'Crimson', linestyle='dashdot', label ='$\\eta_{\\rm burst}=1$')
    ind = np.where(hi_modelvar_nu0p5 > -10)
    ax.plot(us.look_back_time(redshifts[ind]), hi_modelvar_nu0p5[ind], 'Salmon', linestyle='dotted', label ='$\\nu_{\\rm SF}=0.5 \\rm Gyr^{-1}$')

    # Rhee+18 compilation
    redR18,reddR18,reduR18,omegaR18,errdnR18,errupR18 = common.load_observation(obsdir, 'Gas/HI_density_for_Claudia.dat', [1,2,3,7,8,9])

    ax.errorbar(us.look_back_time(redR18),np.log10(omegaR18*1e-3), xerr=[reddR18,reduR18],yerr=[errdnR18,errupR18], ls='None', mfc='None', ecolor = 'grey', mec='grey', marker='o', label="Rhee+18 (comp)")

    common.prepare_legend(ax, ['r','Crimson','Salmon','grey'])
    common.savefig(outdir, fig, "omega_HI.pdf")

def main(modeldir, outdir, subvols, obsdir, snapshot):

    plt = common.load_matplotlib()
    fields = {'global': ('redshifts', 'm_hi', 'm_h2', 'mcold', 'mcold_metals',
                         'mhot_halo', 'mejected_halo', 'mstars', 'mstars_bursts_mergers', 'mstars_bursts_diskinstabilities',
                         'm_bh', 'sfr_quiescent', 'sfr_burst', 'm_dm', 'mcold_halo')}

    # Read data from each subvolume at a time and add it up
    # rather than appending it all together
    for idx, subvol in enumerate(subvols):
        subvol_data = common.read_data(modeldir, snapshot, fields, [subvol])
        if idx == 0:
            hdf5_data = subvol_data
        else:
            for subvol_datum, hdf5_datum in zip(subvol_data[3:], hdf5_data[3:]):
                hdf5_datum += subvol_datum

    # Also make sure that the total volume takes into account the number of subvolumes read
    hdf5_data[1] = hdf5_data[1] * len(subvols)

    h0, redshifts = hdf5_data[0], hdf5_data[2]

    (mstar_plot, mcold_plot, mhot_plot, meje_plot,
     mstar_dm_plot, mcold_dm_plot, mhot_dm_plot, meje_dm_plot, mbar_dm_plot,
     sfr, sfrd, sfrb, mstarden, mstarbden_mergers, mstarbden_diskins, sfre, sfreH2, mhrat,
     mHI_plot, mH2_plot, mH2den, mdustden, omegaHI, mdustden_mol, mcoldden, mhotden, mejeden) = prepare_data(hdf5_data, redshifts)

    plot_mass_densities(plt, outdir, obsdir, h0, redshifts, mstar_plot, mcold_plot, mhot_plot, meje_plot, mstarden, mcoldden, mhotden, mejeden)
    plot_baryon_fractions(plt, outdir, redshifts, mstar_dm_plot, mcold_dm_plot, mhot_dm_plot, meje_dm_plot, mbar_dm_plot)
    plot_cosmic_sfr(plt, outdir, obsdir, redshifts, h0, sfr, sfrd, sfrb)
    plot_stellar_mass_cosmic_density(plt, outdir, obsdir, redshifts, h0, mstarden, mstarbden_mergers, mstarbden_diskins)
    plot_sft_efficiency(plt, outdir, redshifts, sfre, sfreH2, mhrat)
    plot_mass_cosmic_density(plt, outdir, redshifts, mcold_plot, mHI_plot, mH2_plot)
    plot_omega_h2(plt, outdir, obsdir, redshifts, h0, mH2den)
    plot_cosmic_dust(plt, outdir, obsdir, redshifts, h0, mdustden, mdustden_mol)
    plot_omega_HI(plt, outdir, obsdir, redshifts, h0, omegaHI)

if __name__ == '__main__':
    main(*common.parse_args())
