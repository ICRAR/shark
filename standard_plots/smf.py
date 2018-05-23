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
"""SMF plots"""

import collections
import functools
import math

import numpy as np

import common
import utilities_statistics as us


observation = collections.namedtuple('observation', 'label x y yerrup yerrdn err_absolute')

##################################
# Constants
GyrToYr = 1e9
Zsun = 0.0127
XH = 0.72

##################################
# Mass function initialization
mlow = 5
mupp = 14
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0
imf   = 'cha'


def plot_stellarmf_z(plt, outdir, obsdir, h0, plotz, hist_smf, hist_smf_cen, hist_smf_sat, hist_smf_err):

    # z0 obs
    z0obs = []
    lm, p, dpdn, dpup = common.load_observation(obsdir, 'mf/GAMAII_BBD_GSMFs.data', [0,1,2,3])
    xobs = lm
    indx = np.where(p > 0)
    yobs = np.log10(p[indx])
    ydn = yobs - np.log10(p[indx]-dpdn[indx])
    yup = np.log10(p[indx]+dpup[indx]) - yobs
    z0obs.append((observation("Wright+2017", xobs[indx], yobs, ydn, yup, err_absolute=False), 'o'))

    lm, p3, dp3 = common.load_observation(obsdir, 'mf/SMF_Bernardi2013_SerExp.data', [0,1,2])
    if imf == 'ken':
        lm = lm - 0.09
    xobs = lm
    yobs = np.full(xobs.shape, -999.)
    lerr = np.full(xobs.shape, -999.)
    herr = np.full(xobs.shape, 999.)
    indx = np.where( p3 < 1)
    yobs[indx] = p3[indx]
    lerr[indx] = dp3[indx]
    herr[indx] = dp3[indx]
    z0obs.append((observation("Bernardi+2013", xobs, yobs, lerr, herr, err_absolute=False), 'x'))

    # Moustakas (Chabrier IMF), ['Moustakas+2013, several redshifts']
    zdnM13, lmM13, pM13, dp_dn_M13, dp_up_M13 = common.load_observation(obsdir, 'mf/SMF_Moustakas2013.data', [0,3,5,6,7])
    if imf == 'ken':
        lmM13 = lmM13 - 0.09
    xobsM13 = lmM13

    yobsM13 = np.full(xobsM13.shape, -999.)
    lerrM13 = np.full(xobsM13.shape, -999.)
    herrM13 = np.full(xobsM13.shape, 999.)
    indx = np.where( pM13 < 1)
    yobsM13[indx] = (pM13[indx])
    indx = np.where( dp_dn_M13 > 0)
    lerrM13[indx]  = dp_dn_M13[indx] 
    indx = np.where( dp_up_M13 > 0)
    herrM13[indx]  = dp_up_M13[indx]

    # Muzzin (Chabrier IMF), ['Moustakas+2013, several redshifts']
    zdnMu13,zupMu13,lmMu13,pMu13,dp_dn_Mu13,dp_up_Mu13 = common.load_observation(obsdir, 'mf/SMF_Muzzin2013.data', [0,1,2,4,5,5])
    if(imf == 'ken'):
        lmMu13 = lmMu13 - 0.09
    xobsMu13 = lmMu13

    yobsMu13 = np.full(xobsMu13.shape, -999.)
    lerrMu13 = np.full(xobsMu13.shape, -999.)
    herrMu13 = np.full(xobsMu13.shape, 999.)
    indx = np.where( pMu13 < 1)
    yobsMu13[indx] = (pMu13[indx])
    indx = np.where( dp_dn_Mu13 > 0)
    lerrMu13[indx]  = dp_dn_Mu13[indx] 
    indx = np.where( dp_up_Mu13 > 0)
    herrMu13[indx]  = dp_up_Mu13[indx]

    # Santini 2012
    zdnS12, lmS12, pS12, dp_dn_S12, dp_up_S12 = common.load_observation(obsdir, 'mf/SMF_Santini2012.data', [0,2,3,4,5])
    hobs = 0.7
    xobsS12 = lmS12 - 0.24 +  np.log10(hobs/h0)
    yobsS12 = np.full(xobsS12.shape, -999.)
    lerrS12 = np.full(xobsS12.shape, -999.)
    herrS12 = np.full(xobsS12.shape, 999.)
    indx = np.where( pS12 < 1)
    yobsS12[indx] = (pS12[indx]) + np.log10(pow(h0/hobs,3.0))
    indx = np.where( dp_dn_S12 > 0)
    lerrS12[indx]  = dp_dn_S12[indx]
    indx = np.where( dp_up_S12 > 0)
    herrS12[indx]  = dp_up_S12[indx]

    # z0.5 obs
    z05obs = []
    in_redshift = np.where(zdnM13 == 0.4)
    z05obs.append((observation("Moustakas+2013", xobsM13[in_redshift], yobsM13[in_redshift], lerrM13[in_redshift], herrM13[in_redshift], err_absolute=False), 'o'))
    in_redshift = np.where(zdnMu13 == 0.5)
    z05obs.append((observation("Muzzin+2013", xobsMu13[in_redshift], yobsMu13[in_redshift], lerrMu13[in_redshift], herrMu13[in_redshift], err_absolute=False), '+'))

    # z1 obs
    z1obs = []
    in_redshift = np.where(zdnM13 == 0.8)
    z1obs.append((observation("Moustakas+2013", xobsM13[in_redshift], yobsM13[in_redshift], lerrM13[in_redshift], herrM13[in_redshift], err_absolute=False), 'o'))
    in_redshift = np.where(zdnMu13 == 1)
    z1obs.append((observation("Muzzin+2013", xobsMu13[in_redshift], yobsMu13[in_redshift], lerrMu13[in_redshift], herrMu13[in_redshift], err_absolute=False), '+'))

    #z2 obs
    z2obs = []
    in_redshift = np.where(zupMu13 == 2.5)
    z2obs.append((observation("Muzzin+2013", xobsMu13[in_redshift], yobsMu13[in_redshift], lerrMu13[in_redshift], herrMu13[in_redshift], err_absolute=False), '+'))
    in_redshift = np.where(zdnS12 == 1.8)
    z2obs.append((observation("Santini+2012", xobsS12[in_redshift], yobsS12[in_redshift], lerrS12[in_redshift], herrS12[in_redshift], err_absolute=False), 'o'))

    # z3 obs
    z3obs = []
    in_redshift = np.where(zupMu13 == 3.0)
    z3obs.append((observation("Muzzin+2013", xobsMu13[in_redshift], yobsMu13[in_redshift], lerrMu13[in_redshift], herrMu13[in_redshift], err_absolute=False), '+'))
    in_redshift = np.where(zdnS12 == 2.5)
    z3obs.append((observation("Santini+2012", xobsS12[in_redshift], yobsS12[in_redshift], lerrS12[in_redshift], herrS12[in_redshift], err_absolute=False), 'o'))

    # z4 obs
    z4obs = []
    in_redshift = np.where(zupMu13 == 4.0)
    z4obs.append((observation("Muzzin+2013", xobsMu13[in_redshift], yobsMu13[in_redshift], lerrMu13[in_redshift], herrMu13[in_redshift], err_absolute=False), '+'))
    in_redshift = np.where(zdnS12 == 3.5)
    z4obs.append((observation("Santini+2012", xobsS12[in_redshift], yobsS12[in_redshift], lerrS12[in_redshift], herrS12[in_redshift], err_absolute=False), 'o'))

    fig = plt.figure(figsize=(9.7,11.7))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm star}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\Phi/dlog{\\rm M_{\\rm star}}/{\\rm Mpc}^{-3} )$"
    xmin, xmax, ymin, ymax = 8, 13, -6, -1
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (321, 322, 323, 324, 325, 326)
    indeces = (0, 1, 2, 3, 4, 5)
    zs = (0, 0.5, 1, 2, 3, 4)
    observations = (z0obs, z05obs, z1obs, z2obs, z3obs, z4obs)
    for subplot, idx, z, obs_and_markers in zip(subplots, indeces, zs, observations):

        ax = fig.add_subplot(subplot)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
        ax.text(xleg, yleg, 'z=%s' % str(z))

        # Observations
        for obs, marker in obs_and_markers:
            common.errorbars(ax, obs.x, obs.y, obs.yerrdn, obs.yerrup, 'grey',
                             marker, err_absolute=obs.err_absolute, label=obs.label)

        # Predicted SMF
        if plotz[idx]:
            y = hist_smf[idx,:]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],'r', label='all galaxies' if idx == 0 else None)
            y = hist_smf_cen[idx,:]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],'b', linestyle='dotted', label ='centrals' if idx == 0 else None)
            y = hist_smf_sat[idx,:]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],'g', linestyle='dashed', label ='satellites' if idx == 0 else None)

            if z >= 3:
                y = hist_smf_err[4,:]
                ind = np.where(y < 0.)
                ax.plot(xmf[ind],y[ind],'r', linestyle='dashdot', linewidth=2, label ='0.25 dex error')

        colors = []
        if idx == 0:
            colors = ['r','b','g']
        if z >= 3:
            colors = ['r']
        colors += ['grey', 'grey']

        common.prepare_legend(ax, colors)

    common.savefig(outdir, fig, 'stellarmf_z.pdf')


def plot_HImf_z0(plt, outdir, obsdir, h0, plotz_HImf, hist_HImf, hist_HImf_cen, hist_HImf_sat):

    fig = plt.figure(figsize=(5,5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm HI}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\Phi/dlog{\\rm M_{\\rm HI}}/{\\rm Mpc}^{-3} )$"
    xmin, xmax, ymin, ymax = 7.1, 12, -6, 0
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg, yleg, 'z=0')

    #HIPASS
    lmHI, pHI, dpHIdn, dpHIup = common.load_observation(obsdir, 'mf/HIMF/HIMF_Zwaan2005.data', [0,1,2,3])

    #correct data for their choice of cosmology
    hobs = 0.75
    xobs = lmHI + np.log10(pow(hobs,2)/pow(h0,2))
    yobs = pHI + np.log10(pow(h0,3)/pow(hobs,3))
    ax.errorbar(xobs, yobs, yerr=[dpHIdn,dpHIup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Zwaan+2005")

    #ALFALFA.40
    lmHI, pHI, dpHI = common.load_observation(obsdir, 'mf/HIMF/HIMF_Martin2010.data', [0,1,2])

    #correct data for their choice of cosmology
    hobs = 0.7
    xobs = lmHI + np.log10(pow(hobs,2)/pow(h0,2))
    yobs = pHI + np.log10(pow(h0,3)/pow(hobs,3))

    ax.errorbar(xobs, yobs, yerr=[dpHI,dpHI], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='x',label="Martin+2010")

    # Predicted HIMF
    if plotz_HImf[0]:
        y = hist_HImf[0,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'r', label ='all galaxies') 
        y = hist_HImf_cen[0,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'b', linestyle='dotted',label ='centrals') 
        y = hist_HImf_sat[0,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'g', linestyle='dashed',label ='satellites') 

    common.prepare_legend(ax, ['r','b','g','grey'])
    common.savefig(outdir, fig, 'HImf_z0.pdf')


def plot_H2mf_z0(plt, outdir, obsdir, h0, plotz_HImf, hist_H2mf, hist_H2mf_cen, hist_H2mf_sat):

    fig = plt.figure(figsize=(5,5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm H_2}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\Phi/dlog{\\rm M_{\\rm H_2}}/{\\rm Mpc}^{-3} )$"
    xmin, xmax, ymin, ymax = 7.1, 12, -6, 0
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg, yleg, 'z=0')

    # H2 mass function
    lmCO, pCO, dpCOdn, dpCOup = common.load_observation(obsdir, 'Keres03_LCOLF60m.dat', [0,1,2,3])

    # correct data for their choice of cosmology
    hobs = 0.7
    dpCOdn = np.abs(dpCOdn-pCO)
    dpCOup = np.abs(dpCOup-pCO)
    xobs = lmCO + np.log10(pow(hobs,2)/pow(h0,2))
    yobs = pCO + np.log10(pow(h0,3)/pow(hobs,3))

    X = 2.0
    corr_fac_H2 = np.log10(580.*X)+2.*np.log10(2.6)-np.log10(4.*math.pi)
    ax.errorbar(xobs+corr_fac_H2, yobs, yerr=[dpCOdn,dpCOup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='^',label="Keres+03")

    #H2 mass function
    lm,p,dpdn,dpup = common.load_observation(obsdir, 'mf/B14_MH2MF.dat', [0,1,2,3])
    #correct data for their choice of cosmology
    logbin = np.log10(0.4)
    ax.errorbar(lm[2:8],p[2:8]-logbin,yerr=[p[2:8]-dpdn[2:8],dpup[2:8]-p[2:8]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s',label="Boselli+14")

    #Predicted H2MF
    if plotz_HImf[0]:
        y = hist_H2mf[0,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'r', label ='all galaxies') 
        y = hist_H2mf_cen[0,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'b', linestyle='dotted',label ='centrals') 
        y = hist_H2mf_sat[0,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'g', linestyle='dashed',label ='satellites') 

    common.prepare_legend(ax, ['r','b','g','grey','grey'])
    common.savefig(outdir, fig, 'H2mf_z0.pdf')


def plot_SSFR_Mstars(plt, outdir, mainseq, mainseq_cen, mainseq_sat):

    fig = plt.figure(figsize=(9.5,9.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm star}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\\rm SSFR/Gyr^{-1})$"
    xmin, xmax, ymin, ymax = 7.1, 12, -3, 2
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    z0obs = observation('Gilbank+10', [8,11.2], np.log10([0.3,0.02]), None, None, False)
    z1obs = observation('Karim+11', [9.6,11.1], np.log10([1,0.3]), None, None, False)
    z2obs = observation('Karim+11', [9.6,11.1], np.log10([3.5,1.5]), None, None, False)

    indeces = (0, 1, 2, 4)
    zs = (0, 0.5, 1, 2)
    subplots = (221, 222, 223, 224)
    observations = (z0obs, None, z1obs, z2obs)

    for subplot, z, idx, obs in zip(subplots, zs, indeces, observations):

        ax = fig.add_subplot(subplot)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
        ax.text(xleg, yleg, 'z=%s' % str(z))

        #Predicted relation
        ind = np.where(mainseq[idx,0,:] != 0)
        xplot = xmf[ind]
        yplot = mainseq[idx,0,ind]
        errdn = mainseq[idx,1,ind]
        errup = mainseq[idx,2,ind]

        ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'b', mec='b',marker='o',label="all galaxies")

        ind = np.where(mainseq_cen[idx,0,:] != 0)
        xplot = xmf[ind]
        yplot = mainseq_cen[idx,0,ind]
        errdn = mainseq_cen[idx,1,ind]
        errup = mainseq_cen[idx,2,ind]

        ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'g', mec='g',marker='o',markersize=2,label="centrals")
        ind = np.where(mainseq_sat[idx,0,:] != 0)
        xplot = xmf[ind]
        yplot = mainseq_sat[idx,0,ind]
        errdn = mainseq_sat[idx,1,ind]
        errup = mainseq_sat[idx,2,ind]

        ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'r', mec='r',marker='o',markersize=2,label="satellites")

        # observations approximate of Gilbank et al. (2010)
        if obs:
            ax.plot(obs.x, obs.y, 'k', label=obs.label)

        colors = ['b', 'g', 'r']
        if obs:
            colors.insert(0, 'k')
        common.prepare_legend(ax, colors)

    common.savefig(outdir, fig, 'SSFR_Mstars.pdf')



def plot_mzr(plt, outdir, obsdir, mzr, mzr_cen, mzr_sat):

    fig = plt.figure(figsize=(9.5,9.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm star}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\\rm Z_{\\rm gas}/Z_{\odot})$"
    xmin, xmax, ymin, ymax = 7.1, 12, -3, 1
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    # Observations at z = 0
    # Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    lm, mz, mzdn, mzup = common.load_observation(obsdir, 'MMAdrews13.dat', [0,1,2,3])
    corrzsun = 8.69
    z0obs = observation('Andrews & Martini (2013)', lm, mz - corrzsun, mzdn - corrzsun, mzup - corrzsun, True)

    # observations approximate of Maiolino et al. (2008)
    z05obs = observation('Maiolino+08', [8.25,10.5,11.5], [-0.59,0.26,0.36], None, None, None)
    z2obs = observation('Maiolino+08', [8.25,10.5,11.5], [-1.19,-0.04,0.21], None, None, None)

    indeces = (0, 1, 2, 3)
    subplots = (221, 222, 223, 224)
    zs = (0, 0.5, 1, 2)
    ytitles = (True, False, True, False)
    observations = (z0obs, z05obs, None, z2obs)

    for idx, subplot, z, ytitle, obs in zip(indeces, subplots, zs, ytitles, observations):

        ax = fig.add_subplot(subplot)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit=xtit,
                          ytit=ytit if ytitle else None, locators=(0.1, 1, 0.1))
        ax.text(xleg, yleg, 'z=%s' % str(z))

        #Predicted relation
        ind = np.where(mzr[idx,0,:] != 0)
        yplot = (mzr[idx,0,ind])
        errdn = (mzr[idx,1,ind])
        errup = (mzr[idx,2,ind])
        xplot = xmf[ind]
        ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'b', mec='b',marker='o',label="all galaxies")

        ind = np.where(mzr_cen[idx,0,:] != 0)
        yplot = (mzr_cen[idx,0,ind])
        errdn = (mzr_cen[idx,1,ind])
        errup = (mzr_cen[idx,2,ind])
        xplot = xmf[ind]
        ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'g', mec='g',marker='o',markersize=2,label="centrals")

        ind = np.where(mzr_sat[idx,0,:] != 0)
        yplot = (mzr_sat[idx,0,ind])
        errdn = (mzr_sat[idx,1,ind])
        errup = (mzr_sat[idx,2,ind])
        xplot = xmf[ind]

        ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'r', mec='r',marker='o',markersize=2,label="satellites")

        colors = ['b','g','r']
        if obs:
            # The observation has error bars or not
            if obs.yerrup is not None:
                common.errorbars(ax, obs.x, obs.y, obs.yerrdn, obs.yerrup, 'grey', 's', err_absolute=obs.err_absolute)
                colors.append('grey')
            else:
                ax.plot(obs.x, obs.y, 'k', label=obs.label)
                colors.insert(0, 'k')

        common.prepare_legend(ax, colors, loc=4)

    common.savefig(outdir, fig, 'mzr.pdf')



def plot_SFR_Mstars(plt, outdir, mainseqsf, mainseqsf_cen, mainseqsf_sat):


    fig = plt.figure(figsize=(9.5,9.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm star}/M_{\odot})$"
    ytit="$\\rm log_{10}(\\rm SFR/M_{\odot} yr^{-1})$"
    xmin, xmax, ymin, ymax = 7.1, 12, -3, 3
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    z0obs = observation('Brinchmann+04', [7,11], [-2,1], [-2.3,0.7], [-1.7,1.3], True)
    z1obs = observation('Elbaz+07', [9,10.7], [-0.3,1.5], [-0.6,1.2], [0,1.8], True)

    indeces = (0, 1, 2, 4)
    zs = (0, 0.5, 1, 2)
    subplots = (221, 222, 223, 224)
    observations = (z0obs, None, z1obs, None)

    for subplot, z, idx, obs in zip(subplots, zs, indeces, observations):

        ax = fig.add_subplot(subplot)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
        ax.text(xleg, yleg, 'z=%s' % str(z))

        #Predicted relation
        ind = np.where(mainseqsf[idx,0,:] != 0)
        xplot = xmf[ind]
        yplot = mainseqsf[idx,0,ind]
        errdn = mainseqsf[idx,1,ind]
        errup = mainseqsf[idx,2,ind]
        ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'b', mec='b',marker='o',label="all galaxies")

        ind = np.where(mainseqsf_cen[idx,0,:] != 0)
        xplot = xmf[ind]
        yplot = mainseqsf_cen[idx,0,ind]
        errdn = mainseqsf_cen[idx,1,ind]
        errup = mainseqsf_cen[idx,2,ind]
        ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'g', mec='g',marker='o',markersize=2,label="centrals")

        ind = np.where(mainseqsf_sat[idx,0,:] != 0)
        xplot = xmf[ind]
        yplot = mainseqsf_sat[idx,0,ind]
        errdn = mainseqsf_sat[idx,1,ind]
        errup = mainseqsf_sat[idx,2,ind]
        ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'r', mec='r',marker='o',markersize=2,label="satellites")

        if obs:
            ax.plot(obs.x, obs.y, 'k', label=obs.label)
            ax.plot(obs.x, obs.yerrup, 'k', linestyle='dotted')
            ax.plot(obs.x, obs.yerrdn, 'k', linestyle='dotted')

        colors = ['b', 'g', 'r']
        if obs:
            colors.insert(0, 'k')
        common.prepare_legend(ax, colors)

    common.savefig(outdir, fig, 'SFR_Mstars.pdf')


def plot_SFE_Mstars(plt, outdir, sfe, sfe_cen, sfe_sat):

    fig = plt.figure(figsize=(9.5,9.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm star}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\\tau_{H_2}/Gyr)$"
    xmin, xmax, ymin, ymax = 7.1, 12, -1, 1
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)


    indeces = (0, 1, 2, 4)
    zs = (0, 0.5, 1, 2)
    subplots = (221, 222, 223, 224)

    for subplot, z, idx, in zip(subplots, zs, indeces):

        ax = fig.add_subplot(subplot)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
        ax.text(xleg, yleg, 'z=%s' % str(z))

        #Predicted relation
        ind = np.where(sfe[idx,0,:] != 0)
        xplot = xmf[ind]
        yplot = sfe[idx,0,ind]
        errdn = sfe[idx,1,ind]
        errup = sfe[idx,2,ind]
        ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'b', mec='b',marker='o',label="all galaxies")

        ind = np.where(sfe_cen[idx,0,:] != 0)
        xplot = xmf[ind]
        yplot = sfe_cen[idx,0,ind]
        errdn = sfe_cen[idx,1,ind]
        errup = sfe_cen[idx,2,ind]
        ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'g', mec='g',marker='o',markersize=2,label="centrals")

        ind = np.where(sfe_sat[idx,0,:] != 0)
        xplot = xmf[ind]
        yplot = sfe_sat[idx,0,ind]
        errdn = sfe_sat[idx,1,ind]
        errup = sfe_sat[idx,2,ind]
        ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'r', mec='r',marker='o',markersize=2,label="satellites")

        common.prepare_legend(ax, ['b', 'g', 'r'])

    common.savefig(outdir, fig, 'SFE_Mstars.pdf')



def plot_fmzr(plt, outdir, fmzr):

    fig = plt.figure(figsize=(5,5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm star}/M_{\odot}) - 0.66 log_{10}(\\rm SFR/M_{\odot} yr^{-1})$"
    ytit = "$\\rm log_{10}(\\rm Z_{\\rm gas}/Z_{\odot})$"
    xmin, xmax, ymin, ymax = 7.1, 12, -3, 1
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg, yleg, 'z=0')

    #Predicted relation
    ind = np.where(fmzr[0,0,:] != 0)
    yplot = (fmzr[0,0,ind])
    errdn = (fmzr[0,1,ind])
    errup = (fmzr[0,2,ind])

    xplot = xmf[ind]
    ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'b', mec='b',marker='o',label="all galaxies")

    #observations approximate of Tremonti et al. (2004)
    xG = [7.3,10.5]
    yG = [-0.97,0.41]
    yGu = [-0.77,0.61]
    yGd = [-1.17,0.21]

    ax.plot(xG,yG,'k',label="Andrews & Martini (2013)")
    ax.plot(xG,yGu,'k',linestyle='dotted')
    ax.plot(xG,yGd,'k',linestyle='dotted')

    common.prepare_legend(ax, ['k','b'], loc=4)
    common.savefig(outdir, fig, 'fmzr.pdf')


def plot_mzr_z0(plt, outdir, obsdir, mzr_cen, mzr_sat):

    fig = plt.figure(figsize=(5,5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm star}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\\rm Z_{\\rm gas}/Z_{\odot})$"
    xmin, xmax, ymin, ymax = 8, 12, -3, 1

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #MZR z=0
    corrzsun = 8.69
    lm, mz, mzdn, mzup = common.load_observation(obsdir, 'MMAdrews13.dat', [0,1,2,3])
    common.errorbars(ax, lm, mz - corrzsun, mzdn - corrzsun, mzup - corrzsun, 'grey', 's', label='Andrews+13')
    lm, mz, mzdn, mzup = common.load_observation(obsdir, 'Tremonti04.dat', [0,1,2,3])
    common.errorbars(ax, lm, mz - corrzsun, mzdn - corrzsun, mzup - corrzsun, 'grey', 'o', label="Tremonti+04")

    ind = np.where(mzr_cen[0,0,:] != 0)
    yplot = (mzr_cen[0,0,ind])
    errdn = (mzr_cen[0,1,ind])
    errup = (mzr_cen[0,2,ind])
    xplot = xmf[ind]
    ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'g', mec='g',color='g',linestyle='solid',label="centrals")

    ind = np.where(mzr_sat[0,0,:] != 0)
    yplot = (mzr_sat[0,0,ind])
    errdn = (mzr_sat[0,1,ind])
    errup = (mzr_sat[0,2,ind])
    xplot = xmf[ind]
    ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'r', mec='r',color='r',linestyle='dashed',label="satellites")

    common.prepare_legend(ax, ['grey','grey','g','r'], loc=4)
    common.savefig(outdir, fig, 'mzr_z0.pdf')



def plot_sfr_mstars_z0(plt, outdir, obsdir, sfr_seq, mainseqsf):

    fig = plt.figure(figsize=(5,5))
    xtit="$\\rm log_{10} (\\rm M_{\\rm star}/M_{\odot})$"
    ytit="$\\rm log_{10}(\\rm SFR/M_{\odot} yr^{-1})$"

    xmin, xmax, ymin, ymax = 8, 12, -5, 3
    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #predicted relation
    ind = np.where((sfr_seq[0,:] > 0) & (sfr_seq[1,:] != 0) )
    xdata = sfr_seq[0,ind]
    ydata = sfr_seq[1,ind]
    us.density_contour(xdata[0], ydata[0], 30, 30, ax=ax) #, **contour_kwargs)

    ind = np.where(mainseqsf[0,0,:] != 0)
    xplot = xmf[ind]
    yplot = mainseqsf[0,0,ind]
    ax.plot(xplot,yplot[0],color='k',linestyle='solid', linewidth = 1, label="SHArk all galaxies")

    #SFR relation z=0
    lm, SFR = common.load_observation(obsdir, 'Brinchmann04.dat', (0, 1))
    ax.plot(lm[0:35], SFR[0:35], color='k', linewidth = 3, linestyle='dashed', label='Brinchmann+04')
    ax.plot(lm[36:70], SFR[36:70], color='k',linewidth = 5, linestyle='dotted')
    ax.plot(lm[71:len(SFR)], SFR[71:len(SFR)], color='k',linewidth = 5, linestyle='dotted')

    xdataD16 = [9.3, 10.6]
    ydataD16 = [-0.39, 0.477]
#     yrrupD16 = [-0.19, 0.677]
#     yrrdnD16 = [-0.59, 0.277]
    ax.plot(xdataD16,ydataD16, color='b',linestyle='dashdot',linewidth = 4, label='Davies+16')

    # Legend
    common.prepare_legend(ax, ['k','k', 'b'], loc=2)
    common.savefig(outdir, fig, 'SFR_Mstars_z0.pdf')


def prepare_data(hdf5_data, index, hist_smf, hist_smf_err, hist_smf_cen, hist_smf_sat,
                 hist_HImf, hist_HImf_cen, hist_HImf_sat, hist_H2mf, hist_H2mf_cen,
                 hist_H2mf_sat, mainseq, mainseqsf, sfe, mainseq_cen, mainseqsf_cen,
                 sfe_cen, mainseq_sat, mainseqsf_sat, sfe_sat, mzr, fmzr, mzr_cen,
                 mzr_sat, plotz, plotz_HImf):

    (h0, volh, sfr_disk, sfr_burst, mdisk, mbulge, rdisk_star, mBH, mHI, mH2,
     mgas_disk, mHI_bulge, mH2_bulge, mgas_bulge, mgas_metals_disk,
     mgas_metals_bulge, mstars_metals_disk, mstars_metals_bulge, typeg) = hdf5_data

    mgas = mgas_disk+mgas_bulge
    mgas_metals = mgas_metals_disk+mgas_metals_bulge

    mass = np.zeros(shape = len(mdisk))
    mass_atom = np.zeros(shape = len(mdisk))
    mass_atom_cen = np.zeros(shape = len(mdisk))
    mass_atom_sat = np.zeros(shape = len(mdisk))

    mass_mol = np.zeros(shape = len(mdisk))
    mass_mol_cen = np.zeros(shape = len(mdisk))
    mass_mol_sat = np.zeros(shape = len(mdisk))

    ind = np.where((mdisk+mbulge) > 0.0)
    mass[ind] = np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0))

    print('number of galaxies with mstars>0 and max mass: %d, %d' % (len(mass[ind]), max(mass[ind])))

    H, _ = np.histogram(mass,bins=np.append(mbins,mupp))
    hist_smf[index,:] = hist_smf[index,:] + H
    ran_err = np.random.normal(0.0, 0.25, len(mass))
    mass_err = mass + ran_err
    H, _ = np.histogram(mass_err,bins=np.append(mbins,mupp))
    hist_smf_err[index,:] = hist_smf_err[index,:] + H

    ind = np.where(typeg == 0)
    H, _ = np.histogram(mass[ind],bins=np.append(mbins,mupp))
    hist_smf_cen[index,:] = hist_smf_cen[index,:] + H
    ind = np.where(typeg > 0)
    H, _ = np.histogram(mass[ind],bins=np.append(mbins,mupp))
    hist_smf_sat[index,:] = hist_smf_sat[index,:] + H

    ind = np.where((mHI+mHI_bulge) > 0)
    mass_atom[ind] = np.log10(mHI[ind]+mHI_bulge[ind]) - np.log10(float(h0)) + np.log10(XH)

    H_HI, _ = np.histogram(mass_atom,bins=np.append(mbins,mupp))
    hist_HImf[index,:] = hist_HImf[index,:] + H_HI

    ind = np.where(((mHI+mHI_bulge) > 0) & (typeg == 0))
    mass_atom_cen[ind] = np.log10(mHI[ind]+mHI_bulge[ind]) - np.log10(float(h0)) + np.log10(XH)
    H_HI, _ = np.histogram(mass_atom_cen,bins=np.append(mbins,mupp))
    hist_HImf_cen[index,:] = hist_HImf_cen[index,:] + H_HI

    ind = np.where(((mHI+mHI_bulge) > 0) & (typeg > 0))
    mass_atom_sat[ind] = np.log10(mHI[ind]+mHI_bulge[ind]) - np.log10(float(h0)) + np.log10(XH)
    H_HI, _ = np.histogram(mass_atom_sat,bins=np.append(mbins,mupp))
    hist_HImf_sat[index,:] = hist_HImf_sat[index,:] + H_HI

    ind = np.where((mH2+mH2_bulge) > 0)
    mass_mol[ind] = np.log10(mH2[ind]+mH2_bulge[ind]) - np.log10(float(h0)) + np.log10(XH)
    H_H2, _ = np.histogram(mass_mol,bins=np.append(mbins,mupp))
    hist_H2mf[index,:] = hist_H2mf[index,:] + H_H2

    ind = np.where(((mH2+mH2_bulge) > 0) & (typeg == 0))
    mass_mol_cen[ind] = np.log10(mH2[ind]+mH2_bulge[ind]) - np.log10(float(h0)) + np.log10(XH)
    H_H2, _ = np.histogram(mass_mol_cen,bins=np.append(mbins,mupp))
    hist_H2mf_cen[index,:] = hist_H2mf_cen[index,:] + H_H2

    ind = np.where(((mH2+mH2_bulge) > 0) & (typeg > 0))
    mass_mol_sat[ind] = np.log10(mH2[ind]+mH2_bulge[ind]) - np.log10(float(h0)) + np.log10(XH)
    H_H2, _ = np.histogram(mass_mol_sat,bins=np.append(mbins,mupp))
    hist_H2mf_sat[index,:] = hist_H2mf_sat[index,:] + H_H2

    bin_it = functools.partial(us.wmedians, xbins=xmf)
    ind = np.where((sfr_disk+sfr_burst > 0) & (mdisk+mbulge > 0))
    mainseq[index,:] = bin_it(x=mass[ind], y=np.log10((sfr_disk[ind]+sfr_burst[ind])/(mdisk[ind]+mbulge[ind])))
    mainseqsf[index,:] = bin_it(x=mass[ind], y=np.log10((sfr_disk[ind]+sfr_burst[ind])/h0/GyrToYr))

    ind = np.where((sfr_disk+sfr_burst > 0) & (mH2+mH2_bulge > 0) & (mass > 0))
    sfe[index,:] = bin_it(x=mass[ind], y=np.log10((mH2[ind]+mH2_bulge[ind])/(sfr_disk[ind]+sfr_burst[ind])))

    ind = np.where((sfr_disk+sfr_burst > 0) & (typeg == 0) & (mdisk+mbulge > 0))
    mainseq_cen[index,:] = bin_it(x=mass[ind], y=np.log10((sfr_disk[ind]+sfr_burst[ind])/(mdisk[ind]+mbulge[ind])))
    mainseqsf_cen[index,:] = bin_it(x=mass[ind], y=np.log10((sfr_disk[ind]+sfr_burst[ind])/h0/GyrToYr))
    sfe_cen[index,:] = bin_it(x=mass[ind], y=np.log10((mH2[ind]+mH2_bulge[ind])/(sfr_disk[ind]+sfr_burst[ind])))

    ind = np.where((sfr_disk+sfr_burst > 0) & (typeg > 0) & (mdisk+mbulge > 0))
    mainseq_sat[index,:] = bin_it(x=mass[ind], y=np.log10((sfr_disk[ind]+sfr_burst[ind])/(mdisk[ind]+mbulge[ind])))
    mainseqsf_sat[index,:] = bin_it(x=mass[ind], y=np.log10((sfr_disk[ind]+sfr_burst[ind])/h0/GyrToYr))
    sfe_sat[index,:] = bin_it(x=mass[ind], y=np.log10((mH2[ind]+mH2_bulge[ind])/(sfr_disk[ind]+sfr_burst[ind])))

    ind = np.where((mgas_metals > 0.0) & (mgas > 0))
    mzr[index,:] = bin_it(x=mass[ind], y=np.log10((mgas_metals[ind]/mgas[ind]/Zsun)))

    ind = np.where((mgas_metals > 0.0) & (mgas > 0) & (sfr_disk+sfr_burst > 0))
    fmzr[index,:] = bin_it(x=mass[ind]-0.66*np.log10((sfr_disk[ind]+sfr_burst[ind])/h0/GyrToYr),
                           y=np.log10((mgas_metals[ind]/mgas[ind]/Zsun)))

    ind = np.where((mgas_metals > 0.0) & (typeg == 0) & (mgas > 0))
    mzr_cen[index,:] = bin_it(x=mass[ind], y=np.log10((mgas_metals[ind]/mgas[ind]/Zsun)))

    ind = np.where((mgas_metals > 0.0) & (typeg > 0) & (mgas > 0))
    mzr_sat[index,:] = bin_it(x=mass[ind], y=np.log10((mgas_metals[ind]/mgas[ind]/Zsun)))

    if volh > 0:
        vol = volh/pow(h0,3.)  # In Mpc^3
        hist_smf[index,:]  = hist_smf[index,:]/vol/dm
        hist_smf_err[index,:]  = hist_smf_err[index,:]/vol/dm
        hist_smf_cen[index,:]  = hist_smf_cen[index,:]/vol/dm
        hist_smf_sat[index,:]  = hist_smf_sat[index,:]/vol/dm

        hist_HImf[index,:] = hist_HImf[index,:]/vol/dm
        hist_HImf_cen[index,:] = hist_HImf_cen[index,:]/vol/dm
        hist_HImf_sat[index,:] = hist_HImf_sat[index,:]/vol/dm

        hist_H2mf[index,:] = hist_H2mf[index,:]/vol/dm
        hist_H2mf_cen[index,:] = hist_H2mf_cen[index,:]/vol/dm
        hist_H2mf_sat[index,:] = hist_H2mf_sat[index,:]/vol/dm

        plotz[index]     = True
        plotz_HImf[index]= True
    else:
        plotz[index]     = False
        plotz_HImf[index]= False

    return mass

def main():

    zlist = ["199","174", "156", "131", "113", "99"]

    plt = common.load_matplotlib()

    mainseq     = np.zeros(shape = (len(zlist), 3, len(xmf)))
    mainseq_cen = np.zeros(shape = (len(zlist), 3, len(xmf)))
    mainseq_sat = np.zeros(shape = (len(zlist), 3, len(xmf)))

    mainseqsf     = np.zeros(shape = (len(zlist), 3, len(xmf)))
    mainseqsf_cen = np.zeros(shape = (len(zlist), 3, len(xmf)))
    mainseqsf_sat = np.zeros(shape = (len(zlist), 3, len(xmf)))

    mzr         = np.zeros(shape = (len(zlist), 3, len(xmf)))
    mzr_cen     = np.zeros(shape = (len(zlist), 3, len(xmf)))
    mzr_sat     = np.zeros(shape = (len(zlist), 3, len(xmf)))
    fmzr        = np.zeros(shape = (len(zlist), 3, len(xmf)))

    sfe         = np.zeros(shape = (len(zlist), 3, len(xmf)))
    sfe_cen     = np.zeros(shape = (len(zlist), 3, len(xmf)))
    sfe_sat     = np.zeros(shape = (len(zlist), 3, len(xmf)))

    # Histograms
    hist_smf = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_err = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_cen = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_sat = np.zeros(shape = (len(zlist), len(mbins)))

    plotz = np.empty(shape=(len(zlist)), dtype=np.bool_)
    hist_HImf = np.zeros(shape = (len(zlist), len(mbins)))
    hist_HImf_cen = np.zeros(shape = (len(zlist), len(mbins)))
    hist_HImf_sat = np.zeros(shape = (len(zlist), len(mbins)))
    plotz_HImf = np.empty(shape=(len(zlist)), dtype=np.bool_)
    hist_H2mf = np.zeros(shape = (len(zlist), len(mbins)))
    hist_H2mf_cen = np.zeros(shape = (len(zlist), len(mbins)))
    hist_H2mf_sat = np.zeros(shape = (len(zlist), len(mbins)))

    fields = {'Galaxies': ('sfr_disk', 'sfr_burst', 'mstars_disk', 'mstars_bulge',
                           'rdisk_star', 'mBH', 'matom_disk', 'mmol_disk', 'mgas_disk',
                           'matom_bulge', 'mmol_bulge', 'mgas_bulge',
                           'mgas_metals_disk', 'mgas_metals_bulge',
                           'mstars_metals_disk', 'mstars_metals_bulge', 'type')}

    modeldir, subvolumes, outdir, obsdir = common.parse_args(requires_snapshot=False)
    for index in range(0,len(zlist)):
        hdf5_data = common.read_data(modeldir, zlist[index], fields, subvolumes)
        mass = prepare_data(hdf5_data, index, hist_smf, hist_smf_err, hist_smf_cen,
                             hist_smf_sat, hist_HImf, hist_HImf_cen, hist_HImf_sat,
                             hist_H2mf, hist_H2mf_cen, hist_H2mf_sat, mainseq, mainseqsf,
                             sfe, mainseq_cen, mainseqsf_cen, sfe_cen, mainseq_sat,
                             mainseqsf_sat, sfe_sat, mzr, fmzr, mzr_cen, mzr_sat, plotz,
                             plotz_HImf)

        h0 = hdf5_data[0]
        if index == 0:
            (sfr_disk, sfr_burst, mdisk, mbulge) = hdf5_data[2:6]
            sfr_seq = np.zeros(shape = (2, len(mdisk)))
            ind  = np.where((sfr_disk + sfr_burst > 0) & (mdisk + mbulge > 0))
            sfr_seq[0,ind] = mass[ind]
            sfr_seq[1,ind] = np.log10((sfr_disk[ind] + sfr_burst[ind]) / h0 / GyrToYr)


    # This should be the same in all HDF5 files

    # Take logs
    ind = np.where(hist_smf > 0.)
    hist_smf[ind] = np.log10(hist_smf[ind])
    ind = np.where(hist_smf_cen > 0.)
    hist_smf_cen[ind] = np.log10(hist_smf_cen[ind])
    ind = np.where(hist_smf_sat > 0.)
    hist_smf_sat[ind] = np.log10(hist_smf_sat[ind])
    ind = np.where(hist_smf_err > 0.)
    hist_smf_err[ind] = np.log10(hist_smf_err[ind])

    ind = np.where(hist_HImf > 0.)
    hist_HImf[ind] = np.log10(hist_HImf[ind])
    ind = np.where(hist_HImf_cen > 0.)
    hist_HImf_cen[ind] = np.log10(hist_HImf_cen[ind])
    ind = np.where(hist_HImf_sat > 0.)
    hist_HImf_sat[ind] = np.log10(hist_HImf_sat[ind])

    ind = np.where(hist_H2mf > 0.)
    hist_H2mf[ind] = np.log10(hist_H2mf[ind])
    ind = np.where(hist_H2mf_cen > 0.)
    hist_H2mf_cen[ind] = np.log10(hist_H2mf_cen[ind])
    ind = np.where(hist_H2mf_sat > 0.)
    hist_H2mf_sat[ind] = np.log10(hist_H2mf_sat[ind])

    plot_stellarmf_z(plt, outdir, obsdir, h0, plotz, hist_smf, hist_smf_cen, hist_smf_sat, hist_smf_err)
    plot_HImf_z0(plt, outdir, obsdir, h0, plotz_HImf, hist_HImf, hist_HImf_cen, hist_HImf_sat)
    plot_H2mf_z0(plt, outdir, obsdir, h0, plotz_HImf, hist_H2mf, hist_H2mf_cen, hist_H2mf_sat)
    plot_SSFR_Mstars(plt, outdir, mainseq, mainseq_cen, mainseq_sat)
    plot_mzr(plt, outdir, obsdir, mzr, mzr_cen, mzr_sat)
    plot_SFR_Mstars(plt, outdir, mainseqsf, mainseqsf_cen, mainseqsf_sat)
    plot_SFE_Mstars(plt, outdir, sfe, sfe_cen, sfe_sat)
    plot_fmzr(plt, outdir, fmzr)
    plot_mzr_z0(plt, outdir, obsdir, mzr_cen, mzr_sat)
    plot_sfr_mstars_z0(plt, outdir, obsdir, sfr_seq, mainseqsf)



if __name__ == '__main__':
    main()
