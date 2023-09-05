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
"""SMF plots"""

import collections
import functools
import logging
import math

import numpy as np

import common
import utilities_statistics as us


observation = collections.namedtuple('observation', 'label x y yerrup yerrdn err_absolute')

logger = logging.getLogger(__name__)

##################################
# Constants
GyrToYr = 1e9
Zsun = 0.0127
XH = 0.72
MpcToKpc = 1e3

##################################
# Mass function initialization
mlow = 5
mupp = 14
dm = 0.125
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0
imf   = 'cha'

mlow2 = 5
mupp2 = 14
dm2 = 0.3
mbins2 = np.arange(mlow2,mupp2,dm2)
xmf2 = mbins2 + dm2/2.0

mlow3 = 5
mupp3 = 14
dm3 = 0.2
mbins3 = np.arange(mlow3,mupp3,dm3)
xmf3 = mbins3 + dm3/2.0

ssfrlow = -6
ssfrupp = 4
dssfr = 0.2
ssfrbins = np.arange(ssfrlow,ssfrupp,dssfr)
xssfr    = ssfrbins + dssfr/2.0

sfrlow = -3
sfrupp = 1.5
dsfr = 0.2
sfrbins = np.arange(sfrlow,sfrupp,dsfr)
xsfr    = sfrbins + dsfr/2.0


def load_smf_observations(obsdir, h0):
    # Thorne et al. (2021)
    def add_thorne21_data(zappend, file_read='z0.4'):
        lm, pD, dn, du = np.loadtxt(obsdir+'/mf/SMF/Thorne21/SMFvals_'+file_read+'.csv', delimiter=',', skiprows=1, usecols = [0,1,2,3], unpack = True)
        hobd = 0.7
        pDlog = np.log10(pD[::3]) +  3.0 * np.log10(hobs/h0)
        dnlog = np.log10(pD[::3]) - np.log10(dn[::3])
        dulog = np.log10(du[::3]) - np.log10(pD[::3])
        lm = lm[::3] -  2.0 * np.log10(hobs/h0)
        zappend.append((observation("Thorne+2021", lm, pDlog, dnlog, dulog, err_absolute=False), 's'))

    # Weaver et al. (2022; COSMOS2020)
    def add_weaver22_data(zappend, file_read='0.2z0.5'):
        lm, pD, dn, du = np.loadtxt(obsdir+'/mf/SMF/COSMOS2020/SMF_Farmer_v2.1_' + file_read + '_total.txt', delimiter=' ', skiprows=0, usecols = [0,2,3,4], unpack = True)
        hobd = 0.7
        pDlog = np.log10(pD) +  3.0 * np.log10(hobs/h0)
        dnlog = np.log10(pD) - np.log10(dn)
        dulog = np.log10(du) - np.log10(pD)
        lm = lm -  2.0 * np.log10(hobs/h0)
        zappend.append((observation("Weaver+2023", lm, pDlog, dnlog, dulog, err_absolute=False), '*'))


    # Driver al. (2022, z=0). Chabrier IMF
    z0obs = []
    lm, p, dp = common.load_observation(obsdir, 'mf/SMF/GAMAIV_Driver22.dat', [0,1,2])
    hobs = 0.7
    xobs = lm + 2.0 * np.log10(hobs/h0)
    yobs = p - 3.0 * np.log10(hobs/h0)
    z0obs.append((observation("Driver+2022", xobs, yobs, dp, dp, err_absolute=False), 'o'))

    lm, p, dpdn, dpup = common.load_observation(obsdir, 'mf/SMF/SMF_Bernardi2013_SerExp.data', [0,1,2,3])
    xobs = lm + 2.0 * np.log10(hobs/h0)
    yobs = np.log10(p) - 3.0 * np.log10(hobs/h0)
    ydn = np.log10(p) - np.log10(p-dpdn)
    yup = np.log10(p+dpup) - np.log10(p) 
    z0obs.append((observation("Bernardi+2013", xobs, yobs, ydn, yup, err_absolute=False), 's'))


    lm, p, dpdn, dpup = common.load_observation(obsdir, 'mf/SMF/SMF_Li2009.dat', [0,1,2,3])
    xobs = lm - 2.0 * np.log10(hobs) + 2.0 * np.log10(hobs/h0)
    yobs = p + 3.0 * np.log10(hobs) - 3.0 * np.log10(hobs/h0)
    z0obs.append((observation("Li&White+2009", xobs, yobs, abs(dpdn), dpup, err_absolute=False), 'd'))


    # Moustakas (Chabrier IMF), ['Moustakas+2013, several redshifts']
    zdnM13, lmM13, pM13, dp_dn_M13, dp_up_M13 = common.load_observation(obsdir, 'mf/SMF/SMF_Moustakas2013.dat', [0,3,5,6,7])
    xobsM13 = lmM13 + 2.0 * np.log10(hobs/h0)

    yobsM13 = np.full(xobsM13.shape, -999.) - 3.0 * np.log10(hobs/h0)
    lerrM13 = np.full(xobsM13.shape, -999.)
    herrM13 = np.full(xobsM13.shape, 999.)
    indx = np.where( pM13 < 1)
    yobsM13[indx] = (pM13[indx])
    indx = np.where( dp_dn_M13 > 0)
    lerrM13[indx]  = dp_dn_M13[indx] 
    indx = np.where( dp_up_M13 > 0)
    herrM13[indx]  = dp_up_M13[indx]

    # Muzzin (Kroupa IMF), ['Moustakas+2013, several redshifts']
    zdnMu13,zupMu13,lmMu13,pMu13,dp_dn_Mu13,dp_up_Mu13 = common.load_observation(obsdir, 'mf/SMF/SMF_Muzzin2013.dat', [0,1,2,4,5,5])
    # -0.09 corresponds to the IMF correction
    xobsMu13 = lmMu13 - 0.09 + 2.0 * np.log10(hobs/h0) 
    yobsMu13 = np.full(xobsMu13.shape, -999.) - 3.0 * np.log10(hobs/h0)
    lerrMu13 = np.full(xobsMu13.shape, -999.)
    herrMu13 = np.full(xobsMu13.shape, 999.)
    indx = np.where( pMu13 < 1)
    yobsMu13[indx] = (pMu13[indx])
    indx = np.where( dp_dn_Mu13 > 0)
    lerrMu13[indx]  = dp_dn_Mu13[indx] 
    indx = np.where( dp_up_Mu13 > 0)
    herrMu13[indx]  = dp_up_Mu13[indx]

    # z0.5 obs
    z05obs = []
    #in_redshift = np.where(zdnM13 == 0.4)
    #z05obs.append((observation("Moustakas+2013", xobsM13[in_redshift], yobsM13[in_redshift], lerrM13[in_redshift], herrM13[in_redshift], err_absolute=False), 'o'))
    in_redshift = np.where(zdnMu13 == 0.5)
    z05obs.append((observation("Muzzin+2013", xobsMu13[in_redshift], yobsMu13[in_redshift], lerrMu13[in_redshift], herrMu13[in_redshift], err_absolute=False), 'o'))
    add_thorne21_data(z05obs, file_read='z0.51')
    add_weaver22_data(z05obs, file_read='0.2z0.5')


    # z1 obs
    z1obs = []
    #in_redshift = np.where(zdnM13 == 0.8)
    #z1obs.append((observation("Moustakas+2013", xobsM13[in_redshift], yobsM13[in_redshift], lerrM13[in_redshift], herrM13[in_redshift], err_absolute=False), 'o'))
    in_redshift = np.where(zdnMu13 == 1)
    z1obs.append((observation("Muzzin+2013", xobsMu13[in_redshift], yobsMu13[in_redshift], lerrMu13[in_redshift], herrMu13[in_redshift], err_absolute=False), 'o'))
    add_thorne21_data(z1obs, file_read='z1.1')
    add_weaver22_data(z1obs, file_read='0.8z1.1')

    #z2 obs
    z2obs = []
    in_redshift = np.where(zupMu13 == 2.5)
    z2obs.append((observation("Muzzin+2013", xobsMu13[in_redshift], yobsMu13[in_redshift], lerrMu13[in_redshift], herrMu13[in_redshift], err_absolute=False), 'o'))
    #in_redshift = np.where(zdnS12 == 1.8)
    #z2obs.append((observation("Santini+2012", xobsS12[in_redshift], yobsS12[in_redshift], lerrS12[in_redshift], herrS12[in_redshift], err_absolute=False), 'o'))
    add_thorne21_data(z2obs, file_read='z2')
    add_weaver22_data(z2obs, file_read='2.0z2.5')

    # z3 obs
    z3obs = []
    in_redshift = np.where(zupMu13 == 3.0)
    z3obs.append((observation("Muzzin+2013", xobsMu13[in_redshift], yobsMu13[in_redshift], lerrMu13[in_redshift], herrMu13[in_redshift], err_absolute=False), 'o'))
    #in_redshift = np.where(zdnS12 == 2.5)
    #z3obs.append((observation("Santini+2012", xobsS12[in_redshift], yobsS12[in_redshift], lerrS12[in_redshift], herrS12[in_redshift], err_absolute=False), 'o'))
    add_thorne21_data(z3obs, file_read='z3')
    add_weaver22_data(z3obs, file_read='3.0z3.5')

    # z4 obs
    z4obs = []
    in_redshift = np.where(zupMu13 == 4.0)
    z4obs.append((observation("Muzzin+2013", xobsMu13[in_redshift], yobsMu13[in_redshift], lerrMu13[in_redshift], herrMu13[in_redshift], err_absolute=False), 'o'))
    #in_redshift = np.where(zdnS12 == 3.5)
    #z4obs.append((observation("Santini+2012", xobsS12[in_redshift], yobsS12[in_redshift], lerrS12[in_redshift], herrS12[in_redshift], err_absolute=False), 'o'))
    add_thorne21_data(z4obs, file_read='z4')
    add_weaver22_data(z4obs, file_read='3.5z4.5')

    return (z0obs, z05obs, z1obs, z2obs, z3obs, z4obs)

def plot_stellarmf_z(plt, outdir, obsdir, h0, plotz, hist_smf, hist_smf_cen, hist_smf_sat, hist_smf_offset, hist_smf_30kpc):

    (z0obs, z05obs, z1obs, z2obs, z3obs, z4obs) = load_smf_observations(obsdir, h0)
    
    PlotLagos18 = True
    def plot_lagos18_smf(ax, z):
        sm, z0, z0p5, z1, z2, z3, z4 = common.load_observation(obsdir, 'Models/SharkVariations/SMF_Lagos18.dat', [0,1,2,3,4,5,6])
        if z == 0:
           y = z0
        elif z == 1:
           y = z0p5
        elif z == 2:
           y = z1
        elif z == 3:
           y = z2
        elif z == 4:
           y = z3
        elif z == 5:
           y = z4
        ax.plot(sm, y, linestyle='dashed', color='black',label='Shark v1.1 (L18)' if z == 0 else None)

    
    fig = plt.figure(figsize=(9.7,11.7))
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\Phi/dlog_{10}{\\rm M_{\\star}}/{\\rm Mpc}^{-3} )$"
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
                             marker, err_absolute=obs.err_absolute, label=obs.label, markersize=4)

        # Predicted SMF
        if plotz[idx]:
            y = hist_smf[idx,:]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],'r', label='Shark v2.0' if idx == 0 else None)
            #print('#stellar mass function redshift', z)
            #for a,b in zip(xmf[ind],y[ind]):
            #    print(a,b)
            #y = hist_smf_cen[idx,:]
            #ind = np.where(y < 0.)
            #ax.plot(xmf[ind],y[ind],'b', linestyle='dotted', label ='centrals' if idx == 0 else None)
            #y = hist_smf_sat[idx,:]
            #ind = np.where(y < 0.)
            #ax.plot(xmf[ind],y[ind],'g', linestyle='dotted', label ='Shark v2.0 (sats)' if idx == 0 else None)
            if z < 1:
                y = hist_smf_30kpc[idx,:]
                ind = np.where(y < 0.)
                ax.plot(xmf[ind],y[ind],'r', linestyle='dotted', linewidth=1, label ='Shark v2.0 (30kpc)'  if idx == 0 else None)
            if z >= 1:
                y = hist_smf_offset[idx,:]
                ind = np.where(y < 0.)
                ax.plot(xmf[ind],y[ind],'r', linestyle='dashdot', linewidth=2, label ='0.3dex error')
            if PlotLagos18 == True:
               plot_lagos18_smf(ax, idx)
       

        colors = []
        if idx == 0:
            colors = ['r','r','k']
        #elif idx == 1:
        #    colors += ['r']
        elif idx > 1:
            colors = ['r']
        colors += ['grey', 'grey','grey']

        common.prepare_legend(ax, colors)

    common.savefig(outdir, fig, 'stellarmf_z.pdf')


def plot_stellarmf_galcomponents(plt, outdir, obsdir, h0, plotz, hist_smf, hist_smf_comp):

    (z0obs, z05obs, z1obs, z2obs, z3obs, z4obs) = load_smf_observations(obsdir, h0)

    #Model comparison plot
    fig = plt.figure(figsize=(4.5,8))
    ytit = "$\\rm log_{10}(\Phi/dlog_{10}{\\rm M_{\\star}}/{\\rm Mpc}^{-3} )$"
    xmin, xmax, ymin, ymax = 8, 12.8, -5, -1
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (211, 212)
    indeces = (0, 2)
    zs = (0, 1)
    observations = (z0obs, z1obs)

    for subplot, idx, z, obs_and_markers in zip(subplots, indeces, zs, observations):

        ax = fig.add_subplot(subplot)
        if(idx != 3):
            xtit = ""
        else:
            xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
        plt.subplots_adjust(left=0.2)
        ax.text(xleg, yleg, 'z=%s' % str(z))

        # Observations
        for obs, marker in obs_and_markers:
            common.errorbars(ax, obs.x, obs.y, obs.yerrdn, obs.yerrup, 'grey',
                             marker, err_absolute=obs.err_absolute, label=obs.label)

        # Predicted SMF
        if plotz[idx]:
            y = hist_smf[idx,:]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],'k', label='Shark' if idx == 0 else None)
            y = hist_smf_comp[idx,0,:]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],'b', label='disks' if idx == 0 else None)
            y = hist_smf_comp[idx,1,:]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],'r', label='bulges' if idx == 0 else None)

            y = hist_smf_comp[idx,2,:]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],'PaleVioletRed', label='bulges by mergers' if idx == 0 else None)
            y = hist_smf_comp[idx,3,:]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],'cyan', label='bulges by diskins' if idx == 0 else None)

        colors = []
        if idx == 0:
            colors = ['k','b','r','PaleVioletRed','cyan']
            common.prepare_legend(ax, colors)

    common.savefig(outdir, fig, 'stellarmf_z_bycomponent.pdf')

def plot_stellarmf_z_molcomp(plt, outdir, obsdir, h0, plotz, hist_smf):

    hist_smf_modelvar = np.zeros(shape = (6, 360))
    hist_smf_modelvar[0,:], hist_smf_modelvar[1,:],hist_smf_modelvar[2,:],hist_smf_modelvar[3,:],hist_smf_modelvar[4,:],hist_smf_modelvar[5,:] = common.load_observation(obsdir, 'Models/SharkVariations/SMF_OtherModels.dat', [0,1,2,3,4,5])

    hist_smf_resolution = np.zeros(shape = (6, 135))
    hist_smf_resolution[0,:], hist_smf_resolution[1,:],hist_smf_resolution[2,:],hist_smf_resolution[3,:],hist_smf_resolution[4,:],hist_smf_resolution[5,:] = common.load_observation(obsdir, 'Models/SharkVariations/SMF_Resolution.dat', [0,1,2,3,4,5])

    #Model comparison plot
    fig = plt.figure(figsize=(4.5,8))
    ytit = "$\\rm log_{10}(\Phi/dlog_{10}{\\rm M_{\\star}}/{\\rm Mpc}^{-3} )$"
    xmin, xmax, ymin, ymax = 8, 12.8, -5, -1
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (211, 212)
    indeces = (0, 3)
    zs = (0, 2)
    for subplot, idx, z in zip(subplots, indeces, zs):

        ax = fig.add_subplot(subplot)
        if(idx != 3):
            xtit = ""
        else:
            xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
        plt.subplots_adjust(left=0.2)
        ax.text(xleg, yleg, 'z=%s' % str(z))

        # Predicted SMF
        if plotz[idx]:
            y = hist_smf[idx,:]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],'k', label='Shark' if idx == 0 else None)
            y = hist_smf_modelvar[idx,0:44]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],'r', linestyle='dotted', label='$z_{\\rm P} = 0$ (SN)' if idx == 0 else None)
            y = hist_smf_modelvar[idx,135:179]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],color='Crimson', linestyle='dashed', label='$\\beta = 2$ (SN)' if idx == 0 else None)
            y = hist_smf_modelvar[idx,90:134]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],'b', linestyle='dashdot', label='$f_{\\rm df} = 0$' if idx == 0 else None)
            y = hist_smf_modelvar[idx,180:225]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],'g', linestyle='dashed', label='$\\tau_{\\rm reinc} = 0$' if idx == 0 else None)
            y = hist_smf_modelvar[idx,270:314]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],color= 'BurlyWood', linestyle='dashed', label='$\\kappa_{\\rm r} = 0.0002$' if idx == 0 else None)

        colors = []
        if idx == 0:
            colors = ['k','r','Crimson','b','g','BurlyWood']
            common.prepare_legend(ax, colors)

    common.savefig(outdir, fig, 'stellarmf_z_modelcomparison.pdf')

    #Resolution test plot
    fig = plt.figure(figsize=(4.5,8))
    subplots = (211, 212)
    indeces = (0, 3)
    zs = (0, 2)
    for subplot, idx, z in zip(subplots, indeces, zs):

        ax = fig.add_subplot(subplot)
        if(idx != 3):
            xtit = ""
        else:
            xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
        plt.subplots_adjust(left=0.2)
        ax.text(xleg, yleg, 'z=%s' % str(z))

        # Predicted SMF
        if plotz[idx]:
            y = hist_smf[idx,:]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],'k', label='L210N1536' if idx == 0 else None)
            y = hist_smf_resolution[idx,0:44]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],'r', linestyle='dotted', label='L210N512' if idx == 0 else None)
            y = hist_smf_resolution[idx,45:89]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],color='g', linestyle='dashed', label='L40N512' if idx == 0 else None)
            y = hist_smf_resolution[idx,90:134]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],'b', linestyle='dashdot', label='L210N1024' if idx == 0 else None)

        colors = []
        if idx == 0:
            colors = ['k','r','g','b']
            common.prepare_legend(ax, colors)

    common.savefig(outdir, fig, 'stellarmf_z_resolutioncomparison.pdf')

def plot_HImf_z0(plt, outdir, obsdir, h0, plotz_HImf, hist_HImf, hist_HImf_cen, hist_HImf_sat, hist_HImf2):

    fig = plt.figure(figsize=(5,4.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm HI}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\Phi/dlog_{10}{\\rm M_{\\rm HI}}/{\\rm Mpc}^{-3} )$"
    xmin, xmax, ymin, ymax = 7.1, 12, -6, 0
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #HIPASS
    lmHI, pHI, dpHIdn, dpHIup = common.load_observation(obsdir, 'mf/GasMF/HIMF_Zwaan2005.dat', [0,1,2,3])

    #correct data for their choice of cosmology
    hobs = 0.75
    xobs = lmHI + np.log10(pow(hobs,2)/pow(h0,2))
    yobs = pHI + np.log10(pow(h0,3)/pow(hobs,3))
    ax.errorbar(xobs, yobs, yerr=[dpHIdn,dpHIup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Zwaan+2005")

    #ALFALFA.40
    lmHI, pHI, pdnHI, pduHI = common.load_observation(obsdir, 'mf/GasMF/HIMF_Jones18.dat', [0,1,2,3])

    #correct data for their choice of cosmology
    dpdnHI = pHI - pdnHI
    dpduHI = pduHI - pHI
    hobs = 0.7
    xobs = lmHI + np.log10(pow(hobs,2)/pow(h0,2))
    yobs = pHI + np.log10(pow(h0,3)/pow(hobs,3))

    ax.errorbar(xobs, yobs, yerr=[dpdnHI,dpduHI], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='x',label="Jones+2018")


    # Predicted HIMF
    if plotz_HImf[0]:
        #y = hist_HImf2[0,:]
        #ind = np.where(y < 0.)
        #ax.plot(xmf[ind],y[ind],'r', linestyle='solid',  label ='Shark v2.0')
        y = hist_HImf[0,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'r', linestyle='solid',  label ='Shark v2.0')

    xm, ym = common.load_observation(obsdir, 'Models/SharkVariations/HIMF_Lagos18.dat', [0,1])
    ax.plot(xm,ym,'k', linestyle='dashed', label = 'Shark v1.1 (L18)')

    common.prepare_legend(ax, ['r','k','grey','grey'])
    common.savefig(outdir, fig, 'HImf_z0.pdf')

    # Plot different resolutions
    fig = plt.figure(figsize=(5,4.5))

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)
    xmin = 6
    ymax = -0.5
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    # Predicted HIMF
    if plotz_HImf[0]:
        y = hist_HImf[0,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'k',  label ='L210N1536')

    hist_hi_resolution = common.load_observation(obsdir, 'Models/SharkVariations/HIH2MassFunctions_Resolution.dat', [0])

    y = hist_hi_resolution[0:44]
    ind = np.where(y < 0.)
    ax.plot(xmf[ind],y[ind],'r', linestyle='dotted', label='L210N512')
    y = hist_hi_resolution[45:89]
    ind = np.where(y < 0.)
    ax.plot(xmf[ind],y[ind],color='g', linestyle='dashed', label='L40N512')
    y = hist_hi_resolution[90:134]
    ind = np.where(y < 0.)
    ax.plot(xmf[ind],y[ind],'b', linestyle='dashdot', label='L210N1024')
    y = hist_hi_resolution[135:178]
    ind = np.where(y < 0.)
    ax.plot(xmf[ind],y[ind],'MediumAquamarine', linestyle='dotted', label='L40N512, $v_{\\rm cut}=30\\rm km\\ s^{-1}$')
    ax.text(8.4, -5.9, '$\\alpha_{\\rm V}=-0.85$', fontsize=12, color='MediumAquamarine')

    common.prepare_legend(ax, ['k','r','g','b','MediumAquamarine'])

    common.savefig(outdir, fig, 'HImf_z0_resolution.pdf')

    fig = plt.figure(figsize=(5,4.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm HI}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\Phi/dlog_{10}{\\rm M_{\\rm HI}}/{\\rm Mpc}^{-3} )$"
    xmin, xmax, ymin, ymax = 7.1, 12, -6, 0
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    labels=('z=0','z=0.5','z=1','z=2','z=3','z=4')
    cols=('red','LightSalmon','LimeGreen','DarkGreen','DarkTurquoise','blue')
    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #HIPASS
    lmHI, pHI, dpHIdn, dpHIup = common.load_observation(obsdir, 'mf/GasMF/HIMF_Zwaan2005.dat', [0,1,2,3])

    #correct data for their choice of cosmology
    hobs = 0.75
    xobs = lmHI + np.log10(pow(hobs,2)/pow(h0,2))
    yobs = pHI + np.log10(pow(h0,3)/pow(hobs,3))
    ax.errorbar(xobs, yobs, yerr=[dpHIdn,dpHIup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Zwaan+2005")

    #ALFALFA.40
    lmHI, pHI, pdnHI, pduHI = common.load_observation(obsdir, 'mf/GasMF/HIMF_Jones18.dat', [0,1,2,3])

    #correct data for their choice of cosmology
    dpdnHI = pHI - pdnHI
    dpduHI = pduHI - pHI
    hobs = 0.7
    xobs = lmHI + np.log10(pow(hobs,2)/pow(h0,2))
    yobs = pHI + np.log10(pow(h0,3)/pow(hobs,3))

    ax.errorbar(xobs, yobs, yerr=[dpdnHI,dpduHI], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='x',label="Jones+2018")


    # Predicted HIMF
    for z in range(0,len(hist_HImf[:,0])):
        y = hist_HImf[z,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],color=cols[z],  label =labels[z])

    common.prepare_legend(ax, cols)
    common.savefig(outdir, fig, 'HImf_evo.pdf')

def plot_H2mf_z0(plt, outdir, obsdir, h0, plotz_HImf, hist_H2mf, hist_H2mf_cen, hist_H2mf_sat):

    fig = plt.figure(figsize=(5,4.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm H_2}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\Phi/dlog_{10}{\\rm M_{\\rm H_2}}/{\\rm Mpc}^{-3} )$"
    xmin, xmax, ymin, ymax = 7.1, 12, -6, 0
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    # H2 mass function
    lmCO, pCO, dpCOdn, dpCOup = common.load_observation(obsdir, 'mf/GasMF/Keres03_LCOLF60m.dat', [0,1,2,3])

    # correct data for their choice of cosmology
    hobs = 0.75
    dpCOdn = np.abs(dpCOdn-pCO)
    dpCOup = np.abs(dpCOup-pCO)
    xobs = lmCO + np.log10(pow(hobs,2)/pow(h0,2))
    yobs = pCO + np.log10(pow(h0,3)/pow(hobs,3))

    # apply a constant MW conversion factor.
    X = 2.0
    corr_fac_H2 = np.log10(580.*X)+2.*np.log10(2.6)-np.log10(4.*math.pi)
    ax.errorbar(xobs+corr_fac_H2, yobs, yerr=[dpCOdn,dpCOup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='^',label="Keres+03")

    #H2 mass function
    lm,p,dpdn,dpup = common.load_observation(obsdir, 'mf/GasMF/H2MF_Fletcher21_Estimated.dat', [0,1,2,3])
    #correct data for their choice of cosmology
    #add bin to the data.
    hobs = 0.7
    lm = lm + np.log10(pow(hobs,2)/pow(h0,2))
    yobs = p + np.log10(pow(h0,3)/pow(hobs,3))
    ax.errorbar(lm, yobs, yerr=[p-dpdn,dpup-p], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s',label="Fletcher+21")

    #Predicted H2MF
    if plotz_HImf[0]:
        y = hist_H2mf[0,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'r', linestyle='solid', label='Shark v2.0')
        print("Will print the H2 MF")
        for a,b in zip(xmf[ind],y[ind]):
            print(a,b)

        #y = hist_H2mf_cen[0,:]
        #ind = np.where(y < 0.)
        #ax.plot(xmf[ind],y[ind],'b', linestyle='dotted')
        #y = hist_H2mf_sat[0,:]
        #ind = np.where(y < 0.)
        #ax.plot(xmf[ind],y[ind],'r', linestyle='dashed')

    #pH2_GD14 = common.load_observation(obsdir, 'Models/SharkVariations/HIH2MassFunctions_OtherModels.dat', [1])

    #ind = np.where(pH2_GD14 < 0)
    #ax.plot(xmf[ind],pH2_GD14[ind],'BurlyWood', linestyle='dashdot') 

    xm, ym = common.load_observation(obsdir, 'Models/SharkVariations/HIMF_Lagos18.dat', [0,2])
    ax.plot(xm,ym,'k', linestyle='dashed')

    common.prepare_legend(ax, ['red','grey','grey'])
    common.savefig(outdir, fig, 'H2mf_z0.pdf')

def plot_SSFR_Mstars(plt, outdir, mainseq, mainseq_cen, mainseq_sat):

    fig = plt.figure(figsize=(9.5,9.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
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
        #ind = np.where(mainseq[idx,0,:] != 0)
        #xplot = xmf[ind]
        #yplot = mainseq[idx,0,ind]
        #errdn = mainseq[idx,1,ind]
        #errup = mainseq[idx,2,ind]

        #ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'b', mec='b',marker='o',label="all galaxies")

        ind = np.where(mainseq_cen[idx,0,:] != 0)
        xplot = xmf[ind]
        yplot = mainseq_cen[idx,0,ind]
        errdn = mainseq_cen[idx,1,ind]
        errup = mainseq_cen[idx,2,ind]

        ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'g', mec='g',marker='o',markersize=2,label="centrals")
        #ind = np.where(mainseq_sat[idx,0,:] != 0)
        #xplot = xmf[ind]
        #yplot = mainseq_sat[idx,0,ind]
        #errdn = mainseq_sat[idx,1,ind]
        #errup = mainseq_sat[idx,2,ind]

        #ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'r', mec='r',marker='o',markersize=2,label="satellites")

        # observations approximate of Gilbank et al. (2010)
        if obs:
            ax.plot(obs.x, obs.y, 'k', label=obs.label)

        colors = ['g']
        if obs:
            colors.insert(0, 'k')
        common.prepare_legend(ax, colors)

    common.savefig(outdir, fig, 'SSFR_Mstars.pdf')


def plot_mzr(plt, outdir, obsdir, h0, mzr, mzr_cen, mzr_sat):

    fig = plt.figure(figsize=(9.5,9.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\\rm Z_{\\rm gas}/Z_{\odot})$"
    xmin, xmax, ymin, ymax = 7.1, 12, -3, 1
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    # Observations at z=0
    lm, mz, mzdn, mzup = common.load_observation(obsdir, 'MZR/MMAdrews13.dat', [0,1,2,3])
    corrzsun = 8.69
    hobs = 0.7
    #add cosmology correction plus IMF correction that goes into the stellar mass.
    corr_cos = np.log10(pow(hobs,2)/pow(h0,2)) - 0.09
    z0obs = observation('Andrews & Martini (2013)', lm + corr_cos, mz - corrzsun, mzdn - corrzsun, mzup - corrzsun, True)

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
        common.errorbars(ax, xplot, yplot[0], errdn[0], errup[0], 'b', 'o', label="all galaxies", err_absolute=False)

        ind = np.where(mzr_cen[idx,0,:] != 0)
        yplot = (mzr_cen[idx,0,ind])
        errdn = (mzr_cen[idx,1,ind])
        errup = (mzr_cen[idx,2,ind])
        xplot = xmf[ind]
        common.errorbars(ax, xplot, yplot[0], errdn[0], errup[0], 'g', 'o', markersize=2, label="centrals", err_absolute=False)

        ind = np.where(mzr_sat[idx,0,:] != 0)
        yplot = (mzr_sat[idx,0,ind])
        errdn = (mzr_sat[idx,1,ind])
        errup = (mzr_sat[idx,2,ind])
        xplot = xmf[ind]
        common.errorbars(ax, xplot, yplot[0], errdn[0], errup[0], 'r', 'o', markersize=2, label="satellites", err_absolute=False)

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


def plot_SFR_Mstars(plt, outdir, obsdir, mainseqsf, mainseqsf_cen, mainseqsf_sat, mainseqsf_1s, mainseqHI, mainseqH2):


    fig = plt.figure(figsize=(4.5,9))
    plt.subplots_adjust(left=0.15,bottom=0.1)

    xtit=' '
    ytit="$\\rm log_{10}(\\rm SFR/M_{\odot} yr^{-1})$"
    xmin, xmax, ymin, ymax = 8., 12, -3, 2
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymin + 0.1 * (ymax - ymin)


    mainseqsf_modelvar = np.zeros(shape = (5, 270))
    mainseqsf_modelvar[0,:], mainseqsf_modelvar[1,:], mainseqsf_modelvar[2,:], mainseqsf_modelvar[3,:], mainseqsf_modelvar[4,:] = common.load_observation(obsdir, 'Models/SharkVariations/SFRMstars_OtherModels.dat', [0,1,2,3,4])

    mainseqsf_GD14  = np.zeros(shape = (2, 5, len(xmf3)))
    mainseqsf_KMT09 = np.zeros(shape = (2, 5, len(xmf3)))
    mainseqsf_K13   = np.zeros(shape = (2, 5, len(xmf3)))

    for z in range(0,2):
        st = len(xmf3) * z
        for j in range(0,5):
            for i in range(0,len(xmf3)):
                mainseqsf_GD14[z,j,i] = mainseqsf_modelvar[j,st+i]

    for z in range(0,2):
        st = len(xmf3) * (z+2)
        for j in range(0,5):
            for i in range(0,len(xmf3)):
                mainseqsf_KMT09[z,j,i] = mainseqsf_modelvar[j,st+i]

    for z in range(0,2):
        st = len(xmf3) * (z+4)
        for j in range(0,5):
            for i in range(0,len(xmf3)):
                mainseqsf_K13[z,j,i] = mainseqsf_modelvar[j,st+i]

    idx = 0
    idx_modelvar = 0
    z = 0
    subplot = (311)

    ax = fig.add_subplot(subplot)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
 
    #Predicted relation
    ind = np.where(mainseqsf[idx,0,:] != 0)
    xplot = xmf[ind]
    yplot = mainseqsf[idx,0,ind]
    errdn = mainseqsf_1s[idx,1,ind]
    errup = mainseqsf_1s[idx,2,ind]
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='grey', alpha=1,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='grey', alpha=1,interpolate=True)
    ax.plot(xplot,yplot[0],'k', linestyle='solid', label="Shark")
 
    ind = np.where(mainseqsf_GD14[idx_modelvar,0,:] != 0)
    xplot = xmf3[ind]
    yplot = mainseqsf_GD14[idx_modelvar,0,ind]
    errdn = mainseqsf_GD14[idx_modelvar,3,ind]
    errup = mainseqsf_GD14[idx_modelvar,4,ind]
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='r', alpha=0.3,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='r', alpha=0.3,interpolate=True)
    ax.plot(xplot,yplot[0],'r', linestyle='dashed', label="GD14")
 
    ind = np.where(mainseqsf_KMT09[idx_modelvar,0,:] != 0)
    xplot = xmf3[ind]
    yplot = mainseqsf_KMT09[idx_modelvar,0,ind]
    errdn = mainseqsf_KMT09[idx_modelvar,3,ind]
    errup = mainseqsf_KMT09[idx_modelvar,4,ind]
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='b', alpha=0.3,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='b', alpha=0.3,interpolate=True)
    ax.plot(xplot,yplot[0],'b', linestyle='dotted', label="KMT09")
 
    ind = np.where(mainseqsf_K13[idx_modelvar,0,:] != 0)
    xplot = xmf3[ind]
    yplot = mainseqsf_K13[idx_modelvar,0,ind]
    errdn = mainseqsf_K13[idx_modelvar,3,ind]
    errup = mainseqsf_K13[idx_modelvar,4,ind]
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='gold', alpha=0.5,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='gold', alpha=0.5,interpolate=True)
    ax.plot(xplot,yplot[0],'gold', linestyle='dotted', label="K13")
 
    if(idx == 0):
            colors = ['k','r','b','gold']
            common.prepare_legend(ax, colors, bbox_to_anchor=(0.6, -0.05))
 
    #HI gas fraction plots
    xtit = "" #$\\rm log_{10} (\\rm M_{\\rm star}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\\rm M_{\\rm HI}/M_{\\star})$"
    xmin, xmax, ymin, ymax = 8, 12, -3, 2
    xleg = xmax - 0.3 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(312)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))

    mainseqgas_modelvar = np.zeros(shape = (6, 135))
    mainseqgas_modelvar[0,:], mainseqgas_modelvar[1,:], mainseqgas_modelvar[2,:], mainseqgas_modelvar[3,:], mainseqgas_modelvar[4,:], mainseqgas_modelvar[5,:] = common.load_observation(obsdir, 'Models/SharkVariations/GasMstars_OtherModels.dat', [0,1,2,3,4,5])

    mainseqgas_GD14  = np.zeros(shape = (6, len(xmf3)))
    mainseqgas_KMT09 = np.zeros(shape = (6, len(xmf3)))
    mainseqgas_K13   = np.zeros(shape = (6, len(xmf3)))

    st = 0
    for j in range(0,6):
        for i in range(0,len(xmf3)):
            mainseqgas_GD14[j,i] = mainseqgas_modelvar[j,st+i]

    st = len(xmf3) 
    for j in range(0,6):
        for i in range(0,len(xmf3)):
            mainseqgas_KMT09[j,i] = mainseqgas_modelvar[j,st+i]

    st = len(xmf3) * 2
    for j in range(0,6):
        for i in range(0,len(xmf3)):
            mainseqgas_K13[j,i] = mainseqgas_modelvar[j,st+i]

    #Predicted relation
    idx = 0
   
    ind = np.where(mainseqHI[idx,0,:] != 0)
    xplot = xmf[ind]
    yplot = mainseqHI[idx,0,ind]
    errdn = mainseqHI[idx,1,ind]
    errup = mainseqHI[idx,2,ind]
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='grey', alpha=1,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='grey', alpha=1,interpolate=True)
    ax.plot(xplot,yplot[0],'k', linestyle='solid')

    idx_modelvar = 0
    ind = np.where(mainseqgas_GD14[0+idx_modelvar,:] != 0)
    xplot = xmf3[ind]
    yplot = mainseqgas_GD14[0+idx_modelvar,ind]
    errdn = mainseqgas_GD14[1+idx_modelvar,ind]
    errup = mainseqgas_GD14[2+idx_modelvar,ind]
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='r', alpha=0.3,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='r', alpha=0.3,interpolate=True)
    ax.plot(xplot,yplot[0],'r', linestyle='dashed')

    ind = np.where(mainseqgas_KMT09[0+idx_modelvar,:] != 0)
    xplot = xmf3[ind]
    yplot = mainseqgas_KMT09[0+idx_modelvar,ind]
    errdn = mainseqgas_KMT09[1+idx_modelvar,ind]
    errup = mainseqgas_KMT09[2+idx_modelvar,ind]
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='b', alpha=0.3,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='b', alpha=0.3,interpolate=True)
    ax.plot(xplot,yplot[0],'b', linestyle='dotted')

    ind = np.where((mainseqgas_K13[0+idx_modelvar,:] != 0) & (xmf3 < 11))
    xplot = xmf3[ind]
    yplot = mainseqgas_K13[0+idx_modelvar,ind]
    errdn = mainseqgas_K13[1+idx_modelvar,ind]
    errup = mainseqgas_K13[2+idx_modelvar,ind]
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='gold', alpha=0.5,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='gold', alpha=0.5,interpolate=True)
    ax.plot(xplot,yplot[0],'gold', linestyle='dotted')

    #HI gas fraction plots
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\\rm M_{\\rm H_2}/M_{\\star})$"

    xmin, xmax, ymin, ymax = 8., 12, -3, 0.5
    xleg = xmax - 0.3 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(313)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))


    #Predicted relation
    ind = np.where(mainseqH2[idx,0,:] != 0)
    xplot = xmf[ind]
    yplot = mainseqH2[idx,0,ind]
    errdn = mainseqH2[idx,1,ind]
    errup = mainseqH2[idx,2,ind]
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='grey', alpha=1,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='grey', alpha=1,interpolate=True)
    ax.plot(xplot,yplot[0],'k', linestyle='solid')

    idx_modelvar = 3
    ind = np.where(mainseqgas_GD14[0+idx_modelvar,:] != 0)
    xplot = xmf3[ind]
    yplot = mainseqgas_GD14[0+idx_modelvar,ind]
    errdn = mainseqgas_GD14[1+idx_modelvar,ind]
    errup = mainseqgas_GD14[2+idx_modelvar,ind]
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='r', alpha=0.3,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='r', alpha=0.3,interpolate=True)
    ax.plot(xplot,yplot[0],'r', linestyle='dashed')

    ind = np.where(mainseqgas_KMT09[0+idx_modelvar,:] != 0)
    xplot = xmf3[ind]
    yplot = mainseqgas_KMT09[0+idx_modelvar,ind]
    errdn = mainseqgas_KMT09[1+idx_modelvar,ind]
    errup = mainseqgas_KMT09[2+idx_modelvar,ind]
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='b', alpha=0.3,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='b', alpha=0.3,interpolate=True)
    ax.plot(xplot,yplot[0],'b', linestyle='dotted')

    ind = np.where((mainseqgas_K13[0+idx_modelvar,:] != 0))
    xplot = xmf3[ind]
    yplot = mainseqgas_K13[0+idx_modelvar,ind]
    errdn = mainseqgas_K13[1+idx_modelvar,ind]
    errup = mainseqgas_K13[2+idx_modelvar,ind]
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='gold', alpha=0.5,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='gold', alpha=0.5,interpolate=True)
    ax.plot(xplot,yplot[0],'gold', linestyle='dotted')


    common.savefig(outdir, fig, 'SFR_Mstars.pdf')


def plot_SFE_Mstars(plt, outdir, sfe, sfe_cen, sfe_sat):

    fig = plt.figure(figsize=(9.5,9.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
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
        common.errorbars(ax, xplot, yplot[0], errdn[0], errup[0], 'b', 'o', label="all galaxies", err_absolute=False)

        ind = np.where(sfe_cen[idx,0,:] != 0)
        xplot = xmf[ind]
        yplot = sfe_cen[idx,0,ind]
        errdn = sfe_cen[idx,1,ind]
        errup = sfe_cen[idx,2,ind]
        common.errorbars(ax, xplot, yplot[0], errdn[0], errup[0], 'g', 'o', markersize=2, label="centrals", err_absolute=False)

        ind = np.where(sfe_sat[idx,0,:] != 0)
        xplot = xmf[ind]
        yplot = sfe_sat[idx,0,ind]
        errdn = sfe_sat[idx,1,ind]
        errup = sfe_sat[idx,2,ind]
        common.errorbars(ax, xplot, yplot[0], errdn[0], errup[0], 'r', 'o', markersize=2, label="satellites", err_absolute=False)

        common.prepare_legend(ax, ['b', 'g', 'r'])

    common.savefig(outdir, fig, 'SFE_Mstars.pdf')


def plot_fmzr(plt, outdir, fmzr):

    fig = plt.figure(figsize=(5,4.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot}) - 0.66 log_{10}(\\rm SFR/M_{\odot} yr^{-1})$"
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


def plot_mzr_z0(plt, outdir, obsdir, h0, mzr_cen, mzr_sat, mszr, mszr_cen, mszr_sat, mzr):

    fig = plt.figure(figsize=(4.5,8))
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\\rm Z_{\\rm gas}/Z_{\odot})$"
    xmin, xmax, ymin, ymax = 8, 12, -2, 1

    ax = fig.add_subplot(211)
    plt.subplots_adjust(bottom=0.15, left=0.18)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    ind = np.where(mzr[0,0,:] != 0)
    yplot = (mzr[0,0,ind])
    errdn = (mzr[0,1,ind])
    errup = (mzr[0,2,ind])
    xplot = xmf[ind]

    #print("mass metallicity relation at z=0 for the gas")
    #for a,b,c,d in zip(xplot,yplot[0],errdn[0],errup[0]):
    #    print(a,b,c,d)

    ind = np.where(mzr_cen[0,0,:] != 0)
    yplot = (mzr_cen[0,0,ind])
    errdn = (mzr_cen[0,1,ind])
    errup = (mzr_cen[0,2,ind])
    xplot = xmf[ind]
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='b', alpha=0.4,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='b', alpha=0.4,interpolate=True)
    ax.plot(xplot,yplot[0],'b', linestyle='solid')

    ind = np.where(mzr_sat[0,0,:] != 0)
    yplot = (mzr_sat[0,0,ind])
    errdn = (mzr_sat[0,1,ind])
    errup = (mzr_sat[0,2,ind])
    xplot = xmf[ind]
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='r', alpha=0.4,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='r', alpha=0.4,interpolate=True)
    ax.plot(xplot,yplot[0],'r', linestyle='dashed')

    #MZR z=0
    corrzsun = 8.69 #solar oxygen abundance in units of 12 + log(O/H)
    hobs = 0.72
    #add cosmology correction plus IMF correction that goes into the stellar mass.
    corr_cos = np.log10(pow(hobs,2)/pow(h0,2)) - 0.09
    lm, mz  = common.load_observation(obsdir, 'MZR/MMR-Kewley08.dat', [0,1])
    #ax.plot(lm[0:11] + corr_cos,mz[0:11]- corrzsun,'grey', linestyle='solid', linewidth = 0.8, label='Kewley & Ellison 08')
    #ax.plot(lm[13:27] + corr_cos,mz[13:27]- corrzsun,'grey', linestyle='solid', linewidth = 0.6)
    #ax.plot(lm[119:134] + corr_cos,mz[119:134]- corrzsun,'grey', linestyle='solid', linewidth = 0.6)
    #ax.plot(lm[136:148] + corr_cos,mz[136:148]- corrzsun,'grey', linestyle='solid', linewidth = 0.6)

    lm, mz, mzdn, mzup = common.load_observation(obsdir, 'MZR/MMAdrews13.dat', [0,1,2,3])
    hobs = 0.7
    #add cosmology correction plus IMF correction that goes into the stellar mass.
    corr_cos = np.log10(pow(hobs,2)/pow(h0,2)) 
    #common.errorbars(ax, lm+ corr_cos, mz - corrzsun, mzdn - corrzsun, mzup - corrzsun, 'grey', 's', label='Andrews+13')
    #correction for Tremonti is the same.
    lm, mz, mzdn, mzup = common.load_observation(obsdir, 'MZR/Curti2020.dat', [0,1,2,3])
    common.errorbars(ax, lm+ corr_cos, mz - corrzsun, mzdn - corrzsun, mzup - corrzsun, 'grey', 'o', label="Curti+20")
 
    lm, mz = common.load_observation(obsdir, 'MZR/Curti2020_Compilation.dat', [0,1])
    ax.plot(lm+ corr_cos, mz - corrzsun, color='darkgreen', marker='*', ls='None', mfc='None', label="Davies+17")
   
    #ZTe = 0.293±0.023 log(M∗) + 5.575±0.192
    #ZTe = 0.332±0.021 log(M∗) + 5.242±0.171
    ind = np.where(xmf < 9.87)
    ax.plot(xmf[ind], 0.332 * xmf[ind] + 5.242 - corrzsun, linestyle = 'solid', color='darkblue', label='Yates+19')
    ax.plot(xmf[ind], 0.332 * xmf[ind] + 5.242 - corrzsun - 0.171, linestyle = 'solid', color='darkblue')
    ax.plot(xmf[ind], 0.332 * xmf[ind] + 5.242 - corrzsun + 0.171, linestyle = 'solid', color='darkblue')

    common.prepare_legend(ax, ['darkgreen','darkblue','grey','grey'], loc=4)


    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\\rm Z_{\\star}/Z_{\odot})$"
    xmin, xmax, ymin, ymax = 8, 12, -2, 1

    ax = fig.add_subplot(212)
    plt.subplots_adjust(bottom=0.15, left=0.18)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))

    #MZR z=0
    lm, mz, mzdn, mzup = common.load_observation(obsdir, 'MZR/MSZR-Gallazzi05.dat', [0,1,2,3])
    common.errorbars(ax, lm[0:7], mz[0:7], mzdn[0:7], mzup[0:7], 'grey', 'D', label='Kirby+13')
    common.errorbars(ax, lm[7:22], mz[7:22], mzdn[7:22], mzup[7:22], 'grey', 'o', label='Gallazzi+05')

    ind = np.where(mszr_cen[0,0,:] != 0)
    yplot = (mszr_cen[0,0,ind])
    errdn = (mszr_cen[0,1,ind])
    errup = (mszr_cen[0,2,ind])
    xplot = xmf[ind]
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='b', alpha=0.4,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='b', alpha=0.4,interpolate=True)
    ax.plot(xplot,yplot[0],'b', linestyle='solid',label="Shark star-forming")

    ind = np.where(mszr_sat[0,0,:] != 0)
    yplot = (mszr_sat[0,0,ind])
    errdn = (mszr_sat[0,1,ind])
    errup = (mszr_sat[0,2,ind])
    xplot = xmf[ind]
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='r', alpha=0.4,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='r', alpha=0.4,interpolate=True)
    ax.plot(xplot,yplot[0],'r', linestyle='dashed',label="Shark all")

    common.prepare_legend(ax, ['b', 'r', 'grey','grey'], loc=4)

    common.savefig(outdir, fig, 'mzr_z0.pdf')


def plot_sfr_mstars_z0(plt, outdir, obsdir, h0, sfr_seq, mainseqsf, sfr_hi):

    bin_it = functools.partial(us.wmedians, xbins=xmf, nmin=10)

    fig = plt.figure(figsize=(5,5))
    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm log_{10}(\\rm SFR/M_{\odot} yr^{-1})$"

    xmin, xmax, ymin, ymax = 8, 12.5, -3.2, 3.3
    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    ind = np.where(sfr_seq[1,:] < -3)
    sfr_seq[1,ind] = -3
    #predicted relation
    ind = np.where((sfr_seq[0,:] > 7) & (sfr_seq[0,:] < 13) & (sfr_seq[1,:] >-10) & (sfr_seq[1,:] < 10))
    xdata = sfr_seq[0,ind]
    ydata = sfr_seq[1,ind]
    us.density_contour(ax, xdata[0], ydata[0], 30, 30, cmap = 'cividis') #, **contour_kwargs)

    ind = np.where(sfr_seq[0,:] > 0)
    toplot = bin_it(x=sfr_seq[0,ind], y=sfr_seq[1,ind])
    ind = np.where(toplot[0,:] != 0)
    yp = toplot[0,ind] 
    ax.plot(xmf[ind], yp[0],color='k',linestyle='solid', linewidth = 2, label="Shark v2.0")

    xm, ym = common.load_observation(obsdir, 'Models/SharkVariations/SFRMstars_Lagos18.dat', [0,1])
    ax.plot(xm, ym, linestyle='dashed', linewidth = 2, color='black',label='Shark v1.1 (L18)')

    #print("will print median MS")
    #for a,b in zip(xmf, toplot[0]):
    #    print(a,b)

    #SFR relation z=0
    lm, SFR = common.load_observation(obsdir, 'SFR/Brinchmann04.dat', (0, 1))
    hobs = 0.7
    #add cosmology correction plus IMF correction that goes into the stellar mass.
    corr_cos = np.log10(pow(hobs,2)/pow(h0,2)) - 0.09
    #apply correction to both stellar mass and SFRs.
    ax.plot(lm[0:35] + corr_cos, SFR[0:35] + corr_cos, color='SandyBrown', linewidth = 4, linestyle='dashed', label='Brinchmann+04')
    #ax.plot(lm[36:70] + corr_cos, SFR[36:70] + corr_cos, color='PaleVioletRed',linewidth = 5, linestyle='dotted')
    #ax.plot(lm[71:len(SFR)] + corr_cos, SFR[71:len(SFR)] + corr_cos, color='PaleVioletRed',linewidth = 5, linestyle='dotted')

    xdataD16 = [9.3, 10.6]
    ydataD16 = [-0.39, 0.477]
    ax.plot(xdataD16,ydataD16, color='Crimson',linestyle='dashdot',linewidth = 4, label='Davies+16')

    #GAMA data at z<0.06
    #CATAID StellarMass_bestfit StellarMass_50 StellarMass_16 StellarMass_84 SFR_bestfit SFR_50 SFR_16 SFR_84 Zgas_bestfit Zgas_50 Zgas_16 Zgas_84 DustMass_bestfit DustMass_50 DustMass_16 DustMass_84 DustLum_50 DustLum_16 DustLum_84 uberID redshift
    ms_gama, sfr_gama = common.load_observation(obsdir, 'GAMA/ProSpect_Claudia.txt', [2,6])
    ind = np.where(sfr_gama < 1e-3)
    sfr_gama[ind] = 1e-3
    #ax.hexbin(np.log10(ms_gama), np.log10(sfr_gama), gridsize=(20,20), mincnt=5) #, cmap = 'plasma') #, **contour_kwargs)
    us.density_contour_reduced(ax, np.log10(ms_gama), np.log10(sfr_gama), 25, 25) #, **contour_kwargs)

    toplot = bin_it(x=np.log10(ms_gama), y=np.log10(sfr_gama))
    ind = np.where(toplot[0,:] != 0)
    yp = toplot[0,ind]
    yup = toplot[2,ind]
    ydn = toplot[1,ind]
    ax.plot(xmf[ind], yp[0],color='Maroon',linestyle='dashed', linewidth = 5, label="Bellstedt+20")
#ax.plot(xmf[ind], yp[0]+yup[0],color='PaleVioletRed',linestyle='dotted', linewidth = 5)
    #ax.plot(xmf[ind], yp[0]-ydn[0],color='PaleVioletRed',linestyle='dotted', linewidth = 5)

    # individual massive galaxies from Terrazas+17
    ms, sfr, upperlimflag = common.load_observation(obsdir, 'BHs/MBH_host_gals_Terrazas17.dat', [0,1,2])
    ind = np.where(ms > 11.3)
    ax.errorbar(ms[ind], sfr[ind], xerr=0.2, yerr=0.3, ls='None', mfc='None', ecolor = 'r', mec='r',marker='s',label="Terrazas+17")
    ind = np.where((upperlimflag == 1) & (ms > 11.3))
    for a,b in zip (ms[ind], sfr[ind]):
        ax.arrow(a, b, 0, -0.3, head_width=0.05, head_length=0.1, fc='r', ec='r')

    # Legend
    common.prepare_legend(ax, ['k','k','SandyBrown','Crimson','Maroon','r'], loc=2)
    plt.tight_layout()
    common.savefig(outdir, fig, 'SFR_Mstars_z0.pdf')


    #plot SFR vs Hneutral for central, massive galaxies
    fig = plt.figure(figsize=(5,5))
    ytit="$\\rm log_{10} (\\rm M_{\\rm neutral}/M_{\odot})$"
    xtit="$\\rm log_{10}(\\rm SFR/M_{\odot} yr^{-1})$"

    xmin, xmax, ymin, ymax = -3, 1.5, 8, 12.5
    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #predicted relation
    ind = np.where(sfr_hi[0,0,:] != 0)
    xplot = xsfr[ind]
    yplot = sfr_hi[0,0,ind]
    ydn = sfr_hi[0,1,ind]
    yup = sfr_hi[0,2,ind]

    ax.fill_between(xplot,yplot[0]+yup[0],yplot[0]-ydn[0], facecolor='grey', alpha=1,interpolate=True)
    ax.plot(xplot,yplot[0],color='k',linestyle='solid', linewidth = 1, label="Shark cens")


    # load observations
    SFR, Mneutral, Mneutralerr = common.load_observation(obsdir, 'SFR/UVIR_SFR_HI_Shi22.dat', (0, 5, 6))
    ax.errorbar(SFR, Mneutral, yerr=Mneutralerr, ls='None', mfc='None', ecolor = 'OrangeRed', mec='OrangeRed',marker='+',label="Shi+22 (UVIR)")
    SFR, Mneutral, Mneutralerr = common.load_observation(obsdir, 'SFR/Halpha_D400_SFR_HI_Shi22.dat', (0, 5, 6))
    ax.errorbar(SFR, Mneutral, yerr=Mneutralerr, ls='None', mfc='None', ecolor = 'gold', mec='gold',marker='+',label="Shi+22 (H$\\alpha$)")

 
    # Legend
    common.prepare_legend(ax, ['k', 'OrangeRed', 'gold'], loc=2)
    common.savefig(outdir, fig, 'SFR_HI_Centrals_z0.pdf')


def plot_passive_fraction(plt, outdir, obsdir, passive_fractions, hist_ssfr):

    fig = plt.figure(figsize=(5,4.5))
    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm passive\,fraction$"

    xmin, xmax, ymin, ymax = 8, 12, 0, 1
    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    ind = np.where(passive_fractions[0,0,:] >= 0)
    xplot = xmf2[ind]
    yplot = passive_fractions[0,0,ind]
    ax.plot(xplot,yplot[0],color='k',linestyle='solid', linewidth = 1, label="Shark all galaxies")
    ind = np.where(passive_fractions[0,1,:] >= 0)
    xplot = xmf2[ind]
    yplot = passive_fractions[0,1,ind]
    ax.plot(xplot,yplot[0],color='b',linestyle='solid', linewidth = 1, label="$\\rm M_{\\rm halo}<10^{12}\,M_{\odot}$")
    ind = np.where(passive_fractions[0,2,:] >= 0)
    xplot = xmf2[ind]
    yplot = passive_fractions[0,2,ind]
    ax.plot(xplot,yplot[0],color='r',linestyle='solid', linewidth = 1, label="$\\rm M_{\\rm halo}> 10^{12}\,M_{\odot}$")

    #passive fraction z=0
    lm, frac = common.load_observation(obsdir, 'SFR/PassiveFraction_Halpha_Davies16.dat', (0, 1))
    ax.plot(lm[0:6], frac[0:6], color='k', linewidth = 2, linestyle='dotted', label='Davies+16 all galaxies')
    ax.plot(lm[8:14], frac[8:14], color='b',linewidth = 2, linestyle='dotted', label='Davies+16 isolated')
    ax.plot(lm[16:20], frac[16:20], color='r',linewidth = 2, linestyle='dotted', label='Davies+16 groups')

    # Legend
    common.prepare_legend(ax, ['k','b','r', 'k', 'b', 'r'], loc=2)
    common.savefig(outdir, fig, 'passive_fraction_z0.pdf')

    #plot SSFR functions
    fig = plt.figure(figsize=(10,9))
    xtit="$\\rm log_{10} (\\rm SSFR/yr^{-1})$"
    ytit="$\\rm log_{10} (\\rm \\phi/Mpc^{-3}\\, dex^{-1})$"

    bins = ['$\\rm log_{10} (\\rm M_{\\rm star}/M_{\\odot})=[9.5,10]$', '$\\rm log_{10} (\\rm M_{\\rm star}/M_{\\odot})=[10,10.5]$','$\\rm log_{10} (\\rm M_{\\rm star}/M_{\\odot})=[10.5,11]$','$\\rm log_{10} (\\rm M_{\\rm star}/M_{\\odot})=[11,11.5]$']
    xmin, xmax, ymin, ymax = -12.9, -8, -6, -1

    xleg = xmin + 0.1 * (xmax - xmin)
    yleg = ymax - 0.2 * (ymax - ymin)

    def plot_katsianis_data(ax, i):
        ppos = (i + 1) * 2 - 1         
        perrpos = (i + 1) * 2
        lm, p, dp = common.load_observation(obsdir, 'SFR/Katsianis21_SSFR_functions.dat', [0,ppos,perrpos]) 
        ind = np.where(p > 0)
        xplot = lm[ind]
        yplot = np.log10(p[ind])
        yerrdn = yplot - np.log10(p[ind] - dp[ind])
        yerrup = np.log10(p[ind] + dp[ind]) - yplot
        ax.errorbar(xplot, yplot, yerr=[yerrdn, yerrup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Katsianis+21")


    subplots = (221, 222, 223, 224)
    for i,s in enumerate(subplots):

        ax = fig.add_subplot(s)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 0.5, 0.5))
        ax.text(xleg, yleg, bins[i])

        ind = np.where(hist_ssfr[0,i,:] != 0)
        xplot = xssfr[ind]-9
        yplot = hist_ssfr[0,i,ind]
        ax.plot(xplot,yplot[0],color='k',linestyle='solid', linewidth = 1, label="Shark")
        plot_katsianis_data(ax, i)
        # Legend
        common.prepare_legend(ax, ['k','grey'], loc=2)
    plt.tight_layout()
    common.savefig(outdir, fig, 'SSFR_functions_z0.pdf')

def calculate_HImass_above_threshold(MHI, rgas, thresh):

    rlow = np.log10(0.01*rgas) #from 0.1*rgas 
    rupp = np.log10(10*rgas) #up to 10 times r50 gas
    dr = 0.2
    rbins = np.arange(rlow,rupp,dr)
    xr    = rbins + dr/2.0
    re    = rgas/1.67 #for an exponential disk

    sigma_HI = MHI/(2.0 * np.pi * re**2.0) * np.exp(-10**xr/re)
    ind = np.where(sigma_HI > thresh)
    if(len(sigma_HI[ind] > 0)):
       rin = 10**xr[ind]
       MHIout = MHI * (1 - ( 1 + max(rin)/re) * np.exp(-max(rin)/re))
    else:
       MHIout = 0
    return MHIout
    

def calculate_SFR_inside(SFR, rs, thresh):

    return SFR * (1 - ( 1 + thresh/(rs/1.67)) * np.exp(-thresh/(rs/1.67)))


def prepare_data(hdf5_data, index, hist_smf, hist_smf_offset, hist_smf_cen, hist_smf_sat, 
                 hist_smf_30kpc, hist_HImf, hist_HImf2, hist_HImf_cen, hist_HImf_sat, hist_H2mf, 
                 hist_H2mf_cen, hist_H2mf_sat, mainseq, mainseqsf, sfe, mainseq_cen, 
                 mainseqsf_cen, sfe_cen, mainseq_sat, mainseqsf_sat, sfe_sat, mzr, 
                 fmzr, mzr_cen, mzr_sat, plotz, plotz_HImf, passive_fractions, hist_ssfr, 
                 mszr, mszr_cen, mszr_sat, mainseqsf_1s, mainseqHI, mainseqH2, hist_smf_err, 
                 hist_HImf_err, hist_smf_comp, sfr_hi):

    (h0, volh, sfr_disk, sfr_burst, mdisk, mbulge, rstar_disk, mBH, mHI, mH2, 
     mgas_disk, mHI_bulge, mH2_bulge, mgas_bulge, mgas_metals_disk, mgas_metals_bulge, 
     mstars_metals_disk, mstars_metals_bulge, typeg, mvir_hosthalo, rstar_bulge, 
     mbulge_mergers, mbulge_diskins, mbulge_mergers_assembly, mbulge_diskins_assembly,
     sAM_atomic_disk, vmax, rgas_disk) = hdf5_data

    rstar = (rstar_disk * mdisk + rstar_bulge * mbulge) / (mdisk + mbulge)

    ind = np.where((sfr_disk + sfr_burst > 0) & (rstar > 0))
    SFRin = calculate_SFR_inside((sfr_disk[ind] + sfr_burst[ind])/h0/1e9, rstar[ind]*1e3  / h0, 30.0)#enclosed in 30kpc
    print("cosmic SFR in 30kpc:", np.log10(sum(SFRin) / (volh / h0**3)), np.log10(sum(sfr_disk[ind] + sfr_burst[ind])/h0/1e9/ (volh / h0**3)))

    MHI_abovethresh = np.zeros(shape = len(mdisk))
    thresh = 1 * 1e6 #Msun/kpc^2
    for j in range(0, len(mdisk)):
        if(mHI[j] > 0):
           MHI_abovethresh[j] = calculate_HImass_above_threshold(mHI[j]/h0, rgas_disk[j]/h0*1e3, thresh)

    print(mHI/h0, MHI_abovethresh)
    zstar = (mstars_metals_disk + mstars_metals_bulge) / (mdisk + mbulge)
    rcomb = (rstar_disk * mdisk + rstar_bulge * mbulge) / (mdisk + mbulge) / h0 * 1e3
    mgas = mgas_disk+mgas_bulge
    mgas_metals = mgas_metals_disk + mgas_metals_bulge

    mass          = np.zeros(shape = len(mdisk))
    mass_30kpc    = np.zeros(shape = len(mdisk))
    massd_30kpc   = np.zeros(shape = len(mdisk))
    massb_30kpc   = np.zeros(shape = len(mdisk))
    mass_atom     = np.zeros(shape = len(mdisk))
    mass_atom_cen = np.zeros(shape = len(mdisk))
    mass_atom_sat = np.zeros(shape = len(mdisk))

    mass_mol = np.zeros(shape = len(mdisk))
    mass_mol_cen = np.zeros(shape = len(mdisk))
    mass_mol_sat = np.zeros(shape = len(mdisk))


    ind = np.where((np.isnan(mdisk+mbulge) == True) | (np.isinf(mdisk+mbulge) == True))
    if (len(mdisk[ind]) > 0):
         print("Number of galaxies with a stellar mass of NaN:", len(mdisk[ind]))
    else:
         print("All galaxies have well defined stellar mass")

    ind = np.where((np.isnan(mBH) == True) | (np.isinf(mBH) == True))
    if (len(mdisk[ind]) > 0):
         print("Number of galaxies with a BH mass of NaN:", len(mdisk[ind]))
    else:
         print("All galaxies have well defined BH mass")


    #select massive centrals
    bin_it_sfr = functools.partial(us.wmedians, xbins=xsfr, low_numbers=False, nmin=10)
    ind = np.where(((mdisk+mbulge)/h0 > 10**10.6) & ((mdisk+mbulge)/h0 > 1e11) & (typeg == 0) & (sfr_disk + sfr_burst > 0))
    sfr_cens = np.log10((sfr_disk[ind] + sfr_burst[ind])/h0) - 9  
    mneutral_cens = np.log10((mgas_disk[ind] + mgas_bulge[ind])/h0)
    sfr_hi[index,:] = bin_it_sfr(x=sfr_cens, y=mneutral_cens)

    ind = np.where((mdisk+mbulge) > 0.0)
    mass[ind] = np.log10(mdisk[ind] + mbulge[ind]) - np.log10(float(h0))
    logger.debug('number of galaxies with mstars>0 and max mass: %d, %d', len(mass[ind]), max(mass[ind]))
    
    H, _ = np.histogram(mass,bins=np.append(mbins,mupp))
    hist_smf[index,:] = hist_smf[index,:] + H
    ran_err = np.random.normal(0.0, 0.3, len(mass))
    mass_err = mass + ran_err
    H, _ = np.histogram(mass_err,bins=np.append(mbins,mupp))
    hist_smf_offset[index,:] = hist_smf_offset[index,:] + H

    #Calculate the stellar mass contained in 30pkpc, assuming an exponential profile for the disk and a Plummer profile for the bulge.
    ind = np.where((mdisk > 0.0)  & (rstar_disk > 0))
    massd_30kpc[ind] = mdisk[ind] * (1.0 - (1.0 + 30.0/(rstar_disk[ind]/1.67/h0 * MpcToKpc)) * np.exp(-30.0/(rstar_disk[ind]/1.67/h0 * MpcToKpc)))
    ind = np.where((mbulge > 0.0)  & (rstar_bulge > 0))
    massb_30kpc[ind] = mbulge[ind] * pow(30.0, 3.0) / pow((pow(30.0, 2.0) + pow(rstar_bulge[ind]/1.3/h0 * MpcToKpc, 2.0)), 3.0/2.0)

    ind = np.where((massd_30kpc + massb_30kpc) > 0)
    mass_30kpc[ind] = np.log10(massd_30kpc[ind] + massb_30kpc[ind]) - np.log10(float(h0))
    H, _ = np.histogram(mass_30kpc,bins=np.append(mbins,mupp))
    hist_smf_30kpc[index,:] = hist_smf_30kpc[index,:] + H

    #stellar mass functions separated into centrals and satellites
    ind = np.where(typeg == 0)
    H, _ = np.histogram(mass[ind],bins=np.append(mbins,mupp))
    hist_smf_cen[index,:] = hist_smf_cen[index,:] + H
    ind = np.where(typeg > 0)
    H, _ = np.histogram(mass[ind],bins=np.append(mbins,mupp))
    hist_smf_sat[index,:] = hist_smf_sat[index,:] + H

    #stellar mass functions by galaxy components (disks, bulges, bulges by mergers, bulges by disk ins)
    ind = np.where(mdisk > 0)
    H, _ = np.histogram(np.log10(mdisk[ind]),bins=np.append(mbins,mupp))
    hist_smf_comp[index,0,:] = hist_smf_comp[index,0,:] + H
    ind = np.where(mbulge > 0)
    H, _ = np.histogram(np.log10(mbulge[ind]),bins=np.append(mbins,mupp))
    hist_smf_comp[index,1,:] = hist_smf_comp[index,1,:] + H
    ind = np.where(mbulge_mergers > 0)
    H, _ = np.histogram(np.log10(mbulge_mergers[ind]),bins=np.append(mbins,mupp))
    hist_smf_comp[index,2,:] = hist_smf_comp[index,2,:] + H
    ind = np.where(mbulge_diskins > 0)
    H, _ = np.histogram(np.log10(mbulge_diskins[ind]),bins=np.append(mbins,mupp))
    hist_smf_comp[index,3,:] = hist_smf_comp[index,3,:] + H


    #gas mass functions and scaling relations
    ind = np.where((mHI+mHI_bulge) > 0)
    mass_atom[ind] = np.log10(mHI[ind]+mHI_bulge[ind]) - np.log10(float(h0)) + np.log10(XH)
    H_HI, _ = np.histogram(mass_atom,bins=np.append(mbins,mupp))
    hist_HImf[index,:] = hist_HImf[index,:] + H_HI
    ind = np.where(MHI_abovethresh+mHI_bulge > 0)
    mass_atom[ind] = np.log10(MHI_abovethresh[ind]+mHI_bulge[ind]) - np.log10(float(h0)) + np.log10(XH)
    H_HI, _ = np.histogram(mass_atom,bins=np.append(mbins,mupp))
    hist_HImf2[index,:] = hist_HImf2[index,:] + H_HI

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

    bin_it = functools.partial(us.wmedians, xbins=xmf, low_numbers=False, nmin=20)
    bin_it_2sigma = functools.partial(us.wmedians_2sigma, xbins=xmf)

    ind = np.where((sfr_disk+sfr_burst > 0) & (mdisk+mbulge > 0))
    mainseq[index,:] = bin_it(x=mass[ind], y=np.log10((sfr_disk[ind]+sfr_burst[ind])/(mdisk[ind]+mbulge[ind])))
    passive_fractions[index,0,:] = us.fractions(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)), y = np.log10((sfr_disk[ind]+sfr_burst[ind])/(mdisk[ind]+mbulge[ind])), xbins=xmf2, ythresh=-2.2)
    passive_fractions[index,0,:] = 1.0 - passive_fractions[index,0,:] 

    smbins = (9.75, 10.25, 10.75, 11.25)
    dsm = 0.25
    for i,sm in enumerate(smbins):
        ind = np.where((sfr_disk+sfr_burst > 0.03162277660168379) & (mass >= sm - dsm) & (mass < sm + dsm))
        H, _ = np.histogram(np.log10((sfr_disk[ind]+sfr_burst[ind])/(mdisk[ind]+mbulge[ind])),bins=np.append(ssfrbins,ssfrupp))
        hist_ssfr[index,i,:] = hist_ssfr[index,i,:] + H

    ind = np.where((sfr_disk+sfr_burst > 0) & (mdisk+mbulge > 0) & (mvir_hosthalo < 1e11))
    passive_fractions[index,1,:] = us.fractions(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)), y = np.log10((sfr_disk[ind]+sfr_burst[ind])/(mdisk[ind]+mbulge[ind])), xbins=xmf2, ythresh=-2.2)
    passive_fractions[index,1,:] = 1.0 - passive_fractions[index,1,:] 

    ind = np.where((sfr_disk+sfr_burst > 0) & (mdisk+mbulge > 0) & (mvir_hosthalo >= 1e11))
    passive_fractions[index,2,:] = us.fractions(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)), y = np.log10((sfr_disk[ind]+sfr_burst[ind])/(mdisk[ind]+mbulge[ind])), xbins=xmf2, ythresh=-2.2)
    passive_fractions[index,2,:] = 1.0 - passive_fractions[index,2,:] 

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

    ind = np.where(mstars_metals_disk+mstars_metals_bulge > 0.0)
    mszr[index,:] = bin_it(x=mass[ind], y=np.log10(((mstars_metals_disk[ind]+mstars_metals_bulge[ind])/(mdisk[ind]+mbulge[ind])/Zsun)))

    ind = np.where((mgas_metals > 0.0) & (mgas > 0) & (sfr_disk+sfr_burst > 0))
    fmzr[index,:] = bin_it(x=mass[ind]-0.66*np.log10((sfr_disk[ind]+sfr_burst[ind])/h0/GyrToYr),
                           y=np.log10((mgas_metals[ind]/mgas[ind]/Zsun)))

    ind = np.where((mgas_metals > 0.0)  & (mgas > 1e5) & ((sfr_disk+sfr_burst)/(mdisk+mbulge) >= 1e-2))
    mzr_cen[index,:] = bin_it(x=mass[ind], y=np.log10((mgas_metals[ind]/mgas[ind]/Zsun)))
    ind = np.where((mstars_metals_disk+mstars_metals_bulge > 0.0) & (typeg == 0) )
    mszr_cen[index,:] = bin_it(x=mass[ind], y=np.log10(((mstars_metals_disk[ind]+mstars_metals_bulge[ind])/(mdisk[ind]+mbulge[ind])/Zsun)))

    ind = np.where((mgas_metals > 0.0) & (mgas > 1e5) & (mass > 8))
    mzr_sat[index,:] = bin_it(x=mass[ind], y=np.log10((mgas_metals[ind]/mgas[ind]/Zsun)))
    ind = np.where((mstars_metals_disk+mstars_metals_bulge > 0.0) & (typeg > 0) & (mass > 8))
    mszr_sat[index,:] = bin_it(x=mass[ind], y=np.log10(((mstars_metals_disk[ind]+mstars_metals_bulge[ind])/(mdisk[ind]+mbulge[ind])/Zsun)))

    ind = np.where((sfr_disk+sfr_burst > 0) & (mdisk+mbulge > 0) & ((sfr_disk+sfr_burst)/(mdisk+mbulge) > 0))
    mainseqsf[index,:] = bin_it_2sigma(x=mass[ind], y=np.log10((sfr_disk[ind]+sfr_burst[ind])/h0/GyrToYr))
    mainseqsf_1s[index,:] = bin_it(x=mass[ind], y=np.log10((sfr_disk[ind]+sfr_burst[ind])/h0/GyrToYr))
    mainseqHI[index,:] = bin_it(x=mass[ind], y=np.log10((mHI[ind]+mHI_bulge[ind])/(mdisk[ind]+mbulge[ind])))
    mainseqH2[index,:] = bin_it(x=mass[ind], y=np.log10((mH2[ind]+mH2_bulge[ind])/(mdisk[ind]+mbulge[ind])))

    if volh > 0:
        vol = volh/pow(h0,3.)  # In Mpc^3
        hist_smf_err[index,:]  = (hist_smf[index,:] - np.sqrt(hist_smf[index,:]))/vol/dm

        hist_smf[index,:]  = hist_smf[index,:]/vol/dm
        hist_smf_comp[index,:]  = hist_smf_comp[index,:]/vol/dm
        hist_smf_30kpc[index,:] = hist_smf_30kpc[index,:]/vol/dm
        hist_smf_offset[index,:] = hist_smf_offset[index,:]/vol/dm
        hist_smf_cen[index,:]  = hist_smf_cen[index,:]/vol/dm
        hist_smf_sat[index,:]  = hist_smf_sat[index,:]/vol/dm

        hist_HImf_err[index,:] = (hist_HImf[index,:] - np.sqrt(hist_HImf[index,:]))/vol/dm
        hist_HImf[index,:] = hist_HImf[index,:]/vol/dm
        hist_HImf2[index,:] = hist_HImf2[index,:]/vol/dm
        hist_HImf_cen[index,:] = hist_HImf_cen[index,:]/vol/dm
        hist_HImf_sat[index,:] = hist_HImf_sat[index,:]/vol/dm

        hist_H2mf[index,:] = hist_H2mf[index,:]/vol/dm
        hist_H2mf_cen[index,:] = hist_H2mf_cen[index,:]/vol/dm
        hist_H2mf_sat[index,:] = hist_H2mf_sat[index,:]/vol/dm
    
        hist_ssfr[index,:] = hist_ssfr[index,:]/vol/dssfr

        plotz[index]     = True
        plotz_HImf[index]= True
    else:
        plotz[index]     = False
        plotz_HImf[index]= False


    return mass

def main(modeldir, outdir, redshift_table, subvols, obsdir):

    zlist = (0, 0.5, 1, 2, 3, 4)
    #zlist = (0.005, 0.2, 0.5 , 0.8 , 1.1 , 1.5 , 2.2 , 2.9 , 3.9, 5.1)

    plt = common.load_matplotlib()

    mainseq     = np.zeros(shape = (len(zlist), 3, len(xmf)))
    mainseq_cen = np.zeros(shape = (len(zlist), 3, len(xmf)))
    mainseq_sat = np.zeros(shape = (len(zlist), 3, len(xmf)))

    mainseqsf     = np.zeros(shape = (len(zlist), 3, len(xmf)))
    mainseqsf_cen = np.zeros(shape = (len(zlist), 3, len(xmf)))
    mainseqsf_sat = np.zeros(shape = (len(zlist), 3, len(xmf)))
    mainseqsf_1s  = np.zeros(shape = (len(zlist), 3, len(xmf)))
    mainseqHI     = np.zeros(shape = (len(zlist), 3, len(xmf)))
    mainseqH2     = np.zeros(shape = (len(zlist), 3, len(xmf)))

    mzr         = np.zeros(shape = (len(zlist), 3, len(xmf)))
    mzr_cen     = np.zeros(shape = (len(zlist), 3, len(xmf)))
    mzr_sat     = np.zeros(shape = (len(zlist), 3, len(xmf)))
    mszr        = np.zeros(shape = (len(zlist), 3, len(xmf)))
    mszr_cen    = np.zeros(shape = (len(zlist), 3, len(xmf)))
    mszr_sat    = np.zeros(shape = (len(zlist), 3, len(xmf)))

    fmzr        = np.zeros(shape = (len(zlist), 3, len(xmf)))

    sfe         = np.zeros(shape = (len(zlist), 3, len(xmf)))
    sfe_cen     = np.zeros(shape = (len(zlist), 3, len(xmf)))
    sfe_sat     = np.zeros(shape = (len(zlist), 3, len(xmf)))

    sfr_hi      = np.zeros(shape = (len(zlist), 3, len(xsfr))) 
    passive_fractions = np.zeros(shape = (len(zlist), 3, len(xmf2)))

    # Histograms
    hist_smf       = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_30kpc = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_offset   = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_cen   = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_sat   = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_err   = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_comp  = np.zeros(shape = (len(zlist), 4, len(mbins)))

    plotz = np.empty(shape=(len(zlist)), dtype=np.bool_)
    hist_HImf = np.zeros(shape = (len(zlist), len(mbins)))
    hist_HImf2 = np.zeros(shape = (len(zlist), len(mbins)))
    hist_HImf_cen = np.zeros(shape = (len(zlist), len(mbins)))
    hist_HImf_sat = np.zeros(shape = (len(zlist), len(mbins)))
    plotz_HImf = np.empty(shape=(len(zlist)), dtype=np.bool_)
    hist_H2mf = np.zeros(shape = (len(zlist), len(mbins)))
    hist_H2mf_cen = np.zeros(shape = (len(zlist), len(mbins)))
    hist_H2mf_sat = np.zeros(shape = (len(zlist), len(mbins)))
    hist_HImf_err = np.zeros(shape = (len(zlist), len(mbins)))

    hist_ssfr = np.zeros(shape = (len(zlist), 4, len(ssfrbins)))

    fields = {'galaxies': ('sfr_disk', 'sfr_burst', 'mstars_disk', 'mstars_bulge',
                           'rstar_disk', 'm_bh', 'matom_disk', 'mmol_disk', 'mgas_disk',
                           'matom_bulge', 'mmol_bulge', 'mgas_bulge',
                           'mgas_metals_disk', 'mgas_metals_bulge',
                           'mstars_metals_disk', 'mstars_metals_bulge', 'type', 
                           'mvir_hosthalo', 'rstar_bulge', 'mstars_burst_mergers', 
                           'mstars_burst_diskinstabilities', 'mstars_bulge_mergers_assembly', 'mstars_bulge_diskins_assembly',
                           'specific_angular_momentum_disk_gas_atom', 'vmax_subhalo', 'rgas_disk')}

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        mass = prepare_data(hdf5_data, index, hist_smf, hist_smf_offset, hist_smf_cen,
                             hist_smf_sat, hist_smf_30kpc, hist_HImf, hist_HImf2, hist_HImf_cen, hist_HImf_sat,
                             hist_H2mf, hist_H2mf_cen, hist_H2mf_sat, mainseq, mainseqsf,
                             sfe, mainseq_cen, mainseqsf_cen, sfe_cen, mainseq_sat,
                             mainseqsf_sat, sfe_sat, mzr, fmzr, mzr_cen, mzr_sat, plotz,
                             plotz_HImf, passive_fractions, hist_ssfr, mszr, mszr_cen, 
                             mszr_sat, mainseqsf_1s, mainseqHI, mainseqH2, hist_smf_err, hist_HImf_err, 
                             hist_smf_comp, sfr_hi)

        h0 = hdf5_data[0]
        volh = hdf5_data[1]
        if index == 0:
            (sfr_disk, sfr_burst, mdisk, mbulge) = hdf5_data[2:6]
            sfr_seq = np.zeros(shape = (2, len(mdisk)))
            ind = np.where(sfr_disk + sfr_burst <= 0)
            sfr_disk[ind] = 2e-10 * GyrToYr * h0 #assume a minimum

            ind  = np.where((sfr_disk + sfr_burst > 0) & (mdisk + mbulge > 0))
            sfr_seq[0,ind] = mass[ind]
            sfr_seq[1,ind] = np.log10((sfr_disk[ind] + sfr_burst[ind]) / h0 / GyrToYr)

            ind = np.where(mass < 10.5)
            sfr_den = np.log10(sum(sfr_disk[ind] + sfr_burst[ind]) /  h0 / GyrToYr) + np.log10(h0**3.0 / volh)
            print("cosmic SFR density at z=0 < 10.5", sfr_den, volh)
            ind = np.where(mass < 11)
            sfr_den = np.log10(sum(sfr_disk[ind] + sfr_burst[ind]) /  h0 / GyrToYr) + np.log10(h0**3.0 / volh)
            print("cosmic SFR density at z=0 < 11", sfr_den, volh)
            ind = np.where(mass < 12)
            sfr_den = np.log10(sum(sfr_disk[ind] + sfr_burst[ind]) /  h0 / GyrToYr) + np.log10(h0**3.0 / volh)
            print("cosmic SFR density at z=0 < 12", sfr_den, volh)

    # This should be the same in all HDF5 files

    # Take logs
    def take_log(array):
        ind = np.where(array > 0.)
        array[ind] = np.log10(array[ind])

    take_log(hist_smf)
    take_log(hist_smf_comp)
    take_log(hist_smf_30kpc)
    take_log(hist_smf_cen)
    take_log(hist_smf_sat)
    take_log(hist_smf_offset)
    take_log(hist_HImf)
    take_log(hist_HImf2)
    take_log(hist_H2mf)
    take_log(hist_HImf_cen)
    take_log(hist_H2mf_cen)
    take_log(hist_HImf_sat)
    take_log(hist_H2mf_sat)
    take_log(hist_ssfr)

    plot_stellarmf_z(plt, outdir, obsdir, h0, plotz, hist_smf, hist_smf_cen, hist_smf_sat, hist_smf_offset, hist_smf_30kpc)
    plot_stellarmf_galcomponents(plt, outdir, obsdir, h0, plotz, hist_smf, hist_smf_comp)
    plot_stellarmf_z_molcomp(plt, outdir, obsdir, h0, plotz, hist_smf)
    plot_HImf_z0(plt, outdir, obsdir, h0, plotz_HImf, hist_HImf, hist_HImf_cen, hist_HImf_sat, hist_HImf2)
    plot_H2mf_z0(plt, outdir, obsdir, h0, plotz_HImf, hist_H2mf, hist_H2mf_cen, hist_H2mf_sat)
    plot_SSFR_Mstars(plt, outdir, mainseq, mainseq_cen, mainseq_sat)
    plot_mzr(plt, outdir, obsdir, h0, mzr, mzr_cen, mzr_sat)
    plot_SFR_Mstars(plt, outdir, obsdir, mainseqsf, mainseqsf_cen, mainseqsf_sat, mainseqsf_1s, mainseqHI, mainseqH2)
    plot_SFE_Mstars(plt, outdir, sfe, sfe_cen, sfe_sat)
    plot_fmzr(plt, outdir, fmzr)
    plot_mzr_z0(plt, outdir, obsdir, h0, mzr_cen, mzr_sat, mszr, mszr_cen, mszr_sat, mzr)
    plot_sfr_mstars_z0(plt, outdir, obsdir, h0, sfr_seq, mainseqsf, sfr_hi)
    plot_passive_fraction(plt, outdir, obsdir, passive_fractions, hist_ssfr) 


if __name__ == '__main__':
    main(*common.parse_args())
