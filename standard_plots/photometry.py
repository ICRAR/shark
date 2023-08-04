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
"""HMF plots"""

import numpy as np
import os

import common


##################################

# Constants
GyrToYr = 1e9
Zsun = 0.0127
XH = 0.72
PI = 3.141592654
MpcToKpc = 1e3
c_light = 299792458.0 #m/s

# Mass function initialization

mlow = -30 + 5.0 * np.log10(0.677)
mupp = -10 + 5.0 * np.log10(0.677)
dm = 0.5
mbins = np.arange(mlow,mupp,dm)
xlf   = mbins + dm/2.0

def plot_flux_contributions(plt, outdir, obsdir, h0, fdisk_emission, fbulge_m_emission, fbulge_d_emission):

    hcorr = 5.0*np.log10(h0)
    xlf_obs  = xlf - hcorr
 
    xtit="$\\rm mag-5log(h) (AB)$"
    ytit="$\\rm log_{10}(f_{\\rm cont})$"
    fig = plt.figure(figsize=(8,8))

    subplots = (331, 332, 333, 334, 335, 336, 337, 338, 339)
    idx = (0, 1, 2, 3, 4, 5, 6, 7, 8)
    bands = (0, 2, 4, 6, 10, 16, 19, 22, 26)
    labels= ('GALEX FUV', 'SDSS u', 'SDSS r', 'SDSS z', 'VISTA K', 'IRAC 5.8', 'P70', 'S250', 'JCTM850')
   
    for subplot, idx, b in zip(subplots, idx, bands):

        if(idx <= 2):
           xmin, xmax, ymin, ymax = -26, -13, -2, 0.05
           xleg = xmax - 0.51 * (xmax-xmin)
           yleg = ymin + 0.1 * (ymax-ymin)
        elif(idx > 2 and idx <= 5):
           xmin, xmax, ymin, ymax = -28, -13, -2, 0.05
           xleg = xmax - 0.51 * (xmax-xmin)
           yleg = ymin + 0.1 * (ymax-ymin)
        elif(idx > 5):
           xmin, xmax, ymin, ymax = -30, -13, -2, 0.05
           xleg = xmax - 0.51 * (xmax-xmin)
           yleg = ymin + 0.1 * (ymax-ymin)

        ax = fig.add_subplot(subplot)
        if (idx == 0 or idx == 3 or idx == 6):
            ytitplot = ytit
        else:
            ytitplot = ' '
        if (idx > 5):
            xtitplot = xtit
        else:
            xtitplot = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytitplot, locators=(4, 4, 1, 1))
        ax.text(xleg,yleg, labels[idx], fontsize=12)

        #Predicted LF
        if(idx == 0):
            ind = np.where(fdisk_emission[0,b,:] > 0.)
            y = np.log10(fdisk_emission[0,b,ind])
            ax.plot(xlf_obs[ind],y[0],'b', linestyle='dotted', label ='disks')

            ind = np.where(fbulge_m_emission[0,b,:] > 0.)
            y = np.log10(fbulge_m_emission[0,b,ind])
            ax.plot(xlf_obs[ind],y[0],'r', linestyle='dashed', label ='merger-driven bulges')

            ind = np.where(fbulge_d_emission[0,b,:] > 0.)
            y = np.log10(fbulge_d_emission[0,b,ind])
            ax.plot(xlf_obs[ind],y[0],'LightSalmon', linestyle='dashdot', label ='disk-ins-driven bulges')

        else:
            ind = np.where(fdisk_emission[0,b,:] > 0.)
            y = np.log10(fdisk_emission[0,b,ind])
            ax.plot(xlf_obs[ind],y[0],'b', linestyle='dotted')

            ind = np.where(fbulge_m_emission[0,b,:] > 0.)
            y = np.log10(fbulge_m_emission[0,b,ind])
            ax.plot(xlf_obs[ind],y[0],'r', linestyle='dashed')

            ind = np.where(fbulge_d_emission[0,b,:] > 0.)
            y = np.log10(fbulge_d_emission[0,b,ind])
            ax.plot(xlf_obs[ind],y[0],'LightSalmon', linestyle='dashdot')

        if idx == 0:
            common.prepare_legend(ax, ['b','r','LightSalmon'], bbox_to_anchor=[0,1.0])

    common.savefig(outdir, fig, "fractions_luminosity_contribution.pdf")

    #z=2
    indexz = 3 
    fig = plt.figure(figsize=(8,8))
    subplots = (331, 332, 333, 334, 335, 336, 337, 338, 339)
    idx = (0, 1, 2, 3, 4, 5, 6, 7, 8)
    bands = (0, 2, 4, 6, 10, 16, 19, 22, 26)
    labels= ('GALEX FUV', 'SDSS u', 'SDSS r', 'SDSS z', 'VISTA K', 'IRAC 5.8', 'P70', 'S250', 'JCTM850')

    for subplot, idx, b in zip(subplots, idx, bands):

        if(idx <= 2):
           xmin, xmax, ymin, ymax = -26, -13, -2, 0.05
           xleg = xmax - 0.51 * (xmax-xmin)
           yleg = ymin + 0.1 * (ymax-ymin)
        elif(idx > 2 and idx <= 5):
           xmin, xmax, ymin, ymax = -28, -13, -2, 0.05
           xleg = xmax - 0.51 * (xmax-xmin)
           yleg = ymin + 0.1 * (ymax-ymin)
        elif(idx > 5):
           xmin, xmax, ymin, ymax = -30, -13, -2, 0.05
           xleg = xmax - 0.51 * (xmax-xmin)
           yleg = ymin + 0.1 * (ymax-ymin)

        ax = fig.add_subplot(subplot)
        if (idx == 0 or idx == 3 or idx == 6):
            ytitplot = ytit
        else:
            ytitplot = ' '
        if (idx > 5):
            xtitplot = xtit
        else:
            xtitplot = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytitplot, locators=(4, 4, 1, 1))
        ax.text(xleg,yleg, labels[idx], fontsize=12)

        #Predicted LF
        if(idx == 0):
            ind = np.where(fdisk_emission[indexz,b,:] > 0.)
            y = np.log10(fdisk_emission[indexz,b,ind])
            ax.plot(xlf_obs[ind],y[0],'b', linestyle='dotted', label ='disks')

            ind = np.where(fbulge_m_emission[indexz,b,:] > 0.)
            y = np.log10(fbulge_m_emission[indexz,b,ind])
            ax.plot(xlf_obs[ind],y[0],'r', linestyle='dashed', label ='merger-driven bulges')

            ind = np.where(fbulge_d_emission[indexz,b,:] > 0.)
            y = np.log10(fbulge_d_emission[indexz,b,ind])
            ax.plot(xlf_obs[ind],y[0],'LightSalmon', linestyle='dashdot', label ='disk-ins-driven bulges')

        else:
            ind = np.where(fdisk_emission[indexz,b,:] > 0.)
            y = np.log10(fdisk_emission[indexz,b,ind])
            ax.plot(xlf_obs[ind],y[0],'b', linestyle='dotted')

            ind = np.where(fbulge_m_emission[indexz,b,:] > 0.)
            y = np.log10(fbulge_m_emission[indexz,b,ind])
            ax.plot(xlf_obs[ind],y[0],'r', linestyle='dashed')

            ind = np.where(fbulge_d_emission[indexz,b,:] > 0.)
            y = np.log10(fbulge_d_emission[indexz,b,ind])
            ax.plot(xlf_obs[ind],y[0],'LightSalmon', linestyle='dashdot')

        if idx == 0:
            common.prepare_legend(ax, ['b','r','LightSalmon'], bbox_to_anchor=[0,1.0])

    common.savefig(outdir, fig, "fractions_luminosity_contribution_z2.pdf")

    #z=3
    indexz = 4
    fig = plt.figure(figsize=(8,8))
    subplots = (331, 332, 333, 334, 335, 336, 337, 338, 339)
    idx = (0, 1, 2, 3, 4, 5, 6, 7, 8)
    bands = (0, 2, 4, 6, 10, 16, 19, 22, 26)
    labels= ('GALEX FUV', 'SDSS u', 'SDSS r', 'SDSS z', 'VISTA K', 'IRAC 5.8', 'P70', 'S250', 'JCTM850')

    for subplot, idx, b in zip(subplots, idx, bands):

        if(idx <= 2):
           xmin, xmax, ymin, ymax = -26, -13, -2, 0.05
           xleg = xmax - 0.51 * (xmax-xmin)
           yleg = ymin + 0.1 * (ymax-ymin)
        elif(idx > 2 and idx <= 5):
           xmin, xmax, ymin, ymax = -28, -13, -2, 0.05
           xleg = xmax - 0.51 * (xmax-xmin)
           yleg = ymin + 0.1 * (ymax-ymin)
        elif(idx > 5):
           xmin, xmax, ymin, ymax = -30, -13, -2, 0.05
           xleg = xmax - 0.51 * (xmax-xmin)
           yleg = ymin + 0.1 * (ymax-ymin)

        ax = fig.add_subplot(subplot)
        if (idx == 0 or idx == 3 or idx == 6):
            ytitplot = ytit
        else:
            ytitplot = ' '
        if (idx > 5):
            xtitplot = xtit
        else:
            xtitplot = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytitplot, locators=(4, 4, 1, 1))
        ax.text(xleg,yleg, labels[idx], fontsize=12)

        #Predicted LF
        if(idx == 0):
            ind = np.where(fdisk_emission[indexz,b,:] > 0.)
            y = np.log10(fdisk_emission[indexz,b,ind])
            ax.plot(xlf_obs[ind],y[0],'b', linestyle='dotted', label ='disks')

            ind = np.where(fbulge_m_emission[indexz,b,:] > 0.)
            y = np.log10(fbulge_m_emission[indexz,b,ind])
            ax.plot(xlf_obs[ind],y[0],'r', linestyle='dashed', label ='merger-driven bulges')

            ind = np.where(fbulge_d_emission[indexz,b,:] > 0.)
            y = np.log10(fbulge_d_emission[indexz,b,ind])
            ax.plot(xlf_obs[ind],y[0],'LightSalmon', linestyle='dashdot', label ='disk-ins-driven bulges')

        else:
            ind = np.where(fdisk_emission[indexz,b,:] > 0.)
            y = np.log10(fdisk_emission[indexz,b,ind])
            ax.plot(xlf_obs[ind],y[0],'b', linestyle='dotted')

            ind = np.where(fbulge_m_emission[indexz,b,:] > 0.)
            y = np.log10(fbulge_m_emission[indexz,b,ind])
            ax.plot(xlf_obs[ind],y[0],'r', linestyle='dashed')

            ind = np.where(fbulge_d_emission[indexz,b,:] > 0.)
            y = np.log10(fbulge_d_emission[indexz,b,ind])
            ax.plot(xlf_obs[ind],y[0],'LightSalmon', linestyle='dashdot')

        if idx == 0:
            common.prepare_legend(ax, ['b','r','LightSalmon'], bbox_to_anchor=[0,1.0])

    common.savefig(outdir, fig, "fractions_luminosity_contribution_z3.pdf")


def plot_lfs(plt, outdir, obsdir, h0, LFs_dust, LFs_nodust):

    hcorr = 5.0*np.log10(h0)
    xlf_obs  = xlf - hcorr
 
    xtit="$\\rm mag-5log(h) (AB)$"
    ytit="$\\rm log_{10}(\Phi/(0.5\\, {\\rm mag})/h^3 {\\rm Mpc}^{-3})$"

    xmin, xmax, ymin, ymax = -25, -13, -5, -1
    xleg = xmin + 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)

    fig = plt.figure(figsize=(12,12))

    subplots = (331, 332, 334, 335, 336, 337, 338)
    idx = (0, 1, 2, 3, 4, 5, 6, 7)
    bands = (0, 1, 2, 3, 4, 5, 6)
    labels= ('GALEX FUV', 'GALEX NUV', 'SDSS u', 'SDSS g', 'SDSS r', 'SDSS i', 'SDSS z')
    obs = ('lf1500', 'lf2300','lfu','lfg','lfr', 'lfi', 'lfz')
   
    for subplot, idx, b in zip(subplots, idx, bands):

        ax = fig.add_subplot(subplot)
        if (idx == 0 or idx == 2 or idx == 5):
            ytitplot = ytit
        else:
            ytitplot = ' '
        if (idx == 5 or idx == 6 or idx == 7):
            xtitplot = xtit
        else:
            xtitplot = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytitplot, locators=(2, 2, 1, 1))
        ax.text(xleg,yleg, labels[idx])

        file = obsdir+'/lf/'+obs[idx]+'_z0_driver12.data'
        lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
        indx = np.where(p > 0)
        yobs = np.log10(p[indx])
        ydn  = np.log10(p[indx]-dp[indx])
        yup  = np.log10(p[indx]+dp[indx])

        if(idx == 6):
            ax.errorbar(lm[indx], yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Driver+2012")
        else:
            ax.errorbar(lm[indx], yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

        #Predicted LF
        if(idx == 6):
            ind = np.where(LFs_dust[0,4,b,:] < 0.)
            y = LFs_dust[0,4,b,ind]
            ax.plot(xlf_obs[ind],y[0],'k', linewidth=3, label ='Shark')
            ind = np.where(LFs_nodust[0,4,b,:] < 0.)
            y = LFs_nodust[0,4,b,ind]
            ax.plot(xlf_obs[ind],y[0],'k', linewidth=1,  label ='Shark intrinsic')

            ind = np.where(LFs_dust[0,3,b,:] < 0.)
            y = LFs_dust[0,3,b,ind]
            ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted', label ='disks')
            ind = np.where(LFs_dust[0,1,b,:] < 0.)
            y = LFs_dust[0,1,b,ind]
            ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed', label ='merger-driven bulges')
            ind = np.where(LFs_dust[0,0,b,:] < 0.)
            y = LFs_dust[0,0,b,ind]
            ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot', label='disk-ins-driven bulges')
        else:
            ind = np.where(LFs_dust[0,4,b,:] < 0.)
            y = LFs_dust[0,4,b,ind]
            ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
            ind = np.where(LFs_nodust[0,4,b,:] < 0.)
            y = LFs_nodust[0,4,b,ind]
            ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)

            ind = np.where(LFs_dust[0,3,b,:] < 0.)
            y = LFs_dust[0,3,b,ind]
            ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
            ind = np.where(LFs_dust[0,1,b,:] < 0.)
            y = LFs_dust[0,1,b,ind]
            ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
            ind = np.where(LFs_dust[0,0,b,:] < 0.)
            y = LFs_dust[0,0,b,ind]
            ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')

        if idx == 6:
            common.prepare_legend(ax, ['k','k','b','r','LightSalmon','grey','grey'], bbox_to_anchor=[1.1,0.1])

    common.savefig(outdir, fig, "SDSS_luminosity_functions.pdf")

    fig = plt.figure(figsize=(12,8.5))

    subplots = (231, 232, 234, 235)
    idx = (0, 1, 2, 3)
    bands = (7, 8, 9, 10)
    labels= ('VISTA Y', 'VISTA J', 'VISTA H', 'VISTA K')
    obs = ('lfy', 'lfj','lfh','lfk')
    
    for subplot, idx, b in zip(subplots, idx, bands):

        ax = fig.add_subplot(subplot)
        if (idx == 0 or idx == 2):
            ytitplot = ytit
        else:
            ytitplot = ' '
        if (idx == 2 or idx == 3):
            xtitplot = xtit
        else:
            xtitplot = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytitplot, locators=(2, 2, 1, 1))
        ax.text(xleg,yleg, labels[idx])

        file = obsdir+'/lf/'+obs[idx]+'_z0_driver12.data'
        lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
        indx = np.where(p > 0)
        yobs = np.log10(p[indx])
        ydn  = np.log10(p[indx]-dp[indx])
        yup  = np.log10(p[indx]+dp[indx])

        if(idx == 0):
            ax.errorbar(lm[indx], yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Driver+2012")
        else:
            ax.errorbar(lm[indx], yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

        #Predicted LF
        if(idx == 0):
            ind = np.where(LFs_dust[0,4,b,:] < 0.)
            y = LFs_dust[0,4,b,ind]
            ax.plot(xlf_obs[ind],y[0],'k', linewidth=3, label ='Shark')
            ind = np.where(LFs_nodust[0,4,b,:] < 0.)
            y = LFs_nodust[0,4,b,ind]
            ax.plot(xlf_obs[ind],y[0],'k', linewidth=1,  label ='Shark intrinsic')

            ind = np.where(LFs_dust[0,3,b,:] < 0.)
            y = LFs_dust[0,3,b,ind]
            ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted', label ='disks')
            ind = np.where(LFs_dust[0,1,b,:] < 0.)
            y = LFs_dust[0,1,b,ind]
            ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed', label ='merger-driven bulges')
            ind = np.where(LFs_dust[0,0,b,:] < 0.)
            y = LFs_dust[0,0,b,ind]
            ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot', label='disk-ins-driven bulges')
        else:
            ind = np.where(LFs_dust[0,4,b,:] < 0.)
            y = LFs_dust[0,4,b,ind]
            ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
            ind = np.where(LFs_nodust[0,4,b,:] < 0.)
            y = LFs_nodust[0,4,b,ind]
            ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)

            ind = np.where(LFs_dust[0,3,b,:] < 0.)
            y = LFs_dust[0,3,b,ind]
            ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
            ind = np.where(LFs_dust[0,1,b,:] < 0.)
            y = LFs_dust[0,1,b,ind]
            ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
            ind = np.where(LFs_dust[0,0,b,:] < 0.)
            y = LFs_dust[0,0,b,ind]
            ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')

        if idx == 0:
            common.prepare_legend(ax, ['k','k','b','r','LightSalmon','grey','grey'], bbox_to_anchor=[2.2,-0.3])

    common.savefig(outdir, fig, "VISTA_luminosity_functions.pdf")

    ###################################################################
    fig = plt.figure(figsize=(12,8.5))

    xtit="$\\rm mag-5log(h) (AB)$"
    ytit="$\\rm log_{10}(\Phi/dex^{-1} h^3 {\\rm Mpc}^{-3})$"

    xmin, xmax, ymin, ymax = -28, -15, -5, -1
    xleg = xmin + 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)

    subplots = (231, 232, 234, 235)
    idx = (0, 1, 2, 3)
    bands = (12, 13, 15, 16)
    labels= ('IRAC 3.6', 'IRAC 4.5', 'IRAC 5.8', 'IRAC 8')
    obs = (3.6, 4.5, 5.8, 8.0) 

    vegacorr = (2.79, 3.26, 3.73, 4.40) #From Gillian Wilson webpage

    file = obsdir+'/lf/lf_IRAC_zLT0p6_dai2009.data'
    bandI,lm,p,dp = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
    indx  = np.where(p > 0)
    bandI = bandI[indx]
    lm    = lm[indx] - hcorr
    yobs  = np.log10(p[indx])
    ydn   = np.log10(p[indx]-dp[indx])
    yup   = np.log10(p[indx]+dp[indx])
   
    for subplot, idx, b in zip(subplots, idx, bands):

        ax = fig.add_subplot(subplot)
        if (idx == 0 or idx == 2):
            ytitplot = ytit
        else:
            ytitplot = ' '
        if (idx == 2 or idx == 3):
            xtitplot = xtit
        else:
            xtitplot = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytitplot, locators=(2, 2, 1, 1))
        ax.text(xleg,yleg, labels[idx])

        iband = np.where(bandI == obs[idx])
        if(idx == 0):
            ax.errorbar(lm[iband]+vegacorr[idx], yobs[iband], yerr=[yobs[iband]-ydn[iband],yup[iband]-yobs[iband]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Dai+2009")
        else:
            ax.errorbar(lm[iband]+vegacorr[idx], yobs[iband], yerr=[yobs[iband]-ydn[iband],yup[iband]-yobs[iband]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

        #Predicted LF
        if(idx == 0):
            ind = np.where(LFs_dust[0,4,b,:] < 0.)
            y = LFs_dust[0,4,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'k', linewidth=3, label ='Shark')
            ind = np.where(LFs_nodust[0,4,b,:] < 0.)
            y = LFs_nodust[0,4,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'k', linewidth=1,  label ='Shark intrinsic')

            #plot same line to get the right labels
            ax.plot(xlf_obs[0:2],xlf_obs[0:2],'SlateGray', linewidth=2, label ='Shark z=0.25')

            ind = np.where(LFs_dust[0,3,b,:] < 0.)
            y = LFs_dust[0,3,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted', label ='disks')
            ind = np.where(LFs_dust[0,1,b,:] < 0.)
            y = LFs_dust[0,1,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed', label ='merger-driven bulges')
            ind = np.where(LFs_dust[0,0,b,:] < 0.)
            y = LFs_dust[0,0,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot', label='disk-ins-driven bulges')
        else:
            ind = np.where(LFs_dust[0,4,b,:] < 0.)
            y = LFs_dust[0,4,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
            ind = np.where(LFs_nodust[0,4,b,:] < 0.)
            y = LFs_nodust[0,4,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)

            if(idx == 3):
               ind = np.where(LFs_dust[1,4,b,:] < 0.)
               y = LFs_dust[1,4,b,ind]-np.log10(dm)
               ax.plot(xlf_obs[ind],y[0],'SlateGray', linewidth=2)

            ind = np.where(LFs_dust[0,3,b,:] < 0.)
            y = LFs_dust[0,3,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
            ind = np.where(LFs_dust[0,1,b,:] < 0.)
            y = LFs_dust[0,1,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
            ind = np.where(LFs_dust[0,0,b,:] < 0.)
            y = LFs_dust[0,0,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')

        if idx == 0:
            common.prepare_legend(ax, ['k','k','SlateGray','b','r','LightSalmon','grey'], bbox_to_anchor=[2.2,-0.3])

    common.savefig(outdir, fig, "IRAC_luminosity_functions.pdf")

    # All NIR bands together
    xtit="$\\rm mag-5log(h) (AB)$"
    ytit="$\\rm log_{10}(\Phi/(0.5\\, {\\rm mag})/h^3 {\\rm Mpc}^{-3})$"
    ytit2="$\\rm log_{10}(\Phi/dex^{-1} h^3 {\\rm Mpc}^{-3})$"

    xmin, xmax, ymin, ymax = -25, -13, -5, -1
    xleg = xmin + 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)

    fig = plt.figure(figsize=(8.5,13))

    subplots = (421, 422, 423, 424, 425, 426, 427, 428)
    idx = (0, 1, 2, 3, 4, 5, 6, 7)
    bands = (7, 8, 9, 10, 12, 13, 15, 16) 
    obs = ('lf1500', 'lf2300','lfu','lfg','lfr', 'lfi', 'lfz')
    labels= ('VISTA Y', 'VISTA J', 'VISTA H', 'VISTA K','IRAC 3.6', 'IRAC 4.5', 'IRAC 5.8', 'IRAC 8')
    obs_vista = ('lfy', 'lfj','lfh','lfk')
    obs_irac = (3.6, 4.5, 5.8, 8.0)

    file = obsdir+'/lf/lf_IRAC_zLT0p6_dai2009.data'
    bandI,lm,p,dp = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
    indx  = np.where(p > 0)
    bandI = bandI[indx]
    lmir    = lm[indx] - hcorr
    yobsir  = np.log10(p[indx])
    ydnir   = np.log10(p[indx]-dp[indx])
    yupir   = np.log10(p[indx]+dp[indx])

    for subplot, idx, b in zip(subplots, idx, bands):
        ax = fig.add_subplot(subplot)
        if (idx == 0 or idx == 2):
            ytitplot = ytit
        elif (idx == 4 or idx == 6):
            ytitplot = ytit2
        else:
            ytitplot = ' '
        if (idx >= 6):
            xtitplot = xtit
        else:
            xtitplot = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytitplot, locators=(2, 2, 1, 1))
        ax.text(xleg,yleg, labels[idx])

        if(idx  <= 3):
            file = obsdir+'/lf/'+obs_vista[idx]+'_z0_driver12.data'
            lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
            indx = np.where(p > 0)
            yobs = np.log10(p[indx])
            ydn  = np.log10(p[indx]-dp[indx])
            yup  = np.log10(p[indx]+dp[indx])
            ax.errorbar(lm[indx], yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Driver+2012")
        else:
            iband = np.where(bandI == obs_irac[idx-4])
            ax.errorbar(lmir[iband]+vegacorr[idx-4], yobsir[iband], yerr=[yobsir[iband]-ydnir[iband],yupir[iband]-yobsir[iband]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='d',label="Dai+2009")

        #Predicted LF
        if(idx <= 3):
              ind = np.where(LFs_dust[0,4,b,:] < 0.)
              y = LFs_dust[0,4,b,ind]
              ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
              ind = np.where(LFs_nodust[0,4,b,:] < 0.)
              y = LFs_nodust[0,4,b,ind]
              ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)
 
              ind = np.where(LFs_dust[0,3,b,:] < 0.)
              y = LFs_dust[0,3,b,ind]
              ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
              ind = np.where(LFs_dust[0,1,b,:] < 0.)
              y = LFs_dust[0,1,b,ind]
              ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
              ind = np.where(LFs_dust[0,0,b,:] < 0.)
              y = LFs_dust[0,0,b,ind]
              ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')
        else:
              ind = np.where(LFs_dust[0,4,b,:] < 0.)
              y = LFs_dust[0,4,b,ind]-np.log10(dm)
              ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
              ind = np.where(LFs_nodust[0,4,b,:] < 0.)
              y = LFs_nodust[0,4,b,ind]-np.log10(dm)
              ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)
              
              ind = np.where(LFs_dust[0,3,b,:] < 0.)
              y = LFs_dust[0,3,b,ind]-np.log10(dm)
              ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
              ind = np.where(LFs_dust[0,1,b,:] < 0.)
              y = LFs_dust[0,1,b,ind]-np.log10(dm)
              ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
              ind = np.where(LFs_dust[0,0,b,:] < 0.)
              y = LFs_dust[0,0,b,ind]-np.log10(dm)
              ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')

        if(idx == 0 or idx == 4):
           common.prepare_legend(ax, ['k','k','b','r','LightSalmon','grey','grey'], loc='lower right')

    common.savefig(outdir, fig, "NIR_luminosity_functions.pdf")


    ########################################################
    fig = plt.figure(figsize=(12,12))
    xmin, xmax, ymin, ymax = -30, -16, -5, -1
    xleg = xmin + 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)


    subplots = (331, 332, 333, 334, 335, 336, 337, 338)
    idx = (0, 1, 2, 3, 4, 5, 6, 7)
    bands = (14, 17, 20, 21, 22, 23, 25, 26)
    labels= ('WISE 2', 'WISE 3', 'P100', 'P160', 'S250', 'S350', 'S500', 'JCTM850')
 
    for subplot, idx, b in zip(subplots, idx, bands):

        ax = fig.add_subplot(subplot)
        if (idx == 0 or idx == 3 or idx == 6):
            ytitplot = ytit
        else:
            ytitplot = ' '
        if (idx == 5 or idx == 6 or idx == 7):
            xtitplot = xtit
        else:
            xtitplot = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytitplot, locators=(2, 2, 1, 1))
        ax.text(xleg,yleg, labels[idx])

        #Predicted LF
        if(idx == 7):
            ind = np.where(LFs_dust[0,4,b,:] < 0.)
            y = LFs_dust[0,4,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'k', linewidth=3, label ='Shark z=0')

            #plot same line to get the right labels
            ind = np.where(LFs_dust[1,4,b,:] < 0.)
            y = LFs_dust[1,4,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'SlateGray', linewidth=2, label ='Shark z=0.25')

            ind = np.where(LFs_dust[2,4,b,:] < 0.)
            y = LFs_dust[2,4,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'Gray', linewidth=2, label ='Shark z=0.5')

            ind = np.where(LFs_dust[0,4,b,:] < 0.)
            y = LFs_dust[0,4,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'k', linewidth=1, label ='Shark intrinsic z=0')

            ind = np.where(LFs_dust[0,3,b,:] < 0.)
            y = LFs_dust[0,3,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted', label ='disks z=0')
            ind = np.where(LFs_dust[0,1,b,:] < 0.)
            y = LFs_dust[0,1,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed', label ='merger-driven bulges z=0')
            ind = np.where(LFs_dust[0,0,b,:] < 0.)
            y = LFs_dust[0,0,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot', label='disk-ins-driven bulges z=0')
        else:
            ind = np.where(LFs_dust[0,4,b,:] < 0.)
            y = LFs_dust[0,4,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
            if(idx == 0 or idx == 1):
               ind = np.where(LFs_nodust[0,4,b,:] < 0.)
               y = LFs_nodust[0,4,b,ind]-np.log10(dm)
               ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)
            if(idx >= 3):
               ind = np.where(LFs_dust[1,4,b,:] < 0.)
               y = LFs_dust[1,4,b,ind]-np.log10(dm)
               ax.plot(xlf_obs[ind],y[0],'SlateGray', linewidth=2)
               ind = np.where(LFs_dust[2,4,b,:] < 0.)
               y = LFs_dust[2,4,b,ind]-np.log10(dm)
               ax.plot(xlf_obs[ind],y[0],'Grey', linewidth=2)
 
            ind = np.where(LFs_dust[0,3,b,:] < 0.)
            y = LFs_dust[0,3,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
            ind = np.where(LFs_dust[0,1,b,:] < 0.)
            y = LFs_dust[0,1,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
            ind = np.where(LFs_dust[0,0,b,:] < 0.)
            y = LFs_dust[0,0,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')

        if idx == 3:
           hobs = 0.7
           corrhobs = np.log10(hobs/h0)

           file = obsdir+'/lf/L160microns.dat'
           lm,p,dpn,dpu = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
           dml       = 0.15
           xm        = [25.0,25.0]
           xm[1]     = xm[1] + dml
           tenpctocm = 3.086e19
           corrpc    = np.log10(tenpctocm)*2.0
           xobs      = -2.5*(lm[0:14]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
           xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
           dmm       = np.abs(xm[0] - xm[1])
           yobs      = np.log10(pow(10.0,p[0:14]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
           ydn       = np.log10(pow(10.0,dpn[0:14]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
           yup       = np.log10(pow(10.0,dpu[0:14]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
          
           ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*')
 
           dml       = 0.4
           xm        = [25.0,25.0]
           xm[1]     = xm[1] + dml
           tenpctocm = 3.086e19
           corrpc    = np.log10(tenpctocm)*2.0
           xobs      = -2.5*(lm[15:21]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
           xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
           dmm       = np.abs(xm[0] - xm[1])
           yobs      = np.log10(pow(10.0,p[15:21]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
           ydn       = np.log10(pow(10.0,dpn[15:21]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
           yup       = np.log10(pow(10.0,dpu[15:21]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
          
           ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='p')

        if idx == 4:
            hobs = 0.7
            corrhobs = np.log10(hobs/h0)

            file = obsdir+'/lf/lf250_dye10.data'
            lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
            lx        = np.log10(lm[0:7])
            dml       = np.abs(lx[1] - lx[0])
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(np.log10(lm[0:7])-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
            dmm       = np.abs(xobs[0] - xobs[1])
            yobs      = np.log10(p[0:7] * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            ydn       = np.log10((p[0:7]- dp[0:7]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            yup       = np.log10((p[0:7]+ dp[0:7]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
           
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')
           
            file = obsdir+'/lf/L250microns.dat'
            lm,p,dpn,dpu = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
            dml       = 0.15
            xm        = [25.0,25.0]
            xm[1]     = xm[1] + dml
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(lm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
            pxm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
            dmm       = np.abs(xm[0] - xm[1])
            yobs      = np.log10(pow(10.0,p) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            ydn       = np.log10(pow(10.0,dpn) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            yup       = np.log10(pow(10.0,dpu) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
           
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*')

        if idx == 5:
            hobs = 0.7
            corrhobs = np.log10(hobs/h0)

            dml       = 0.3
            xm        = [25.0,25.0]
            xm[1]     = xm[1] + dml
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(lm[10:19]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
            xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
            dmm       = np.abs(xm[0] - xm[1])
            yobs      = np.log10(pow(10.0,p[10:19]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            ydn       = np.log10(pow(10.0,dpn[10:19]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            yup       = np.log10(pow(10.0,dpu[10:19]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s')

        if idx == 6:
            hobs = 0.7
            corrhobs = np.log10(hobs/h0)

            file = obsdir+'/lf/L500microns.dat'
            lm,p,dpn,dpu = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
            dml       = 0.2
            xm        = [25.0,25.0]
            xm[1]     = xm[1] + dml
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(lm[0:11]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
            xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
            dmm       = np.abs(xm[0] - xm[1])
            yobs      = np.log10(pow(10.0,p[0:11]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            ydn       = np.log10(pow(10.0,dpn[0:11]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            yup       = np.log10(pow(10.0,dpu[0:11]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs

            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*')

            dml       = 0.3
            xm        = [25.0,25.0]
            xm[1]     = xm[1] + dml
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(lm[11:19]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
            xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
            dmm       = np.abs(xm[0] - xm[1])
            yobs      = np.log10(pow(10.0,p[11:19]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            ydn       = np.log10(pow(10.0,dpn[11:19]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            yup       = np.log10(pow(10.0,dpu[11:19]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs

            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s')

        if idx == 7:
            hobs = 0.75
            corrhobs = np.log10(hobs/h0)
            file = obsdir+'/lf/lf850_dunne2000.data'
            lm,p,dp   = np.loadtxt(file,usecols=[0,1,2],unpack=True)
            dml       = 0.24
            xm        = [25.0,25.0]
            xm[1]     = xm[1] + dml
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(lm[7:13]-corrpc+7.0) - 48.6 - hcorr + 5.0*corrhobs
            xm        = -2.5*(xm-corrpc+7.0) - 48.6 - hcorr  + 5.0*corrhobs
            dmm       = np.abs(xm[0] - xm[1])
            yobs      = np.log10(pow(10.0,np.log10(p[7:13])) * dml/dmm)- 3.0*np.log10(0.677) - 3.0 * corrhobs
            ydn       = np.log10(pow(10.0,np.log10(p[7:13] - dp[7:13])) * dml/dmm)- 3.0*np.log10(0.677)- 3.0 * corrhobs
            yup       = np.log10(pow(10.0,np.log10(p[7:13] + dp[7:13])) * dml/dmm)- 3.0*np.log10(0.677)- 3.0 * corrhobs
            
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='P',label="Dunne+00")
            
            file = obsdir+'/lf/lf850_vlahakis05_OS.data'
            lm,p,dp   = np.loadtxt(file,usecols=[0,1,2],unpack=True)
            dml       = 0.26
            xm        = [25.0,25.0]
            xm[1]     = xm[1] + dml
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(lm-corrpc+7.0) - 48.6 - hcorr + 5.0*corrhobs
            xm        = -2.5*(xm-corrpc+7.0) - 48.6 - hcorr + 5.0*corrhobs
            dmm       = np.abs(xm[0] - xm[1])
            yobs      = np.log10(pow(10.0,np.log10(p)) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            ydn       = np.log10(pow(10.0,np.log10(p - dp)) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            yup       = np.log10(pow(10.0,np.log10(p + dp)) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D',label="Vlahakis+05")

            xobs = np.zeros(shape = 2)
            xobs[0] = 1.0
            xobs[1] = 1.0
            yobs=xobs
            ydn=xobs
            yup=xobs
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='p',label="Patel+2013")
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Dye+2010")
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*',label="Marchetti+2016")
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s',label="Negrello+2013")

        if idx == 7:
            common.prepare_legend(ax, ['k','SlateGray','Gray','k','b','r','LightSalmon','grey','grey','grey','grey','grey','grey'], bbox_to_anchor=[1.1,-0.2])

    common.savefig(outdir, fig, "IR_luminosity_functions.pdf")

    ########################################################
    fig = plt.figure(figsize=(5,10))
    xmin, xmax, ymin, ymax = -30, -16, -5, -1
    xleg = xmin + 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)

    ytit="$\\rm log_{10}(\Phi/dex^{-1} h^3 {\\rm Mpc}^{-3})$"

    subplots = (311, 312, 313)
    idx = (0, 1, 2)
    bands = (21, 22,  26)
    labels= ('P160', 'S250', 'JCMT850')
    
    for subplot, idx, b in zip(subplots, idx, bands):

        ax = fig.add_subplot(subplot)
        ytitplot = ytit
        if (idx == 3):
            xtitplot = xtit
        else:
            xtitplot = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytitplot, locators=(2, 2, 1, 1))
        ax.text(xleg,yleg, labels[idx])

        #Predicted LF
        ind = np.where(LFs_dust[0,4,b,:] < 0.)
        y = LFs_dust[0,4,b,ind]-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
   
        #plot same line to get the right labels
        ind = np.where(LFs_dust[1,4,b,:] < 0.)
        y = LFs_dust[1,4,b,ind]-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)
   
        ind = np.where(LFs_dust[0,3,b,:] < 0.)
        y = LFs_dust[0,3,b,ind]-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
        ind = np.where(LFs_dust[0,1,b,:] < 0.)
        y = LFs_dust[0,1,b,ind]-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
        ind = np.where(LFs_dust[0,0,b,:] < 0.)
        y = LFs_dust[0,0,b,ind]-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')

        if idx == 0:
           hobs = 0.7
           corrhobs = np.log10(hobs/h0)

           file = obsdir+'/lf/L160microns.dat'
           lm,p,dpn,dpu = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)

           dml       = 0.4
           xm        = [25.0,25.0]
           xm[1]     = xm[1] + dml
           tenpctocm = 3.086e19
           corrpc    = np.log10(tenpctocm)*2.0
           xobs      = -2.5*(lm[15:21]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
           xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
           dmm       = np.abs(xm[0] - xm[1])
           yobs      = np.log10(pow(10.0,p[15:21]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
           ydn       = np.log10(pow(10.0,dpn[15:21]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
           yup       = np.log10(pow(10.0,dpu[15:21]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
          
           ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='p')

        if idx == 1:
            hobs = 0.7
            corrhobs = np.log10(hobs/h0)

            file = obsdir+'/lf/lf250_dye10.data'
            lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
            lx        = np.log10(lm[0:7])
            dml       = np.abs(lx[1] - lx[0])
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(np.log10(lm[0:7])-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
            dmm       = np.abs(xobs[0] - xobs[1])
            yobs      = np.log10(p[0:7] * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            ydn       = np.log10((p[0:7]- dp[0:7]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            yup       = np.log10((p[0:7]+ dp[0:7]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
           
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

        if idx == 2:
            hobs = 0.75
            corrhobs = np.log10(hobs/h0)
            file = obsdir+'/lf/lf850_dunne2000.data'
            lm,p,dp   = np.loadtxt(file,usecols=[0,1,2],unpack=True)
            dml       = 0.24
            xm        = [25.0,25.0]
            xm[1]     = xm[1] + dml
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(lm[7:13]-corrpc+7.0) - 48.6 - hcorr + 5.0*corrhobs
            xm        = -2.5*(xm-corrpc+7.0) - 48.6 - hcorr  + 5.0*corrhobs
            dmm       = np.abs(xm[0] - xm[1])
            yobs      = np.log10(pow(10.0,np.log10(p[7:13])) * dml/dmm)- 3.0*np.log10(0.677) - 3.0 * corrhobs
            ydn       = np.log10(pow(10.0,np.log10(p[7:13] - dp[7:13])) * dml/dmm)- 3.0*np.log10(0.677)- 3.0 * corrhobs
            yup       = np.log10(pow(10.0,np.log10(p[7:13] + dp[7:13])) * dml/dmm)- 3.0*np.log10(0.677)- 3.0 * corrhobs
            
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='P',label="Dunne+00")
            
            file = obsdir+'/lf/lf850_vlahakis05_OS.data'
            lm,p,dp   = np.loadtxt(file,usecols=[0,1,2],unpack=True)
            dml       = 0.26
            xm        = [25.0,25.0]
            xm[1]     = xm[1] + dml
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(lm-corrpc+7.0) - 48.6 - hcorr + 5.0*corrhobs
            xm        = -2.5*(xm-corrpc+7.0) - 48.6 - hcorr + 5.0*corrhobs
            dmm       = np.abs(xm[0] - xm[1])
            yobs      = np.log10(pow(10.0,np.log10(p)) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            ydn       = np.log10(pow(10.0,np.log10(p - dp)) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            yup       = np.log10(pow(10.0,np.log10(p + dp)) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D',label="Vlahakis+05")

            xobs = np.zeros(shape = 2)
            xobs[0] = 1.0
            xobs[1] = 1.0
            yobs=xobs
            ydn=xobs
            yup=xobs
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='p',label="Patel+2013")
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Dye+2010")

        if(idx == 2):
           common.prepare_legend(ax, ['grey','grey','grey','grey','grey','grey'], loc='lower left')
       
    common.savefig(outdir, fig, "IR_luminosity_functions_reliable_observations.pdf")

    ########################################################
    fig = plt.figure(figsize=(8.5,12))
    xmin, xmax, ymin, ymax = -30, -16, -5, -1
    xleg = xmin + 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)


    subplots = (321, 322, 323, 324, 325, 326)
    idx = (0, 1, 2, 3, 4, 5, 6, 7)
    bands = (21, 22, 23, 25, 26)
    labels= ('P160', 'S250', 'S350', 'S500', 'JCMT850')
    
    for subplot, idx, b in zip(subplots, idx, bands):

        ax = fig.add_subplot(subplot)
        if (idx == 0 or idx == 2 or idx == 4):
            ytitplot = ytit
        else:
            ytitplot = ' '
        if (idx == 3 or idx == 4):
            xtitplot = xtit
        else:
            xtitplot = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytitplot, locators=(2, 2, 1, 1))
        ax.text(xleg,yleg, labels[idx])

        #Predicted LF
        if(idx == 4):
            ind = np.where(LFs_dust[0,4,b,:] < 0.)
            y = LFs_dust[0,4,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'k', linewidth=3, label ='Shark z=0')

            #plot same line to get the right labels
            ind = np.where(LFs_dust[1,4,b,:] < 0.)
            y = LFs_dust[1,4,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'SlateGray', linewidth=1, label ='Shark z=0.5')

            ind = np.where(LFs_dust[0,3,b,:] < 0.)
            y = LFs_dust[0,3,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted', label ='disks z=0')
            ind = np.where(LFs_dust[0,1,b,:] < 0.)
            y = LFs_dust[0,1,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed', label ='merger-driven bulges z=0')
            ind = np.where(LFs_dust[0,0,b,:] < 0.)
            y = LFs_dust[0,0,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot', label='disk-ins-driven bulges z=0')
        else:
            ind = np.where(LFs_dust[0,4,b,:] < 0.)
            y = LFs_dust[0,4,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
            ind = np.where(LFs_dust[2,4,b,:] < 0.)
            y = LFs_dust[2,4,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],linewidth=1)
 
            ind = np.where(LFs_dust[0,3,b,:] < 0.)
            y = LFs_dust[0,3,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
            ind = np.where(LFs_dust[0,1,b,:] < 0.)
            y = LFs_dust[0,1,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
            ind = np.where(LFs_dust[0,0,b,:] < 0.)
            y = LFs_dust[0,0,b,ind]-np.log10(dm)
            ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')

        if idx == 0:
           hobs = 0.7
           corrhobs = np.log10(hobs/h0)

           file = obsdir+'/lf/L160microns.dat'
           lm,p,dpn,dpu = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
           dml       = 0.15
           xm        = [25.0,25.0]
           xm[1]     = xm[1] + dml
           tenpctocm = 3.086e19
           corrpc    = np.log10(tenpctocm)*2.0
           xobs      = -2.5*(lm[0:14]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
           xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
           dmm       = np.abs(xm[0] - xm[1])
           yobs      = np.log10(pow(10.0,p[0:14]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
           ydn       = np.log10(pow(10.0,dpn[0:14]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
           yup       = np.log10(pow(10.0,dpu[0:14]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
          
           ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*')
 
           dml       = 0.4
           xm        = [25.0,25.0]
           xm[1]     = xm[1] + dml
           tenpctocm = 3.086e19
           corrpc    = np.log10(tenpctocm)*2.0
           xobs      = -2.5*(lm[15:21]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
           xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
           dmm       = np.abs(xm[0] - xm[1])
           yobs      = np.log10(pow(10.0,p[15:21]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
           ydn       = np.log10(pow(10.0,dpn[15:21]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
           yup       = np.log10(pow(10.0,dpu[15:21]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
          
           ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='p')

        if idx == 1:
            hobs = 0.7
            corrhobs = np.log10(hobs/h0)

            file = obsdir+'/lf/lf250_dye10.data'
            lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
            lx        = np.log10(lm[0:7])
            dml       = np.abs(lx[1] - lx[0])
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(np.log10(lm[0:7])-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
            dmm       = np.abs(xobs[0] - xobs[1])
            yobs      = np.log10(p[0:7] * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            ydn       = np.log10((p[0:7]- dp[0:7]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            yup       = np.log10((p[0:7]+ dp[0:7]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
           
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')
         
            file = obsdir+'/lf/L250microns.dat'
            lm,p,dpn,dpu = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
            dml       = 0.15
            xm        = [25.0,25.0]
            xm[1]     = xm[1] + dml
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(lm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
            pxm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
            dmm       = np.abs(xm[0] - xm[1])
            yobs      = np.log10(pow(10.0,p) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            ydn       = np.log10(pow(10.0,dpn) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            yup       = np.log10(pow(10.0,dpu) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
           
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*')
  
        if idx == 2:
            hobs = 0.7
            corrhobs = np.log10(hobs/h0)

            file = obsdir+'/lf/L350microns.dat'
            lm,p,dpn,dpu = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
            dml       = 0.3
            xm        = [25.0,25.0]
            xm[1]     = xm[1] + dml
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(lm[10:19]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
            xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
            dmm       = np.abs(xm[0] - xm[1])
            yobs      = np.log10(pow(10.0,p[10:19]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            ydn       = np.log10(pow(10.0,dpn[10:19]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            yup       = np.log10(pow(10.0,dpu[10:19]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s')

        if idx == 3:
            hobs = 0.7
            corrhobs = np.log10(hobs/h0)

            file = obsdir+'/lf/L500microns.dat'
            lm,p,dpn,dpu = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
            dml       = 0.2
            xm        = [25.0,25.0]
            xm[1]     = xm[1] + dml
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(lm[0:11]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
            xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
            dmm       = np.abs(xm[0] - xm[1])
            yobs      = np.log10(pow(10.0,p[0:11]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            ydn       = np.log10(pow(10.0,dpn[0:11]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            yup       = np.log10(pow(10.0,dpu[0:11]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs

            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*')

            dml       = 0.3
            xm        = [25.0,25.0]
            xm[1]     = xm[1] + dml
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(lm[11:19]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
            xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
            dmm       = np.abs(xm[0] - xm[1])
            yobs      = np.log10(pow(10.0,p[11:19]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            ydn       = np.log10(pow(10.0,dpn[11:19]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            yup       = np.log10(pow(10.0,dpu[11:19]) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs

            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s')

        if idx == 4:
            hobs = 0.75
            corrhobs = np.log10(hobs/h0)
            file = obsdir+'/lf/lf850_dunne2000.data'
            lm,p,dp   = np.loadtxt(file,usecols=[0,1,2],unpack=True)
            dml       = 0.24
            xm        = [25.0,25.0]
            xm[1]     = xm[1] + dml
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(lm[7:13]-corrpc+7.0) - 48.6 - hcorr + 5.0*corrhobs
            xm        = -2.5*(xm-corrpc+7.0) - 48.6 - hcorr  + 5.0*corrhobs
            dmm       = np.abs(xm[0] - xm[1])
            yobs      = np.log10(pow(10.0,np.log10(p[7:13])) * dml/dmm)- 3.0*np.log10(0.677) - 3.0 * corrhobs
            ydn       = np.log10(pow(10.0,np.log10(p[7:13] - dp[7:13])) * dml/dmm)- 3.0*np.log10(0.677)- 3.0 * corrhobs
            yup       = np.log10(pow(10.0,np.log10(p[7:13] + dp[7:13])) * dml/dmm)- 3.0*np.log10(0.677)- 3.0 * corrhobs
            
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='P',label="Dunne+00")
            
            file = obsdir+'/lf/lf850_vlahakis05_OS.data'
            lm,p,dp   = np.loadtxt(file,usecols=[0,1,2],unpack=True)
            dml       = 0.26
            xm        = [25.0,25.0]
            xm[1]     = xm[1] + dml
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(lm-corrpc+7.0) - 48.6 - hcorr + 5.0*corrhobs
            xm        = -2.5*(xm-corrpc+7.0) - 48.6 - hcorr + 5.0*corrhobs
            dmm       = np.abs(xm[0] - xm[1])
            yobs      = np.log10(pow(10.0,np.log10(p)) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            ydn       = np.log10(pow(10.0,np.log10(p - dp)) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            yup       = np.log10(pow(10.0,np.log10(p + dp)) * dml/dmm)-3.0*np.log10(0.677)- 3.0 * corrhobs
            
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D',label="Vlahakis+05")

            xobs = np.zeros(shape = 2)
            xobs[0] = 1.0
            xobs[1] = 1.0
            yobs=xobs
            ydn=xobs
            yup=xobs
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='p',label="Patel+2013")
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Dye+2010")
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s',label="Negrello+2013")
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*',label="Marchetti+2016")

            common.prepare_legend(ax, ['k','k','b','r','LightSalmon','grey','grey','grey','grey','grey','grey'], bbox_to_anchor=[1.1,-0.2])

    common.savefig(outdir, fig, "IR_luminosity_functions_all_with_obs.pdf")


def prepare_data(hdf5_data, phot_data, phot_data_nod, LFs_dust, LFs_nodust, index, nbands, 
                 fdisk_emission, fbulge_m_emission, fbulge_d_emission):
   
    #star_formation_histories and SharkSED have the same number of galaxies in the same order, and so we can safely assume that to be the case.
    #to select the same galaxies in galaxies.hdf5 we need to ask for all of those that have a stellar mass > 0, and then assume that they are in the same order.

    (h0, _, mdisk, mbulge, mhalo, mshalo, typeg, age, 
     sfr_disk, sfr_burst, id_gal) = hdf5_data
   
    #(bulge_diskins_hist, bulge_mergers_hist, disk_hist) = sfh

    #components:
    #(len(my_data), 2, 2, 5, nbands)
    #0: disk instability bulge
    #1: galaxy merger bulge
    #2: total bulge
    #3: disk
    #4: total

    ind = np.where(mdisk + mbulge > 0)
    #print(len(mdisk[ind]))
    #print(phot_data)
    SEDs_dust = np.zeros(shape = (len(mdisk[ind]), 5, nbands))
    SEDs_nodust = np.zeros(shape = (len(mdisk[ind]), 5, nbands))

    p = 0
    for c in range(0,5):
        indust = phot_data[p]
        innodust = phot_data_nod[p]
        for i in range(0,nbands):
            SEDs_dust[:,c,i] = indust[i,:]
            SEDs_nodust[:,c,i] = innodust[i,:]
        p = p + 1

    for i in range(0,nbands):
        for c in range(0,5):
            #calculate LF with bands with dust
            ind = np.where(SEDs_dust[:,c,i] < -1)
            H, bins_edges = np.histogram(SEDs_dust[ind,c,i],bins=np.append(mbins,mupp))
            LFs_dust[index,c,i,:] = LFs_dust[index,c,i,:] + H

            #calculate LF of intrinsic bands 
            ind = np.where(SEDs_nodust[:,c,i] < -1)
            H, bins_edges = np.histogram(SEDs_nodust[ind,c,i],bins=np.append(mbins,mupp))
            LFs_nodust[index,c,i,:] = LFs_nodust[index,c,i,:] + H


    for i in range(0,nbands):

        for m in range(0,len(xlf)):
            ind = np.where((SEDs_dust[:,4,i] > xlf[m]-dm/2.0) & (SEDs_dust[:,4,i] <= xlf[m]+dm/2.0))
            fdisk_emission[index,i,m]    = np.sum(pow(10.0,SEDs_dust[ind,3,i]/(-2.5)))/np.sum(pow(10.0,SEDs_dust[ind,4,i]/(-2.5)))
            fbulge_m_emission[index,i,m] = np.sum(pow(10.0,SEDs_dust[ind,1,i]/(-2.5)))/np.sum(pow(10.0,SEDs_dust[ind,4,i]/(-2.5)))
            fbulge_d_emission[index,i,m] = np.sum(pow(10.0,SEDs_dust[ind,0,i]/(-2.5)))/np.sum(pow(10.0,SEDs_dust[ind,4,i]/(-2.5)))

def main(model_dir, outdir, redshift_table, subvols, obsdir):

    # Loop over redshift and subvolumes
    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mvir_hosthalo',
                           'mvir_subhalo', 'type', 'mean_stellar_age', 
                           'sfr_disk', 'sfr_burst', 'id_galaxy')}

    #sfh_fields = {'bulges_diskins': ('star_formation_rate_histories'),
    #              'bulges_mergers': ('star_formation_rate_histories'),
    #              'disks': ('star_formation_rate_histories')}

    Variable_Ext = True

    fields_sed = {'SED/ab_dust': ('bulge_d','bulge_m','bulge_t','disk','total'),}
    fields_sed_nod = {'SED/ab_nodust': ('bulge_d','bulge_m','bulge_t','disk','total')}

    z = (0, 0.25, 0.5) #, 1.0, 1.5, 2.0, 3.0)
    snapshots = redshift_table[z]

    file_hdf5_sed = "Shark-SED-eagle-rr14.hdf5" 
    # Create histogram
    for index, snapshot in enumerate(snapshots):

        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        if(Variable_Ext == False):
           seds = common.read_photometry_data(model_dir, snapshot, fields_sed, subvols)
           seds_nod = common.read_photometry_data(model_dir, snapshot, fields_sed_nod, subvols)
        else:
           seds = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed, subvols, file_hdf5_sed)
           seds_nod = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed_nod, subvols, file_hdf5_sed)

        nbands = len(seds[0]) 

        if(index == 0):
            LFs_dust     = np.zeros(shape = (len(z), 5, nbands, len(mbins)))
            LFs_nodust   = np.zeros(shape = (len(z), 5, nbands, len(mbins)))
            fdisk_emission    = np.zeros(shape = (len(z), nbands, len(mbins)))
            fbulge_m_emission = np.zeros(shape = (len(z), nbands, len(mbins)))
            fbulge_d_emission = np.zeros(shape = (len(z), nbands, len(mbins)))

        prepare_data(hdf5_data, seds, seds_nod, LFs_dust, LFs_nodust, index, nbands, 
                     fdisk_emission, fbulge_m_emission, fbulge_d_emission)

        h0, volh = hdf5_data[0], hdf5_data[1]
        if(volh > 0.):
            LFs_dust[index,:]   = LFs_dust[index,:]/volh
            LFs_nodust[index,:] = LFs_nodust[index,:]/volh

    # Take logs
    ind = np.where(LFs_dust > 0.)
    LFs_dust[ind] = np.log10(LFs_dust[ind])

    ind = np.where(LFs_nodust > 0.)
    LFs_nodust[ind] = np.log10(LFs_nodust[ind])

    if(Variable_Ext):
       outdir = os.path.join(outdir, 'eagle-rr14')
       #check if directory exists, otherwise create it
       if not os.path.exists(outdir):
          os.makedirs(outdir)


    plot_lfs(plt, outdir, obsdir, h0, LFs_dust, LFs_nodust)
    #plot_flux_contributions(plt, outdir, obsdir, h0, fdisk_emission, fbulge_m_emission, fbulge_d_emission)

if __name__ == '__main__':
    main(*common.parse_args())
