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

def plot_ir_lf_z0(plt, outdir, obsdir, h0, LFs_dust, LFs_nodust,  LFs_dust2, LFs_dust3, LFs_dust4):

    volcorr = 3.0*np.log10(h0)

    hcorr = 5.0*np.log10(h0)
    xlf_obs  = xlf - hcorr
 
    xtit="$\\rm mag-5log(h) (AB)$"
    ytit="$\\rm log_{10}(\Phi/dex^{-1} h^3 {\\rm Mpc}^{-3})$"

    xmin, xmax, ymin, ymax = -30, -16, -5, -1
    xleg = xmin + 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)

    fig = plt.figure(figsize=(6.5,12))

    subplots = (311, 312, 313)
    idxs = (0, 1, 2)
    cs  = (4, 3, 2)
    colors = ('Indigo','Teal','OrangeRed')

    bands = (21, 22, 25)
    bids  = (0, 1, 2)
    labels= ('P160', 'S250', 'S500')

    z = 0
    for band, ib in zip(bands, bids):
        ax = fig.add_subplot(subplots[ib])
        if (ib == 2):
            xtitplot = xtit
        else:
            xtitplot = ' '
        ytitplot = ytit

        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytitplot, locators=(2, 2, 1, 1))
        ax.text(-29, -1.5, labels[ib], fontsize=16, color=colors[0])

        if band == 21:
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

           ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'k', mec='k',marker='*')

           dml       = 0.4
           xm        = [25.0,25.0]
           xm[1]     = xm[1] + dml
           tenpctocm = 3.086e19
           corrpc    = np.log10(tenpctocm)*2.0
           xobs      = -2.5*(lm[15:21]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
           xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr + 5.0*corrhobs
           dmm       = np.abs(xm[0] - xm[1])
           yobs      = np.log10(pow(10.0,p[15:21]) * dml/dmm)-3.0*np.log10(h0)- 3.0 * corrhobs
           ydn       = np.log10(pow(10.0,dpn[15:21]) * dml/dmm)-3.0*np.log10(h0)- 3.0 * corrhobs
           yup       = np.log10(pow(10.0,dpu[15:21]) * dml/dmm)-3.0*np.log10(h0)- 3.0 * corrhobs
           ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'k', mec='k',marker='p')

        elif band == 22:
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
            yobs      = np.log10(p[0:7] * dml/dmm)-3.0*np.log10(h0)- 3.0 * corrhobs
            ydn       = np.log10((p[0:7]- dp[0:7]) * dml/dmm)-3.0*np.log10(h0)- 3.0 * corrhobs
            yup       = np.log10((p[0:7]+ dp[0:7]) * dml/dmm)-3.0*np.log10(h0)- 3.0 * corrhobs
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'k', mec='k',marker='o')

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

            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'k', mec='k',marker='*')

           
        elif band == 25:
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

            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'k', mec='k',marker='*')

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

            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'k', mec='k',marker='s')

            xobs = np.zeros(shape = 2)
            xobs[0] = 1.0
            xobs[1] = 1.0
            yobs=xobs
            ydn=xobs
            yup=xobs
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'k', mec='k',marker='p',label="Patel+2013")
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'k', mec='k',marker='o',label="Dye+2010")
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'k', mec='k',marker='s',label="Negrello+2013")
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'k', mec='k',marker='*',label="Marchetti+2016")

        c = 4
        idx = 0
        #Predicted LF
        ind = np.where(LFs_dust[z,c,band,:] < 0.)
        y = LFs_dust[z,c,band,ind]-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0], colors[idx], linewidth=3, alpha=0.5)#, label='$\\rm EAGLE-\\tau\, RR14-steep$')
        ind = np.where(LFs_dust2[z,c,band,:] < 0.)
        y = LFs_dust2[z,c,band,ind]-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0], colors[idx], linewidth=3, linestyle='dotted')#,label='$\\rm EAGLE-\\tau\, RR14$')
        ind = np.where(LFs_dust3[z,c,band,:] < 0.)
        y = LFs_dust3[z,c,band,ind]-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0], colors[idx], linewidth=3, alpha=0.4, linestyle='dashed')#, label='$\\rm EAGLE-\\tau\,f_{\\rm dust}\, const$')
        ind = np.where(LFs_dust4[z,c,band,:] < 0.)
        y = LFs_dust4[z,c,band,ind]-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0], colors[idx], linewidth=3, alpha=0.3, linestyle='dashdot')#, label='CF00')
    
        if (ib == 2):
           common.prepare_legend(ax, ['k','k','k','k'], loc='lower right', handlelength=5)

    common.savefig(outdir, fig, "IR_luminosity_function_z0_total.pdf")


def plot_uv_lf_z0(plt, outdir, obsdir, h0, LFs_dust, LFs_nodust,  LFs_dust2, LFs_dust3, LFs_dust4):

    volcorr = 3.0*np.log10(h0)

    hcorr = 5.0*np.log10(h0)
    xlf_obs  = xlf - hcorr
 
    xtit="$\\rm mag-5log(h) (AB)$"
    ytit="$\\rm log_{10}(\Phi/(0.5\\, {\\rm mag})/h^3 {\\rm Mpc}^{-3})$"

    xmin, xmax, ymin, ymax = -24, -13, -5, -1
    xleg = xmin + 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)

    fig = plt.figure(figsize=(8.5,12))

    subplots = (321, 323, 325, 322, 324, 326)
    idxs = (0, 1, 2)
    cs  = (4, 3, 2)
    colors = ('Indigo','Teal','OrangeRed')
    bands = (0, 1)

    labels= ('GALEX FUV', 'GALEX NUV')
    obs = ('lf1500', 'lf2300')

    labelsc = ('Total', 'Disk', 'Bulge')
    z = 0
    sp = 0
    for band in bands:

        file = obsdir+'/lf/'+obs[band]+'_z0_driver12.data'
        lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
        indx = np.where(p > 0)
        yobs = np.log10(p[indx])
        ydn  = np.log10(p[indx]-dp[indx])
        yup  = np.log10(p[indx]+dp[indx])

        for c, idx in zip(cs, idxs):
            ax = fig.add_subplot(subplots[sp])
            if (c == 2):
                xtitplot = xtit
            else:
                xtitplot = ' '
            if (band == 0):
                ytitplot = ytit
            else:
                ytitplot = ' '

            common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytitplot, locators=(2, 2, 1, 1))
            if(c == 4):
               ax.text(-23, -2, labels[band], fontsize=16, color=colors[idx])
            if(band == 0): 
               ax.text(-23, -1.5, labelsc[idx], fontsize=16, color=colors[idx])

            if(c == 4):
               ax.errorbar(lm[indx], yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'Black', mec='Black',marker='o')
            else:
               ind = np.where(LFs_dust[0,4,band,:] < 0.)
               y = LFs_dust[0,4,band,ind]
               ax.plot(xlf_obs[ind],y[0], 'Indigo', linewidth=3, alpha=0.5)

            #Predicted LF
            ind = np.where(LFs_nodust[z,c,band,:] < 0.)
            y = LFs_nodust[z,c,band,ind]
            ax.plot(xlf_obs[ind],y[0],colors[idx], linewidth=1, label='intrinsic')
    
            ind = np.where(LFs_dust[z,c,band,:] < 0.)
            y = LFs_dust[z,c,band,ind]
            ax.plot(xlf_obs[ind],y[0], colors[idx], linewidth=3, alpha=0.5, label='$\\rm EAGLE-\\tau\, RR14-steep$')
            ind = np.where(LFs_dust2[z,c,band,:] < 0.)
            y = LFs_dust2[z,c,band,ind]
            ax.plot(xlf_obs[ind],y[0], colors[idx], linewidth=3, linestyle='dotted',label='$\\rm EAGLE-\\tau\, RR14$')
            ind = np.where(LFs_dust3[z,c,band,:] < 0.)
            y = LFs_dust3[z,c,band,ind]
            ax.plot(xlf_obs[ind],y[0], colors[idx], linewidth=3, alpha=0.4, linestyle='dashed', label='$\\rm EAGLE-\\tau\,f_{\\rm dust}\, const$')
            ind = np.where(LFs_dust4[z,c,band,:] < 0.)
            y = LFs_dust4[z,c,band,ind]
            ax.plot(xlf_obs[ind],y[0], colors[idx], linewidth=3, alpha=0.3, linestyle='dashdot', label='CF00')
    
            if (sp == 0):
                common.prepare_legend(ax, [colors[idx], colors[idx], colors[idx], colors[idx], colors[idx]], bbox_to_anchor=[0.6, 1], handlelength=5)

            sp = sp + 1

    common.savefig(outdir, fig, "UV_luminosity_function_z0_total.pdf")


def plot_uv_lf_evo(plt, outdir, obsdir, h0, LFs_dust, LFs_nodust,  LFs_dust2, LFs_dust3, LFs_dust4, nbands):

    volcorr = 3.0*np.log10(h0)
    xlf_obs  = xlf
 
    xtit="$\\rm 1500\\AA\, mag\, (AB)$ top-hat $100\\AA$"
    ytit="$\\rm log_{10}(\Phi/{\\rm dex^{-1}} {\\rm Mpc}^{-3})$"

    xmin, xmax, ymin, ymax = -25, -15, -6, -1
    xleg = xmin + 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)

    fig = plt.figure(figsize=(5,16))

    subplots = (511, 512, 513, 514, 515)
    idxs = (0, 1, 2, 3, 4)
    zs  = (1, 2, 3, 4, 5)
    band = nbands-1
    labels= ('z=3', 'z=4', 'z=6', 'z=8', 'z=10')
  
    corrm_obs = -5.0*np.log10(h0/0.7) 
    corry_obs = 3.0*np.log10(h0/0.7)
    for subplot, idx, z in zip(subplots, idxs, zs):

        ax = fig.add_subplot(subplot)
        ytitplot = ytit
        if (idx == 4):
            xtitplot = xtit
        else:
            xtitplot = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytitplot, locators=(2, 2, 1, 1))
        ax.text(xleg,yleg, labels[idx])

        if(idx == 0):
           ax.text(-24.5, -2.5, 'Total', fontsize=16, color='Indigo')
           file = obsdir+'/lf/lf1700_z3_sawicki06.data'
           lm,p,dp = np.loadtxt(file,usecols=[0,2,3],unpack=True)
           indx = np.where(p > 0)
           yobs = np.log10(p[indx])
           ydn  = np.log10(p[indx]-dp[indx])
           yup  = np.log10(p[indx]+dp[indx])
           ax.errorbar(lm[indx]+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o') #,label="Sawicki+2006")

           file = obsdir+'/lf/lf1700_z3_reddy09.data'
           lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
           indx = np.where(p > 0)
           yobs = np.log10(p[indx])
           ydn  = np.log10(p[indx]-dp[indx])
           yup  = np.log10(p[indx]+dp[indx])
           ax.errorbar(lm[indx]+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s') #,label="Reddy+2009")

        if(idx == 1):
           file = obsdir+'/lf/lf1500_z4_adams19.data'
           lmA19z4,Ap4,Adp4 = np.loadtxt(file,usecols=[0, 1, 2],unpack=True)
           yobs = np.log10(Ap4*1e-4)
           ydn  = np.log10(Ap4*1e-4 - Adp4*1e-4)
           yup  = np.log10(Ap4*1e-4 + Adp4*1e-4)
           ax.errorbar(lmA19z4+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='d',label="Adams+2019")

           file = obsdir+'/lf/lf1600_z4-10_Bouwens2015.data'
           lmB15z4,p4,dp4,lmB15z6,p6,dp6,lmB15z8,p8,dp8 = np.loadtxt(file,usecols=[0, 1, 2, 6, 7, 8, 12, 13, 14],unpack=True)
           yobs = np.log10(p4)
           ydn  = np.log10(p4-dp4)
           yup  = np.log10(p4+dp4)
           ax.errorbar(lmB15z4+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='^',label="Bouwens+2015")

           file = obsdir+'/lf/lf1500_z4-8_Finkelstein2015.data'
           lmF15,pF4,dpuF4,dpdF4,pF6,dpuF6,dpdF6,pF8,dpuF8,dpdF8 = np.loadtxt(file,usecols=[0,1, 2, 3, 7, 8, 9, 13, 14, 15],unpack=True)
           yobs = np.log10(pF4*1e-3)
           ydn  = np.log10(pF4*1e-3-dpdF4*1e-3)
           yup  = np.log10(pF4*1e-3+dpuF4*1e-3)
           ax.errorbar(lmF15+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='v',label="Finkelstein+2015")

        if(idx == 2):
           yobs = np.log10(p6)
           ydn  = np.log10(p6-dp6)
           yup  = np.log10(p6+dp6)
           ax.errorbar(lmB15z6+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='^')

           yobs = np.log10(pF6*1e-3)
           ydn  = np.log10(pF6*1e-3-dpdF6*1e-3)
           yup  = np.log10(pF6*1e-3+dpuF6*1e-3)
           ax.errorbar(lmF15+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='v')

        if(idx == 3):
           yobs = np.log10(p8)
           ydn  = np.log10(p8-dp8)
           yup  = np.log10(p8+dp8)
           ax.errorbar(lmB15z8+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='^')

           yobs = np.log10(pF8*1e-3)
           ydn  = np.log10(pF8*1e-3-dpdF8*1e-3)
           yup  = np.log10(pF8*1e-3+dpuF8*1e-3)
           ax.errorbar(lmF15+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='v')

        if(idx == 4):
           file = obsdir+'/lf/lf1500_z10_oesch2018.data'
           lm,p,dpu,dpd = np.loadtxt(file,usecols=[0,1, 2, 3],unpack=True)
           ax.errorbar(lm+corrm_obs, p+corry_obs, yerr=[p-dpu,dpd-p], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*', label='Oesch+2018')

        #Predicted LF
        ind = np.where(LFs_nodust[z,4,band,:] < 0.)
        y = LFs_nodust[z,4,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'Indigo', linewidth=1, label='intrinsic')

        ind = np.where(LFs_dust[z,4,band,:] < 0.)
        y = LFs_dust[z,4,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'Indigo', linewidth=3, alpha=0.5, label='$\\rm EAGLE-\\tau\, RR14-steep$')
        ind = np.where(LFs_dust2[z,4,band,:] < 0.)
        y = LFs_dust2[z,4,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'Indigo', linewidth=3, linestyle='dotted',label='$\\rm EAGLE-\\tau\, RR14$')
        ind = np.where(LFs_dust3[z,4,band,:] < 0.)
        y = LFs_dust3[z,4,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'Indigo', linewidth=3, alpha=0.4, linestyle='dashed', label='$\\rm EAGLE-\\tau\,f_{\\rm dust}\, const$')
        ind = np.where(LFs_dust4[z,4,band,:] < 0.)
        y = LFs_dust4[z,4,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'Indigo', linewidth=3, alpha=0.3, linestyle='dashdot', label='CF00')

        if (idx == 0):
            common.prepare_legend(ax, ['Indigo','Indigo','Indigo','Indigo','Indigo'], bbox_to_anchor=[0.1, 1], handlelength=5)
        #ind = np.where(LFs_dust[z,3,band,:] < 0.)
        #y = LFs_dust[z,3,band,ind]+volcorr-np.log10(dm)
        #ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
        #ind = np.where(LFs_dust[z,2,band,:] < 0.)
        #y = LFs_dust[z,2,band,ind]+volcorr-np.log10(dm)
        #ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')

    common.savefig(outdir, fig, "UV_luminosity_function_evolution_total.pdf")

    #now plot disks only
    fig = plt.figure(figsize=(5,16))
    for subplot, idx, z in zip(subplots, idxs, zs):

        ax = fig.add_subplot(subplot)
        ytitplot = ytit
        if (idx == 4):
            xtitplot = xtit
        else:
            xtitplot = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytitplot, locators=(2, 2, 1, 1))
        ax.text(xleg,yleg, labels[idx])

        if(idx == 0):
           ax.text(-24.5, -2.5, 'Disk', fontsize=16, color='Teal')

        ind = np.where(LFs_dust[z,4,band,:] < 0.)
        y = LFs_dust[z,4,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'Indigo', linewidth=3, alpha=0.5)

        #Predicted LF
        ind = np.where(LFs_nodust[z,3,band,:] < 0.)
        y = LFs_nodust[z,3,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'Teal', linewidth=1, label='intrinsic')

        ind = np.where(LFs_dust[z,3,band,:] < 0.)
        y = LFs_dust[z,3,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'Teal', linewidth=3, alpha=0.5, label='$\\rm EAGLE-\\tau\, RR14-steep$')
        ind = np.where(LFs_dust2[z,3,band,:] < 0.)
        y = LFs_dust2[z,3,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'Teal', linewidth=3, linestyle='dotted',label='$\\rm EAGLE-\\tau\, RR14$')
        ind = np.where(LFs_dust3[z,3,band,:] < 0.)
        y = LFs_dust3[z,3,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'Teal', linewidth=3, alpha=0.4, linestyle='dashed', label='$\\rm EAGLE-\\tau\,f_{\\rm dust}\, const$')
        ind = np.where(LFs_dust4[z,3,band,:] < 0.)
        y = LFs_dust4[z,3,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'Teal', linewidth=3, alpha=0.3, linestyle='dashdot', label='CF00')

        if (idx == 0):
            common.prepare_legend(ax, ['Teal','Teal','Teal','Teal','Teal'], bbox_to_anchor=[0.1, 1], handlelength=5)
        #ind = np.where(LFs_dust[z,3,band,:] < 0.)
        #y = LFs_dust[z,3,band,ind]+volcorr-np.log10(dm)
        #ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
        #ind = np.where(LFs_dust[z,2,band,:] < 0.)
        #y = LFs_dust[z,2,band,ind]+volcorr-np.log10(dm)
        #ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')

    common.savefig(outdir, fig, "UV_luminosity_function_evolution_disk.pdf")

    #now plot bulges only
    fig = plt.figure(figsize=(5,16))
    for subplot, idx, z in zip(subplots, idxs, zs):

        ax = fig.add_subplot(subplot)
        ytitplot = ytit
        if (idx == 4):
            xtitplot = xtit
        else:
            xtitplot = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytitplot, locators=(2, 2, 1, 1))
        ax.text(xleg,yleg, labels[idx])

        if(idx == 0):
           ax.text(-24.5, -2.5, 'Bulge', fontsize=16, color='OrangeRed')
        ind = np.where(LFs_dust[z,4,band,:] < 0.)
        y = LFs_dust[z,4,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'Indigo', linewidth=3, alpha=0.5)


        #Predicted LF
        ind = np.where(LFs_nodust[z,2,band,:] < 0.)
        y = LFs_nodust[z,2,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'OrangeRed', linewidth=1, label='intrinsic')

        ind = np.where(LFs_dust[z,2,band,:] < 0.)
        y = LFs_dust[z,2,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'OrangeRed', linewidth=3, alpha=0.5, label='$\\rm EAGLE-\\tau\, RR14-steep$')
        ind = np.where(LFs_dust2[z,2,band,:] < 0.)
        y = LFs_dust2[z,2,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'OrangeRed', linewidth=3,  linestyle='dotted',label='$\\rm EAGLE-\\tau\, RR14$')
        ind = np.where(LFs_dust3[z,2,band,:] < 0.)
        y = LFs_dust3[z,2,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'OrangeRed', linewidth=3, alpha=0.4, linestyle='dashed', label='$\\rm EAGLE-\\tau\,f_{\\rm dust}\, const$')
        ind = np.where(LFs_dust4[z,2,band,:] < 0.)
        y = LFs_dust4[z,2,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'OrangeRed', linewidth=3, alpha=0.3, linestyle='dashdot', label='CF00')

        if (idx == 0):
            common.prepare_legend(ax, ['OrangeRed','OrangeRed','OrangeRed','OrangeRed','OrangeRed'], bbox_to_anchor=[0.1, 1], handlelength=5)
        #ind = np.where(LFs_dust[z,3,band,:] < 0.)
        #y = LFs_dust[z,3,band,ind]+volcorr-np.log10(dm)
        #ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
        #ind = np.where(LFs_dust[z,2,band,:] < 0.)
        #y = LFs_dust[z,2,band,ind]+volcorr-np.log10(dm)
        #ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')

    common.savefig(outdir, fig, "UV_luminosity_function_evolution_bulge.pdf")

def prepare_data(hdf5_data, phot_data, phot_data_nod, seds2, seds3, seds4, LFs_dust, LFs_nodust, LFs_dust2, LFs_dust3, LFs_dust4, index, nbands):
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
    SEDs_dust = np.zeros(shape = (len(mdisk[ind]), 5, nbands))
    SEDs_nodust = np.zeros(shape = (len(mdisk[ind]), 5, nbands))
    SEDs_dust2 = np.zeros(shape = (len(mdisk[ind]), 5, nbands))
    SEDs_dust3 = np.zeros(shape = (len(mdisk[ind]), 5, nbands))
    SEDs_dust4 = np.zeros(shape = (len(mdisk[ind]), 5, nbands))

    p = 0
    for c in range(0,5):
        indust = phot_data[p]
        innodust = phot_data_nod[p]
        indust2 = seds2[p]
        indust3 = seds3[p]
        indust4 = seds4[p]
        for i in range(0,nbands):
            SEDs_dust[:,c,i] = indust[i,:]
            SEDs_nodust[:,c,i] = innodust[i,:]
            SEDs_dust2[:,c,i] = indust2[i,:]
            SEDs_dust3[:,c,i] = indust3[i,:]
            SEDs_dust4[:,c,i] = indust4[i,:]
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

            #calculate LF for other models
            ind = np.where(SEDs_dust2[:,c,i] < -1)
            H, bins_edges = np.histogram(SEDs_dust2[ind,c,i],bins=np.append(mbins,mupp))
            LFs_dust2[index,c,i,:] = LFs_dust2[index,c,i,:] + H
            ind = np.where(SEDs_dust3[:,c,i] < -1)
            H, bins_edges = np.histogram(SEDs_dust3[ind,c,i],bins=np.append(mbins,mupp))
            LFs_dust3[index,c,i,:] = LFs_dust3[index,c,i,:] + H
            ind = np.where(SEDs_dust4[:,c,i] < -1)
            H, bins_edges = np.histogram(SEDs_dust4[ind,c,i],bins=np.append(mbins,mupp))
            LFs_dust4[index,c,i,:] = LFs_dust4[index,c,i,:] + H

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

    z = (0.0, 3.0, 4.0, 6.0, 8.0, 10.0) #, 1.0, 1.5, 2.0)
    snapshots = redshift_table[z]

    file_hdf5_sed = "Shark-SED-eagle-rr14-steep.hdf5"
    file_hdf5_sed2 = "Shark-SED-eagle-rr14.hdf5"
    file_hdf5_sed3 = "Shark-SED-eagle-const.hdf5"
    file_hdf5_sed4 = "Shark-SED.hdf5"

    # Create histogram
    for index, snapshot in enumerate(snapshots):

        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        if(Variable_Ext == False):
           seds = common.read_photometry_data(model_dir, snapshot, fields_sed, subvols)
           seds_nod = common.read_photometry_data(model_dir, snapshot, fields_sed_nod, subvols)
        else:
           seds = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed, subvols, file_hdf5_sed)
           seds_nod = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed_nod, subvols, file_hdf5_sed)
           seds2 = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed, subvols, file_hdf5_sed2)
           seds3 = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed, subvols, file_hdf5_sed3)
           seds4 = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed, subvols, file_hdf5_sed4)

        nbands = len(seds[0]) 

        if(index == 0):
            LFs_dust     = np.zeros(shape = (len(z), 5, nbands, len(mbins)))
            LFs_nodust   = np.zeros(shape = (len(z), 5, nbands, len(mbins)))
            LFs_dust2     = np.zeros(shape = (len(z), 5, nbands, len(mbins)))
            LFs_dust3     = np.zeros(shape = (len(z), 5, nbands, len(mbins)))
            LFs_dust4     = np.zeros(shape = (len(z), 5, nbands, len(mbins)))

        prepare_data(hdf5_data, seds, seds_nod, seds2, seds3, seds4, LFs_dust, LFs_nodust, LFs_dust2, LFs_dust3, LFs_dust4, index, nbands)

        h0, volh = hdf5_data[0], hdf5_data[1]
        if(volh > 0.):
            LFs_dust[index,:]    = LFs_dust[index,:]/volh
            LFs_nodust[index,:]  = LFs_nodust[index,:]/volh
            LFs_dust2[index,:]   = LFs_dust2[index,:]/volh
            LFs_dust3[index,:]   = LFs_dust3[index,:]/volh
            LFs_dust4[index,:]   = LFs_dust4[index,:]/volh

    # Take logs
    ind = np.where(LFs_dust > 0.)
    LFs_dust[ind] = np.log10(LFs_dust[ind])

    ind = np.where(LFs_nodust > 0.)
    LFs_nodust[ind] = np.log10(LFs_nodust[ind])

    ind = np.where(LFs_dust2 > 0.)
    LFs_dust2[ind] = np.log10(LFs_dust2[ind])

    ind = np.where(LFs_dust3 > 0.)
    LFs_dust3[ind] = np.log10(LFs_dust3[ind])

    ind = np.where(LFs_dust4 > 0.)
    LFs_dust4[ind] = np.log10(LFs_dust4[ind])

    if(Variable_Ext):
       outdir = os.path.join(outdir, 'eagle-rr14')

    plot_uv_lf_evo(plt, outdir, obsdir, h0, LFs_dust, LFs_nodust, LFs_dust2, LFs_dust3, LFs_dust4, nbands)
    plot_uv_lf_z0(plt, outdir, obsdir, h0, LFs_dust, LFs_nodust, LFs_dust2, LFs_dust3, LFs_dust4)
    plot_ir_lf_z0(plt, outdir, obsdir, h0, LFs_dust, LFs_nodust, LFs_dust2, LFs_dust3, LFs_dust4)

if __name__ == '__main__':
    main(*common.parse_args())
