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

import functools
import numpy as np
import os

import common
import utilities_statistics as us

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


mflow = 4
mfupp = 14
dmf = 0.5
mfbins = np.arange(mflow,mfupp,dmf)
xmf = mfbins + dmf/2.0

mflow2 = 4
mfupp2 = 14
dmf2 = 0.25
mfbins2 = np.arange(mflow2,mfupp2,dmf2)
xmf2 = mfbins2 + dmf2/2.0


def plot_uv_lf_evo(plt, outdir, obsdir, h0, LFs_dust, LFs_nodust, nbands):

    volcorr = 3.0*np.log10(h0)
    xlf_obs  = xlf
 
    xtit="$\\rm 1500\\AA\, mag\, (AB)$"
    ytit="$\\rm log_{10}(\Phi/{\\rm dex^{-1}} {\\rm Mpc}^{-3})$"

    xmin, xmax, ymin, ymax = -25, -15, -7, -1
    xleg = xmin + 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)

    fig = plt.figure(figsize=(5,14))

    subplots = (711, 712, 713, 714, 715, 716, 717)
    idx = (0, 1, 2, 3, 4, 5, 6)
    zs  = (0, 1, 2, 3, 4, 5, 6)
    band = 0 #28
    labels= ('z=3', 'z=4', 'z=6', 'z=8', 'z=10', 'z=13', 'z=17')
  
    corrm_obs = -5.0*np.log10(h0/0.7) 
    corry_obs = 3.0*np.log10(h0/0.7)
    for subplot, idx, z in zip(subplots, idx, zs):

        ax = fig.add_subplot(subplot)
        ytitplot = ytit
        if (idx == 6):
            xtitplot = xtit
        else:
            xtitplot = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytitplot, locators=(2, 2, 1, 1))
        ax.text(xleg,yleg, labels[idx])

        if(idx == 0):
           file = obsdir+'/lf/lf1700_z3_sawicki06.data'
           lm,p,dp = np.loadtxt(file,usecols=[0,2,3],unpack=True)
           indx = np.where(p > 0)
           yobs = np.log10(p[indx])
           ydn  = np.log10(p[indx]-dp[indx])
           yup  = np.log10(p[indx]+dp[indx])
           ax.errorbar(lm[indx]+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Sawicki+2006")

           file = obsdir+'/lf/lf1700_z3_reddy09.data'
           lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
           indx = np.where(p > 0)
           yobs = np.log10(p[indx])
           ydn  = np.log10(p[indx]-dp[indx])
           yup  = np.log10(p[indx]+dp[indx])
           ax.errorbar(lm[indx]+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s',label="Reddy+2009")

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

           file = obsdir+'/lf/lf1500_z8_adams23.data'
           lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
           ax.errorbar(lm+corrm_obs, np.log10(p * 1e-5)+corry_obs, yerr=[np.log10(p)-np.log10(dp),np.log10(p + dp) - np.log10(p)], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D', label='Adams+2023')

        if(idx == 4):
           file = obsdir+'/lf/lf1500_z10_oesch2018.data'
           lm,p,dpu,dpd = np.loadtxt(file,usecols=[0,1, 2, 3],unpack=True)
           ax.errorbar(lm+corrm_obs, p+corry_obs, yerr=[p-dpu,dpd-p], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*', label='Oesch+2018')

           file = obsdir+'/lf/lf1500_z10_adams23.data'
           lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
           ax.errorbar(lm+corrm_obs, np.log10(p * 1e-5)+corry_obs, yerr=[np.log10(p)-np.log10(dp),np.log10(p + dp) - np.log10(p)], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D')

           file = obsdir+'/lf/lf1500_z10_weibel25.data'
           lm,p,dpup,dpdn,flag = np.loadtxt(file,usecols=[0,1,2,3,4],unpack=True)
           ax.errorbar(lm+corrm_obs, p+corry_obs, yerr=[dpup, dpdn], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o', label='Weibel+25')

           file = obsdir+'/lf/lf1500_z10_whitler25.data'
           lm,p,dpdn,dpup = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
           ax.errorbar(lm+corrm_obs, np.log10(p)+corry_obs, yerr=[np.log10(p) - np.log10(p-dpdn), np.log10(p + dpup) - np.log10(p)], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s', label='Whitler+25')

        if(idx == 5):
           file = obsdir+'/lf/lf1500_z13_weibel25.data'
           lm,p,dpup,dpdn,flag = np.loadtxt(file,usecols=[0,1,2,3,4],unpack=True)
           ax.errorbar(lm+corrm_obs, p+corry_obs, yerr=[dpup, dpdn], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

           file = obsdir+'/lf/lf1500_z12p8_whitler25.data'
           lm,p,dpdn,dpup = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
           ax.errorbar(lm+corrm_obs, np.log10(p)+corry_obs, yerr=[np.log10(p) - np.log10(p-dpdn), np.log10(p + dpup) - np.log10(p)], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s')

        if(idx == 6):
           file = obsdir+'/lf/lf1500_z17_weibel25.data'
           lm,p= np.loadtxt(file,usecols=[0,1],unpack=True)
           ax.plot(lm+corrm_obs, p+corry_obs, ls='None',marker='o', fillstyle ='none', color='grey')
           for x,y in zip(lm+corrm_obs, p+corry_obs):
               ax.arrow(x,y,0,-0.75, head_width=0.05, head_length=0.1, color='grey')

           file = obsdir+'/lf/lf1500_z15_whitler25.data'
           lm,p,dpdn,dpup = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
           ax.errorbar(lm+corrm_obs, np.log10(p)+corry_obs, yerr=[np.log10(p) - np.log10(p-dpdn), np.log10(p + dpup) - np.log10(p)], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s')

        #Predicted LF
        ind = np.where(LFs_dust[z,4,band,:] < 0.)
        y = LFs_dust[z,4,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
        ind = np.where(LFs_nodust[z,4,band,:] < 0.)
        y = LFs_nodust[z,4,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)
        if(idx == 1):
           for a,b,c in zip(xlf_obs,LFs_dust[z,4,band,:],LFs_nodust[z,4,band,:]):
               print (a, b+volcorr-np.log10(dm),c+volcorr-np.log10(dm))

        ind = np.where(LFs_dust[z,3,band,:] < 0.)
        y = LFs_dust[z,3,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
        ind = np.where(LFs_dust[z,2,band,:] < 0.)
        y = LFs_dust[z,2,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
        if ((idx == 0) or (idx == 1) or (idx ==4) or (idx == 3)):
            common.prepare_legend(ax, ['grey','grey','grey'], loc=4)

    plt.tight_layout()
    common.savefig(outdir, fig, "UV_luminosity_function_evolution.pdf")


    fig = plt.figure(figsize=(12,8))
    xmin, xmax, ymin, ymax = -25, -13, -7, -1

    subplots = (231, 232, 233, 234, 235)
    idx = (2, 3, 4, 5, 6)
    zs  = (2, 3, 4, 5, 6)
    band = 0 #28

    for subplot, idx, z in zip(subplots, idx, zs):

        ax = fig.add_subplot(subplot)
        if (idx > 3):
            xtitplot = xtit
        else:
            xtitplot = ' '
        if ((idx == 2) | (idx == 5)):
            ytitplot = ytit
        else:
            ytitplot = ' '

        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytitplot, locators=(2, 2, 1, 1))
        ax.text(xleg,yleg, labels[idx])

        if(idx == 2):
           file = obsdir+'/lf/lf1600_z4-10_Bouwens2015.data'
           lmB15z4,p4,dp4,lmB15z6,p6,dp6,lmB15z8,p8,dp8 = np.loadtxt(file,usecols=[0, 1, 2, 6, 7, 8, 12, 13, 14],unpack=True)

           file = obsdir+'/lf/lf1500_z4-8_Finkelstein2015.data'
           lmF15,pF4,dpuF4,dpdF4,pF6,dpuF6,dpdF6,pF8,dpuF8,dpdF8 = np.loadtxt(file,usecols=[0,1, 2, 3, 7, 8, 9, 13, 14, 15],unpack=True)

           yobs = np.log10(p6)
           ydn  = np.log10(p6-dp6)
           yup  = np.log10(p6+dp6)
           ax.errorbar(lmB15z6+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='^', label="Bouwens+2015")

           yobs = np.log10(pF6*1e-3)
           ydn  = np.log10(pF6*1e-3-dpdF6*1e-3)
           yup  = np.log10(pF6*1e-3+dpuF6*1e-3)
           ax.errorbar(lmF15+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='v', label="Finkelstein+2015")

        if(idx == 3):
           yobs = np.log10(p8)
           ydn  = np.log10(p8-dp8)
           yup  = np.log10(p8+dp8)
           ax.errorbar(lmB15z8+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='^')

           yobs = np.log10(pF8*1e-3)
           ydn  = np.log10(pF8*1e-3-dpdF8*1e-3)
           yup  = np.log10(pF8*1e-3+dpuF8*1e-3)
           ax.errorbar(lmF15+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='v')

           file = obsdir+'/lf/lf1500_z8_adams23.data'
           lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
           ax.errorbar(lm+corrm_obs, np.log10(p * 1e-5)+corry_obs, yerr=[np.log10(p)-np.log10(dp),np.log10(p + dp) - np.log10(p)], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D', label='Adams+2023')

        if(idx == 4):
           file = obsdir+'/lf/lf1500_z10_oesch2018.data'
           lm,p,dpu,dpd = np.loadtxt(file,usecols=[0,1, 2, 3],unpack=True)
           ax.errorbar(lm+corrm_obs, p+corry_obs, yerr=[p-dpu,dpd-p], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*', label='Oesch+2018')

           file = obsdir+'/lf/lf1500_z10_adams23.data'
           lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
           ax.errorbar(lm+corrm_obs, np.log10(p * 1e-5)+corry_obs, yerr=[np.log10(p)-np.log10(dp),np.log10(p + dp) - np.log10(p)], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D')

           file = obsdir+'/lf/lf1500_z10_weibel25.data'
           lm,p,dpup,dpdn,flag = np.loadtxt(file,usecols=[0,1,2,3,4],unpack=True)
           ax.errorbar(lm+corrm_obs, p+corry_obs, yerr=[dpup, dpdn], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o', label='Weibel+2025')
           
           file = obsdir+'/lf/lf1500_z10_whitler25.data'
           lm,p,dpdn,dpup = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
           #ax.errorbar(lm+corrm_obs, np.log10(p)+corry_obs, yerr=[np.log10(p) - np.log10(p-dpdn), np.log10(p + dpup) - np.log10(p)], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s', label='Whitler+25')

        if(idx == 5):
           file = obsdir+'/lf/lf1500_z13_weibel25.data'
           lm,p,dpup,dpdn,flag = np.loadtxt(file,usecols=[0,1,2,3,4],unpack=True)
           ax.errorbar(lm+corrm_obs, p+corry_obs, yerr=[dpup, dpdn], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')
           file = obsdir+'/lf/lf1500_z12p8_whitler25.data'
           lm,p,dpdn,dpup = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
           #ax.errorbar(lm+corrm_obs, np.log10(p)+corry_obs, yerr=[np.log10(p) - np.log10(p-dpdn), np.log10(p + dpup) - np.log10(p)], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s')

        if(idx == 6):
           file = obsdir+'/lf/lf1500_z17_weibel25.data'
           lm,p= np.loadtxt(file,usecols=[0,1],unpack=True)
           ax.plot(lm+corrm_obs, p+corry_obs, ls='None',marker='o', fillstyle ='none', color='grey')
           for x,y in zip(lm+corrm_obs, p+corry_obs):
               ax.arrow(x,y,0,-0.75, head_width=0.05, head_length=0.1, color='grey')
           file = obsdir+'/lf/lf1500_z17_perez-gonzalez25.data'
           lm,p,dn,du,flag= np.loadtxt(file,usecols=[0,1,2,3,4],unpack=True)
           flagin = np.where(flag == 0)
           ax.errorbar(lm[flagin]+corrm_obs+0.15, p[flagin]+corry_obs, yerr=[dn[flagin], du[flagin]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s', label = 'Perez-Gonzalez+2025')

           ax.plot(lm[0]+corrm_obs, p[0]+corry_obs, ls='None',marker='s', fillstyle ='none', color='grey')
           ax.arrow(lm[0]+corrm_obs, p[0]+corry_obs,0,-0.75, head_width=0.05, head_length=0.1, color='grey')

           file = obsdir+'/lf/lf1500_z15_whitler25.data'
           lm,p,dpdn,dpup = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
           #ax.errorbar(lm+corrm_obs, np.log10(p)+corry_obs, yerr=[np.log10(p) - np.log10(p-dpdn), np.log10(p + dpup) - np.log10(p)], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s')

        #Predicted LF
        ind = np.where(LFs_dust[z,4,band,:] < 0.)
        y = LFs_dust[z,4,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'k', linewidth=3, label='Total attenuated' if idx == 2 else None)
        ind = np.where(LFs_nodust[z,4,band,:] < 0.)
        y = LFs_nodust[z,4,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'k', linewidth=1, label = 'Total unattenuated' if idx == 2 else None)
        #if(idx == 1):
        #   for a,b,c in zip(xlf_obs,LFs_dust[z,4,band,:],LFs_nodust[z,4,band,:]):
        #       print (a, b+volcorr-np.log10(dm),c+volcorr-np.log10(dm))

        ind = np.where(LFs_dust[z,3,band,:] < 0.)
        y = LFs_dust[z,3,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'LightSeaGreen', linewidth=2, linestyle='dotted', label='star formation disks' if idx == 2 else None)
        ind = np.where(LFs_dust[z,0,band,:] < 0.)
        y = LFs_dust[z,0,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'OrangeRed', linewidth=2, linestyle='dashed', label='starbursts (disk ins.)' if idx == 2 else None)
        ind = np.where(LFs_dust[z,1,band,:] < 0.)
        y = LFs_dust[z,1,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'Orange', linewidth=2, linestyle='dashdot', label='starbursts (mergers)' if idx == 2 else None)

        if(idx == 2):
            common.prepare_legend(ax, ['k','k','LightSeaGreen','OrangeRed','Orange','grey','grey','grey'],  bbox_to_anchor=(2.305, -0.9))
        if idx == 3:
            common.prepare_legend(ax, ['grey','grey','grey'], bbox_to_anchor=(1.11, -0.98))
        if idx == 4:
            common.prepare_legend(ax, ['grey','grey','grey'], bbox_to_anchor=(-0.09, -1.15))
        if idx == 6:
            common.prepare_legend(ax, ['grey','grey','grey'], bbox_to_anchor=(1.11, -0.025))

        #if ((idx ==4) or (idx == 3)):
        #    common.prepare_legend(ax, ['grey','grey','grey'], loc=0)

    #plt.tight_layout()
    common.savefig(outdir, fig, "UV_luminosity_function_evolution_zGT6.pdf")

def plot_uv_mass_evo(plt, outdir, obsdir, uvmhalo, uvmstar, mstarhalo, mstarsfr):

    fig = plt.figure(figsize=(5,7))
    xmin, xmax, ymin, ymax = 5, 12, -12, -24

    ytit="$\\rm1500\\AA\, mag\, (AB)$"
    xtit="$\\rm log_{10}(M_{\\star}/M_{\\odot})$"

    idxs = (2, 3, 4, 5, 6)
    zs  = (2, 3, 4, 5, 6)
    labels= ('z=3', 'z=4', 'z=6', 'z=8', 'z=10', 'z=13', 'z=17')

    ax = fig.add_subplot(211)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))

    cols = ('Indigo','purple','Navy','Aquamarine', 'Green','Gold','red','DarkRed') 

    for idx, z in zip(idxs, zs):
  
        #Predicted LF
        ind = np.where(uvmstar[z,0,:] != 0)
        x = xmf[ind]
        y = uvmstar[z,0,ind]
        ydn = uvmstar[z,1,ind]
        yup = uvmstar[z,2,ind]
        #ax.fill_between(x, ydn[0], yup[0], color=cols[z], alpha = 0.25, linestyle='solid', linewidth = 2)
        ax.plot(x, y[0],color=cols[z], linestyle='solid', linewidth = 2, label=labels[z])
    common.prepare_legend(ax, cols[2:7], loc=0)

    ax = fig.add_subplot(212)

    ytit="$\\rm \\sigma(1500\\AA\, mag)$"

    xmin, xmax, ymin, ymax = 5, 12, -2,2
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))

    for idx, z in zip(idxs, zs):
  
        #Predicted LF
        ind = np.where(uvmstar[z,0,:] != 0)
        x = xmf[ind]
        y = uvmstar[z,0,ind]
        ydn = uvmstar[z,1,ind] 
        yup = uvmstar[z,2,ind] 
        ax.plot(x, -1 * ydn[0],color=cols[z], linestyle='solid', linewidth = 2)
        ax.plot(x, yup[0],color=cols[z], linestyle='solid', linewidth = 2)
    ax.plot([xmin, xmax], [0,0], linestyle='dotted', color='k')

    plt.tight_layout()
    common.savefig(outdir, fig, "UV_StellarMass_zGT6.pdf")
   
    fig = plt.figure(figsize=(5,7))
    xmin, xmax, ymin, ymax = 9.8, 13, -12, -24

    ytit="$\\rm 1500\\AA\, mag\, (AB)$"
    xtit="$\\rm log_{10}(M_{\\rm halo}/M_{\\odot})$"
    #idxs = (2, 3, 4, 5)
    #zs  = (2, 3, 4, 5)

    ax = fig.add_subplot(211)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))

    for idx, z in zip(idxs, zs):
  
        #Predicted LF
        ind = np.where(uvmhalo[z,0,:] != 0)
        x = xmf2[ind]
        y = uvmhalo[z,0,ind]
        ydn = uvmhalo[z,1,ind]
        yup = uvmhalo[z,2,ind]
        ax.plot(x, y[0],color=cols[z], linestyle='solid', linewidth = 2, label=labels[z])
    common.prepare_legend(ax, cols[2:7], loc=0)

    ax = fig.add_subplot(212)
    ytit="$\\rm \\sigma(1500\\AA\, mag)$"

    ymin, ymax = -3, 3
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))
    for idx, z in zip(idxs, zs):
  
        #Predicted LF
        ind = np.where(uvmhalo[z,0,:] != 0)
        x = xmf2[ind]
        y = uvmhalo[z,0,ind]
        ydn = uvmhalo[z,1,ind] 
        yup = uvmhalo[z,2,ind] 
        ax.plot(x, -1 * ydn[0],color=cols[z], linestyle='solid', linewidth = 2)
        ax.plot(x, yup[0],color=cols[z], linestyle='solid', linewidth = 2)
    ax.plot([xmin, xmax], [0,0], linestyle='dotted', color='k')

    plt.tight_layout()
    common.savefig(outdir, fig, "UV_HaloMass_zGT6.pdf")



    fig = plt.figure(figsize=(5,7))
    xmin, xmax, ymin, ymax = 9, 14, 5, 12
    idxs = (2, 3, 4, 5, 6)
    zs  = (2, 3, 4, 5, 6)

    xtit="$\\rm log_{10}(M_{\\rm halo}/M_{\\odot})$"
    ytit="$\\rm log_{10}(M_{\\star}/M_{\\odot})$"

    ax = fig.add_subplot(211)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))

    for idx, z in zip(idxs, zs):
  
        #Predicted LF
        ind = np.where(mstarhalo[z,0,:] != 0)
        x = xmf[ind]
        y = mstarhalo[z,0,ind]
        ax.plot(x, y[0],color=cols[z], linestyle='solid', linewidth = 2, label=labels[z])
    common.prepare_legend(ax, cols[2:7], loc=0)

    ax = fig.add_subplot(212)

    ytit="$\\rm \\sigma(log_{10}(M_{\\star}))$"

    xmin, xmax, ymin, ymax = 5, 12, -2,2
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))

    for idx, z in zip(idxs, zs):
  
        #Predicted LF
        ind = np.where(mstarhalo[z,0,:] != 0)
        x = xmf[ind]
        ydn = mstarhalo[z,1,ind] 
        yup = mstarhalo[z,2,ind] 
        ax.plot(x, -1 * ydn[0],color=cols[z], linestyle='solid', linewidth = 2)
        ax.plot(x, yup[0],color=cols[z], linestyle='solid', linewidth = 2)
    ax.plot([xmin, xmax], [0,0], linestyle='dotted', color='k')

    plt.tight_layout()
    common.savefig(outdir, fig, "Mhalo_StellarMass_zGT6.pdf")
   
    fig = plt.figure(figsize=(5,7))
    xmin, xmax, ymin, ymax = 5, 12, -2, 3.5

    ytit="$\\rm log_{10}(SFR/M_{\\odot}\\, yr^{-1})$"
    xtit="$\\rm log_{10}(M_{\\star}/M_{\\odot})$"

    ax = fig.add_subplot(211)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))

    cols = ('Indigo','purple','Navy','Aquamarine', 'Green','Gold','red','DarkRed') 

    for idx, z in zip(idxs, zs):
  
        #Predicted LF
        ind = np.where(mstarsfr[z,0,:] != 0)
        x = xmf[ind]
        y = mstarsfr[z,0,ind]
        ax.plot(x, y[0],color=cols[z], linestyle='solid', linewidth = 2, label=labels[z])
    common.prepare_legend(ax, cols[2:7], loc=0)

    ax = fig.add_subplot(212)

    ytit="$\\rm \\sigma(log_{10}(SFR))$"

    xmin, xmax, ymin, ymax = 5, 12, -1.5, 1.5
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 0.5, 0.5))

    for idx, z in zip(idxs, zs):
  
        #Predicted LF
        ind = np.where(mstarsfr[z,0,:] != 0)
        x = xmf[ind]
        ydn = mstarsfr[z,1,ind] 
        yup = mstarsfr[z,2,ind] 
        ax.plot(x, -1 * ydn[0],color=cols[z], linestyle='solid', linewidth = 2)
        ax.plot(x, yup[0],color=cols[z], linestyle='solid', linewidth = 2)
    ax.plot([xmin, xmax], [0,0], linestyle='dotted', color='k')

    plt.tight_layout()
    common.savefig(outdir, fig, "StellarMass_SFR_zGT6.pdf")
   


def prepare_data(hdf5_data, phot_data, phot_data_nod, LFs_dust, LFs_nodust, index, nbands, uvmhalo, uvmstar, mstarhalo, mstarsfr):
  
    bin_it = functools.partial(us.wmedians, xbins=xmf, nmin=10)
    bin_it2 = functools.partial(us.wmedians, xbins=xmf2, nmin=10)

    #star_formation_histories and SharkSED have the same number of galaxies in the same order, and so we can safely assume that to be the case.
    #to select the same galaxies in galaxies.hdf5 we need to ask for all of those that have a stellar mass > 0, and then assume that they are in the same order.

    (h0, volh, mdisk, mbulge, mhalo, mshalo, typeg, age,
     sfr_disk, sfr_burst, id_gal) = hdf5_data

    #components:
    #(len(my_data), 2, 2, 5, nbands)
    #0: disk instability bulge
    #1: galaxy merger bulge
    #2: total bulge
    #3: disk
    #4: total
    seds_bulge_dummy = phot_data[0]
    ngals = len(seds_bulge_dummy[0,:])
    SEDs_dust = np.zeros(shape = (5,nbands,ngals))
    SEDs_dust[0,:] = phot_data[0]
    SEDs_dust[1,:] = phot_data[1]
    SEDs_dust[2,:] = phot_data[2]
    SEDs_dust[3,:] = phot_data[3]
    SEDs_dust[4,:] = phot_data[4]
    SEDs_nodust = np.zeros(shape = (5,nbands,ngals))
    SEDs_nodust[0,:] = phot_data_nod[0]
    SEDs_nodust[1,:] = phot_data_nod[1]
    SEDs_nodust[2,:] = phot_data_nod[2]
    SEDs_nodust[3,:] = phot_data_nod[3]
    SEDs_nodust[4,:] = phot_data_nod[4]

    #print(SEDs_dust.shape)
    for i in range(0,nbands):
        for c in range(0,5):
            #calculate LF with bands with dust
            ind = np.where(SEDs_dust[c,i,:] < -1)
            H, bins_edges = np.histogram(SEDs_dust[c,i,:],bins=np.append(mbins,mupp))
            LFs_dust[index,c,i,:] = LFs_dust[index,c,i,:] + H

            #calculate LF of intrinsic bands 
            ind = np.where(SEDs_nodust[c,i,:] < -1)
            H, bins_edges = np.histogram(SEDs_nodust[c,i,:],bins=np.append(mbins,mupp))
            LFs_nodust[index,c,i,:] = LFs_nodust[index,c,i,:] + H

    ind = np.where((mdisk + mbulge > 0) & (SEDs_nodust[4,0,:] < -5) & (SEDs_nodust[4,0,:] > -40))
    uvmstar[index,:] = bin_it(x=np.log10((mdisk[ind] + mbulge[ind])/h0), y=SEDs_dust[4,0,ind][0])
    uvmhalo[index,:] = bin_it2(x=np.log10((mhalo[ind])/h0), y=SEDs_dust[4,0,ind][0])
    mstarhalo[index,:] = bin_it(x=np.log10((mhalo[ind])/h0), y=np.log10((mdisk[ind] + mbulge[ind])/h0))
    mstarsfr[index,:] = bin_it(x=np.log10((mdisk[ind] + mbulge[ind])/h0), y=np.log10((sfr_disk[ind]+sfr_burst[ind])/1e9/h0))

    if index == 6:
        print(np.log10((mdisk[ind] + mbulge[ind])/h0))
        print(uvmstar[index,:], uvmhalo[index,:])

def main(model_dir, outdir, redshift_table, subvols, obsdir):

    # Loop over redshift and subvolumes
    plt = common.load_matplotlib()

    Variable_Ext = True

    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mvir_hosthalo',
                           'mvir_subhalo', 'type', 'mean_stellar_age',
                           'sfr_disk', 'sfr_burst', 'id_galaxy')}

    fields_sed = {'SED/ab_dust': ('bulge_d','bulge_m','bulge_t','disk','total'),}
    fields_sed_nod = {'SED/ab_nodust': ('bulge_d','bulge_m','bulge_t','disk','total')}

    z = (3.0, 4.0, 6.0, 8.0, 10.0, 13, 15)  #(8.0, 9, 10.5, 12.5) #, 1.0, 1.5, 2.0)
    snapshots = redshift_table[z]
    uvmhalo = np.zeros(shape = (len(z), 3, len(xmf2)))
    uvmstar = np.zeros(shape = (len(z), 3, len(xmf)))
    mstarhalo = np.zeros(shape = (len(z), 3, len(xmf)))
    mstarsfr = np.zeros(shape = (len(z), 3, len(xmf)))

    file_hdf5_sed = "Shark-SED-JWST-eagle-rr14-steep.hdf5" 
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

        prepare_data(hdf5_data, seds, seds_nod, LFs_dust, LFs_nodust, index, nbands, uvmhalo, uvmstar, mstarhalo, mstarsfr)

        h0, volh = hdf5_data[0], hdf5_data[1]
        if(volh > 0.):
            LFs_dust[index,:]   = LFs_dust[index,:]/volh
            LFs_nodust[index,:] = LFs_nodust[index,:]/volh

    print ("number of bands %d" % (nbands))
    # Take logs
    ind = np.where(LFs_dust > 0.)
    LFs_dust[ind] = np.log10(LFs_dust[ind])

    ind = np.where(LFs_nodust > 0.)
    LFs_nodust[ind] = np.log10(LFs_nodust[ind])

    if(Variable_Ext):
       outdir = os.path.join(outdir, 'eagle-rr14-steep')

    volcorr = 3.0*np.log10(h0)

    band = 0
    for i,zin in enumerate(z):
       print("#UV LF at redshift:", zin)
       print("#rest_Frame_mag num_Density_dust num_density_nodust")
       for a,b,c in zip(xlf, LFs_dust[i,4,band,:], LFs_nodust[i,4,band,:]):
           print(a,b+volcorr-np.log10(dm),c+volcorr-np.log10(dm))
    plot_uv_lf_evo(plt, outdir, obsdir, h0, LFs_dust, LFs_nodust, nbands)
    plot_uv_mass_evo(plt, outdir, obsdir, uvmhalo, uvmstar, mstarhalo, mstarsfr)

if __name__ == '__main__':
    main(*common.parse_args())
