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


mflow = 5
mfupp = 14
dmf = 0.5
mfbins = np.arange(mflow,mfupp,dmf)
xmf = mfbins + dmf/2.0

mflow2 = 9
mfupp2 = 14
dmf2 = 0.25
mfbins2 = np.arange(mflow2,mfupp2,dmf2)
xmf2 = mfbins2 + dmf2/2.0

def plot_uv_lf_micro(ax, z, all_channels = False, label = True):
    x, y, yun, ydi, ym, yd, redshift = np.loadtxt("../data/Models/SharkVariations/uv_lf_L24HBTbestparams_micro-SURFS_allz.txt", unpack = True, usecols = [0,1,2,3,4,5,6])

    if(all_channels == False):
       if(z != 17):
          ind = np.where((y < -0.5) & (redshift == z))
          ax.plot(x[ind], y[ind], color='Purple', linewidth=3, label = 'Total att. (micro)' if label == True else None)

          ind = np.where((yun < -0.5) & (redshift == z))
          ax.plot(x[ind], yun[ind],color='Purple', linewidth=3, linestyle='dashed', label='Total unatt. (micro)' if label == True else None)
       else:
          ind = np.where((y < -0.5) & (redshift == z))
          ax.plot(x[ind], y[ind], color='Purple', linewidth=1, label = 'Total att. (micro)' if label == True else None)

    if(all_channels):
       ind = np.where((ydi < -0.5) & (redshift == z))
       ax.plot(x[ind], ydi[ind],'Navy', linewidth=2, linestyle='dashed', label = 'SF disks (micro)' if label == True else None)
       ind = np.where((ym < -0.5) & (redshift == z))
       ax.plot(x[ind], ym[ind], 'DarkSeaGreen', linewidth=2, linestyle='dashed', label='SBs disk ins. (micro)' if label == True else None)
       ind = np.where((yd < -0.5) & (redshift == z))
       ax.plot(x[ind], yd[ind],'Firebrick', linewidth=2, linestyle='dashed', label='SBs mergers (micro)' if label == True else None)
   
def plot_uv_lf_l800(ax, z, all_channels = False, label = True):
    x, y, yun, ydi, ym, yd = np.loadtxt("../data/Models/SharkVariations/uv_lf_L24HBTbestparams_L800_z" + z + ".txt", unpack = True, usecols = [0,1,2,3,4,5])

    if(all_channels == False):
       if(z != "17"):
          ind = np.where(y != -10)
          ax.plot(x[ind], y[ind], color='Olive', linewidth=3, label = 'Total att. (L800)' if label == True else None)

          ind = np.where(yun != -10)
          ax.plot(x[ind], yun[ind],color='Olive', linewidth=3, linestyle='dashed', label='Total unatt. (L800)' if label == True else None)
       else:
          ind = np.where(y != -10)
          ax.plot(x[ind], y[ind], color='Olive', linewidth=1, label = 'Total att. (L800)' if label == True else None)


    if(all_channels):
       ind = np.where(ydi != -10)
       ax.plot(x[ind], ydi[ind],'Navy', linewidth=2, linestyle='dotted', label = 'SF disks (L800)' if label == True else None)
       ind = np.where(ym != -10)
       ax.plot(x[ind], ym[ind], 'DarkSeaGreen', linewidth=2, linestyle='dotted', label='SBs disk ins. (L800)' if label == True else None)
       ind = np.where(yd != -10)
       ax.plot(x[ind], yd[ind],'Firebrick', linewidth=2, linestyle='dotted', label='SBs mergers (L800)' if label == True else None)
   


def plot_uv_lf_evo(plt, outdir, obsdir, h0, LFs_dust, LFs_nodust, fracs, nbands):

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
    band = 17 #FUV_Nathan
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

           plot_uv_lf_micro(ax, 6)
           plot_uv_lf_l800(ax, '6')

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
           plot_uv_lf_micro(ax, 8, label = False)
           plot_uv_lf_l800(ax, '8', label = False)

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
           plot_uv_lf_micro(ax, 10, label = False)
           plot_uv_lf_l800(ax, '10', label = False)

        if(idx == 5):
           file = obsdir+'/lf/lf1500_z13_weibel25.data'
           lm,p,dpup,dpdn,flag = np.loadtxt(file,usecols=[0,1,2,3,4],unpack=True)
           ax.errorbar(lm+corrm_obs, p+corry_obs, yerr=[dpup, dpdn], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')
           file = obsdir+'/lf/lf1500_z12p8_whitler25.data'
           lm,p,dpdn,dpup = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
           #ax.errorbar(lm+corrm_obs, np.log10(p)+corry_obs, yerr=[np.log10(p) - np.log10(p-dpdn), np.log10(p + dpup) - np.log10(p)], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s')
           plot_uv_lf_micro(ax, 13, label = False)
           plot_uv_lf_l800(ax, '13', label = False)

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
           plot_uv_lf_micro(ax, 15, label = False)
           plot_uv_lf_micro(ax, 17, label = False)
           plot_uv_lf_l800(ax, '15', label = False)
           plot_uv_lf_l800(ax, '17', label = False)

        #Predicted LF
        ind = np.where(LFs_dust[z,4,band,:] < 0.)
        y = LFs_dust[z,4,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'k', linewidth=3, label='Total att. (medi)' if idx == 2 else None)
        ind = np.where(LFs_nodust[z,4,band,:] < 0.)
        y = LFs_nodust[z,4,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'k', linewidth=3, linestyle='dashed', label = 'Total unatt. (medi)' if idx == 2 else None)
        #if(idx == 1):
        #   for a,b,c in zip(xlf_obs,LFs_dust[z,4,band,:],LFs_nodust[z,4,band,:]):
        #       print (a, b+volcorr-np.log10(dm),c+volcorr-np.log10(dm))

        #ind = np.where(LFs_dust[z,3,band,:] < 0.)
        #y = LFs_dust[z,3,band,ind]+volcorr-np.log10(dm)
        #ax.plot(xlf_obs[ind],y[0],'LightSeaGreen', linewidth=4, linestyle='dotted', label='star formation disks' if idx == 2 else None)
        #ind = np.where(LFs_dust[z,0,band,:] < 0.)
        #y = LFs_dust[z,0,band,ind]+volcorr-np.log10(dm)
        #ax.plot(xlf_obs[ind],y[0],'OrangeRed', linewidth=4, linestyle='dashed', label='starbursts (disk ins.)' if idx == 2 else None)
        #ind = np.where(LFs_dust[z,1,band,:] < 0.)
        #y = LFs_dust[z,1,band,ind]+volcorr-np.log10(dm)
        #ax.plot(xlf_obs[ind],y[0],'Orange', linewidth=4, linestyle='dashdot', label='starbursts (mergers)' if idx == 2 else None)

        if(idx == 2):
            common.prepare_legend(ax, ['Purple','Purple','Olive','Olive','k','k','grey','grey','grey'],  bbox_to_anchor=(2.305, -0.98))
        if idx == 3:
            common.prepare_legend(ax, ['grey','grey','grey'], bbox_to_anchor=(1.11, -1.06))
        if idx == 4:
            common.prepare_legend(ax, ['grey','grey','grey'], bbox_to_anchor=(-0.09, -1.22))
        if idx == 6:
            common.prepare_legend(ax, ['grey','grey','grey'], bbox_to_anchor=(1.11, -0.095))


        #if ((idx ==4) or (idx == 3)):
        #    common.prepare_legend(ax, ['grey','grey','grey'], loc=0)

    #plt.tight_layout()
    common.savefig(outdir, fig, "UV_luminosity_function_evolution_zGT6.pdf")

    fig = plt.figure(figsize=(12,8))
    xmin, xmax, ymin, ymax = -25, -13, -7, -1

    subplots = (231, 232, 233, 234, 235)
    idx = (2, 3, 4, 5, 6)
    zs  = (2, 3, 4, 5, 6)

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
           ax.errorbar(lmB15z6+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='^') #, label="Bouwens+2015")

           yobs = np.log10(pF6*1e-3)
           ydn  = np.log10(pF6*1e-3-dpdF6*1e-3)
           yup  = np.log10(pF6*1e-3+dpuF6*1e-3)
           ax.errorbar(lmF15+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='v')#, label="Finkelstein+2015")

           plot_uv_lf_micro(ax, '6', all_channels = True)
           plot_uv_lf_l800(ax, '6', all_channels = True)

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
           plot_uv_lf_micro(ax, '8', all_channels = True, label = False)
           plot_uv_lf_l800(ax, '8', all_channels = True, label = False)

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
           plot_uv_lf_micro(ax, '10', all_channels = True, label = False)
           plot_uv_lf_l800(ax, '10', all_channels = True, label = False)

        if(idx == 5):
           file = obsdir+'/lf/lf1500_z13_weibel25.data'
           lm,p,dpup,dpdn,flag = np.loadtxt(file,usecols=[0,1,2,3,4],unpack=True)
           ax.errorbar(lm+corrm_obs, p+corry_obs, yerr=[dpup, dpdn], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')
           file = obsdir+'/lf/lf1500_z12p8_whitler25.data'
           lm,p,dpdn,dpup = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
           #ax.errorbar(lm+corrm_obs, np.log10(p)+corry_obs, yerr=[np.log10(p) - np.log10(p-dpdn), np.log10(p + dpup) - np.log10(p)], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s')
           plot_uv_lf_micro(ax, '13', all_channels = True, label = False)

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
           plot_uv_lf_micro(ax, '15', all_channels = True, label = False)
           plot_uv_lf_l800(ax, '15', all_channels = True, label = False)

        #Predicted LF
        ind = np.where(LFs_dust[z,3,band,:] < 0.)
        y = LFs_dust[z,3,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'Navy', linewidth=2, linestyle='solid', label='SF disks (medi)' if idx == 2 else None)
        ind = np.where(LFs_dust[z,0,band,:] < 0.)
        y = LFs_dust[z,0,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'DarkSeaGreen', linewidth=2, linestyle='solid', label='SBs disk ins. (medi)' if idx == 2 else None)
        ind = np.where(LFs_dust[z,1,band,:] < 0.)
        y = LFs_dust[z,1,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'Firebrick', linewidth=2, linestyle='solid', label='SBs mergers (medi)' if idx == 2 else None)

        if(idx == 2):
            common.prepare_legend(ax, ['Navy','DarkSeaGreen','Firebrick', 'Navy','DarkSeaGreen','Firebrick','Navy','DarkSeaGreen','Firebrick', 'grey','grey','grey'],  bbox_to_anchor=(2.305, -1.2))
        #if idx == 3:
        #    common.prepare_legend(ax, ['grey','grey','grey'], bbox_to_anchor=(1.11, -1.06))
        #if idx == 4:
        #    common.prepare_legend(ax, ['grey','grey','grey'], bbox_to_anchor=(-0.09, -1.23))
        #
        #if idx == 6:
        #    common.prepare_legend(ax, ['grey','grey','grey'], bbox_to_anchor=(1.11, -0.105))

        #if ((idx ==4) or (idx == 3)):
        #    common.prepare_legend(ax, ['grey','grey','grey'], loc=0)

    #plt.tight_layout()
    common.savefig(outdir, fig, "UV_luminosity_function_evolution_zGT6_allchannels.pdf")




    fig = plt.figure(figsize=(12,8))
    xmin, xmax, ymin, ymax = -25, -13, 0, 1

    subplots = (231, 232, 233, 234, 235)
    idx = (2, 3, 4, 5, 6)
    zs  = (2, 3, 4, 5, 6)
    ytit="Fractional channel contribution"

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

        common.prepare_ax(ax, -23, xmax, ymin, ymax, xtitplot, ytitplot, locators=(2, 2, 0.1, 0.1))
        ax.text(-15,0.9, labels[idx])
        
        ylow = np.zeros(shape = len(xlf))
        ax.fill_between(xlf, ylow, fracs[idx,:,0], facecolor='navy', alpha=0.15, interpolate=True, label = 'SF disks') 
        ax.fill_between(xlf, fracs[idx,:,0], fracs[idx,:,0] + fracs[idx,:,1], facecolor='DarkSeaGreen', alpha=0.15, interpolate=True, label = 'SBs disk ins.') 
        yhigh = np.zeros(shape = len(xlf))
        yhigh[:] = 1.0
        ax.fill_between(xlf, fracs[idx,:,0] + fracs[idx,:,1], yhigh, facecolor='Firebrick', alpha=0.15, interpolate=True, label = 'SBs mergers') 

        if(idx == 2):
            common.prepare_legend(ax, ['Navy','DarkSeaGreen','Firebrick', 'Navy','DarkSeaGreen','Firebrick','Navy','DarkSeaGreen','Firebrick', 'grey','grey','grey'],  bbox_to_anchor=(2.305, -0.8))
        if(idx == 6):
            ax.text(-10,0.75, "Shark v2.0 (L800)", fontsize=12)

    #plt.tight_layout()
    common.savefig(outdir, fig, "UV_luminosity_function_channels_contribution_zGT6.pdf")

def plot_uv_mass_evo(plt, outdir, obsdir, uvmhalo, uvmstar, mstarhalo, mstarsfr, mstarbeta):

    fig = plt.figure(figsize=(5,7))
    xmin, xmax, ymin, ymax = 4.5, 11, -10, -23

    ytit="$\\rm1500\\AA\, mag\, (AB)$"
    xtit="$\\rm log_{10}(M_{\\star}/M_{\\odot})$"

    idxs = (2, 3, 4, 5, 6, 7)
    zs  = (2, 3, 4, 5, 6, 7)
    labels= ('z=3', 'z=4', 'z=6', 'z=8', 'z=10', 'z=13', 'z=15', 'z=17')
    redshift=(3,4,6,8,10,13,15,17)
    ax = fig.add_subplot(211)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))

    cols = ('Indigo','purple','Navy','Aquamarine', 'Green','Gold','red','DarkRed') 

    sm_micro,luv_micro,suvdn_micro,suvup_micro,z_micro = np.loadtxt("../data/Models/SharkVariations/uv_stellarmass_L24HBTbestparams_micro-SURFS.txt", unpack = True, usecols = [0,1,2,3,4])
    for idx, z in zip(idxs, zs):
  
        #Predicted LF
        ind = np.where(uvmstar[z,0,:] != 0)
        x = xmf[ind]
        y = uvmstar[z,0,ind]
        ydn = uvmstar[z,1,ind]
        yup = uvmstar[z,2,ind]
        #ax.fill_between(x, ydn[0], yup[0], color=cols[z], alpha = 0.25, linestyle='solid', linewidth = 2)
        ax.plot(x, y[0],color=cols[z], linestyle='solid', linewidth = 4, label=labels[z])
        ind=np.where(z_micro == redshift[z])
        ax.plot(sm_micro[ind], luv_micro[ind], color=cols[z],linestyle='solid', linewidth = 1)

    common.prepare_legend(ax, cols[2:7], loc=0)

    ax = fig.add_subplot(212)

    ytit="$\\rm \\sigma(1500\\AA\, mag)$"

    ymin, ymax = -2,2
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))

    for idx, z in zip(idxs, zs):
  
        #Predicted LF
        ind = np.where(uvmstar[z,0,:] != 0)
        x = xmf[ind]
        y = uvmstar[z,0,ind]
        ydn = uvmstar[z,1,ind] 
        yup = uvmstar[z,2,ind] 
        ax.plot(x, -1 * ydn[0],color=cols[z], linestyle='solid', linewidth = 4)
        ax.plot(x, yup[0],color=cols[z], linestyle='solid', linewidth = 4)

        ind=np.where(z_micro == redshift[z])
        ax.plot(sm_micro[ind], suvdn_micro[ind], color=cols[z], linestyle='solid', linewidth = 1)
        ax.plot(sm_micro[ind], suvup_micro[ind], color=cols[z], linestyle='solid', linewidth = 1)

    ax.plot([xmin, xmax], [0,0], linestyle='dotted', color='k')

    plt.tight_layout()
    common.savefig(outdir, fig, "UV_StellarMass_zGT6.pdf")
   
    fig = plt.figure(figsize=(5,7))
    xmin, xmax, ymin, ymax = 9, 13, -12, -24

    ytit="$\\rm 1500\\AA\, mag\, (AB)$"
    xtit="$\\rm log_{10}(M_{\\rm halo}/M_{\\odot})$"
    #idxs = (2, 3, 4, 5)
    #zs  = (2, 3, 4, 5)

    ax = fig.add_subplot(211)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))
    sm_micro,luv_micro,suvdn_micro,suvup_micro,z_micro = np.loadtxt("../data/Models/SharkVariations/uv_halomass_L24HBTbestparams_micro-SURFS.txt", unpack = True, usecols = [0,1,2,3,4])

    for idx, z in zip(idxs, zs):
  
        #Predicted LF
        ind = np.where(uvmhalo[z,0,:] != 0)
        x = xmf2[ind]
        y = uvmhalo[z,0,ind]
        ydn = uvmhalo[z,1,ind]
        yup = uvmhalo[z,2,ind]
        ax.plot(x, y[0],color=cols[z], linestyle='solid', linewidth = 4, label=labels[z])
        ind=np.where(z_micro == redshift[z])
        ax.plot(sm_micro[ind], luv_micro[ind], color=cols[z],linestyle='solid', linewidth = 1)

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
        ax.plot(x, -1 * ydn[0],color=cols[z], linestyle='solid', linewidth = 4)
        ax.plot(x, yup[0],color=cols[z], linestyle='solid', linewidth = 4)
        ind=np.where(z_micro == redshift[z])
        ax.plot(sm_micro[ind], suvdn_micro[ind], color=cols[z], linestyle='solid', linewidth = 1)
        ax.plot(sm_micro[ind], suvup_micro[ind], color=cols[z], linestyle='solid', linewidth = 1)

    ax.plot([xmin, xmax], [0,0], linestyle='dotted', color='k')

    plt.tight_layout()
    common.savefig(outdir, fig, "UV_HaloMass_zGT6.pdf")



    fig = plt.figure(figsize=(5,7))
    xmin, xmax, ymin, ymax = 9, 14, 5, 12
    idxs = (2, 3, 4, 5, 6, 7)
    zs  = (2, 3, 4, 5, 6, 7)

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
    xmin, xmax, ymin, ymax = 4.5, 11, -3, 3.5

    ytit="$\\rm log_{10}(SFR/M_{\\odot}\\, yr^{-1})$"
    xtit="$\\rm log_{10}(M_{\\star}/M_{\\odot})$"

    ax = fig.add_subplot(211)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))

    cols = ('Indigo','purple','Navy','Aquamarine', 'Green','Gold','red','DarkRed') 
    sm_micro,luv_micro,suvdn_micro,suvup_micro,z_micro = np.loadtxt("../data/Models/SharkVariations/sfr_stellarmass_L24HBTbestparams_micro-SURFS.txt", unpack = True, usecols = [0,1,2,3,4])

    for idx, z in zip(idxs, zs):
  
        #Predicted LF
        ind = np.where(mstarsfr[z,0,:] != 0)
        x = xmf[ind]
        y = mstarsfr[z,0,ind]
        ydn = mstarsfr[z,1,ind]
        yup = mstarsfr[z,2,ind]
        ax.plot(x, y[0],color=cols[z], linestyle='solid', linewidth = 4, label=labels[z])
        print("#SFR vs stellar mass relation at redshift", labels[z])
        for a,b,c,d in zip(x, y[0], -1 * ydn[0], yup[0]):
            print(a,b,c,d)
        ind=np.where(z_micro == redshift[z])
        ax.plot(sm_micro[ind], luv_micro[ind], color=cols[z],linestyle='solid', linewidth = 1)
    common.prepare_legend(ax, cols[2:7], loc=0)

    ax = fig.add_subplot(212)

    ytit="$\\rm \\sigma(log_{10}(SFR))$"

    ymin, ymax = -1.5, 1.5
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 0.5, 0.5))

    for idx, z in zip(idxs, zs):
  
        #Predicted LF
        ind = np.where(mstarsfr[z,0,:] != 0)
        x = xmf[ind]
        ydn = mstarsfr[z,1,ind] 
        yup = mstarsfr[z,2,ind] 
        ax.plot(x, -1 * ydn[0],color=cols[z], linestyle='solid', linewidth = 4)
        ax.plot(x, yup[0],color=cols[z], linestyle='solid', linewidth = 4)
        ind=np.where(z_micro == redshift[z])
        ax.plot(sm_micro[ind], suvdn_micro[ind], color=cols[z], linestyle='solid', linewidth = 1)
        ax.plot(sm_micro[ind], suvup_micro[ind], color=cols[z], linestyle='solid', linewidth = 1)

    ax.plot([xmin, xmax], [0,0], linestyle='dotted', color='k')

    plt.tight_layout()
    common.savefig(outdir, fig, "StellarMass_SFR_zGT6.pdf")
   
    fig = plt.figure(figsize=(5,4))
    xmin, xmax, ymin, ymax = 5, 11, 0.1, 500

    ytit="$\\rm \\dot{M}_{\\rm out}/SFR$"
    xtit="$\\rm log_{10}(M_{\\star}/M_{\\odot})$"

    ax = fig.add_subplot(111)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1e2, 1e2))
    ax.set_yscale('log')
    cols = ('Indigo','purple','Navy','Aquamarine', 'Green','Gold','red','DarkRed') 

    for idx, z in zip(idxs, zs):
  
        #Predicted LF
        ind = np.where(mstarbeta[z,0,:] != 0)
        x = xmf[ind]
        y = mstarbeta[z,0,ind]
        ydn = mstarbeta[z,1,ind]
        yup = mstarbeta[z,2,ind]
        ax.plot(x, y[0],color=cols[z], linestyle='solid', linewidth = 4, label=labels[z])
        ax.plot(x, y[0] - ydn[0],color=cols[z], linestyle='dotted', linewidth = 4)
        ax.plot(x, y[0] + yup[0],color=cols[z], linestyle='dotted', linewidth = 4)

    common.prepare_legend(ax, cols[2:7], loc=0)

    plt.tight_layout()
    common.savefig(outdir, fig, "StellarMass_MassLoading_zGT6.pdf")
   
def prepare_data(hdf5_data, phot_data, phot_data_nod, LFs_dust, LFs_nodust, index, nbands, uvmhalo, uvmstar, uvmhalo_in, uvmstar_in, mstarhalo, mstarsfr, mstarbeta, fracs, zlist):
  
    bin_it = functools.partial(us.wmedians, xbins=xmf, nmin=10)
    bin_it2 = functools.partial(us.wmedians, xbins=xmf2, nmin=10)
    bin_it_lf= functools.partial(us.wmedians, xbins=xlf, nmin=10)

    #star_formation_histories and SharkSED have the same number of galaxies in the same order, and so we can safely assume that to be the case.
    #to select the same galaxies in galaxies.hdf5 we need to ask for all of those that have a stellar mass > 0, and then assume that they are in the same order.

    (h0, volh, mdisk, mbulge, mhalo, mshalo, typeg, age,
     sfr_disk, sfr_burst, id_gal, vvir_halo) = hdf5_data
    volcorr = 3.0*np.log10(h0)
    lgvolh = np.log10(volh)

    redshift_power = 0.13405459588358698
    v_sn = 120
    beta_disk = 3.79746174188
    eps_halo = 2.0
    min_beta = 0.104050197191
   
    age_univ =  us.look_back_time(zlist[index])
    vhot = v_sn * (age_univ)**redshift_power
    const_sn =  (vhot/vvir_halo)**beta_disk
    ind = np.where(const_sn < min_beta)
    const_sn[ind] = min_beta

    ind = np.where(mdisk + mbulge > 0)
    mdisk = mdisk[ind] 
    mbulge = mbulge[ind]
    mhalo = mhalo[ind]
    mshalo = mshalo[ind]
    typeg = typeg[ind]
    sfr_disk = sfr_disk[ind]
    sfr_burst = sfr_burst[ind]
    id_gal = id_gal[ind]
    vvir_halo = vvir_halo[ind]

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

   #(0): "F070W_JWST", "F090W_JWST", "F115W_JWST", "F150W_JWST", "F200W_JWST",
   #(5): "F277W_JWST", "F356W_JWST", "F444W_JWST", "F560W_JWST", "F770W_JWST",
   #(10): "F1000W_JWST", "F1130W_JWST", "F1280W_JWST", "F1500W_JWST",
   #(14): "F1800W_JWST", "F2100W_JWST", "F2550W_JWST",
   #(17): "Band_ionising_photons", "FUV_Nathan",
    band = 17 #FUV_Nathan

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
    print("#UV LF at redshift:", zlist[index])
    print("#rest_Frame_mag num_Density_dust num_density_nodust num_Density_dust_diskins num_Density_dust_mergers num_Density_dust_disks")
    for a,b,c,d,e,f in zip(xlf, np.log10(LFs_dust[index,4,band,:]), np.log10(LFs_nodust[index,4,band,:]),  np.log10(LFs_dust[index,0,band,:]),  np.log10(LFs_dust[index,1,band,:]),  np.log10(LFs_dust[index,3,band,:])):
        print(a,b+volcorr-np.log10(dm)-lgvolh,c+volcorr-np.log10(dm)-lgvolh, d+volcorr-np.log10(dm)-lgvolh, e+volcorr-np.log10(dm)-lgvolh, f+volcorr-np.log10(dm)-lgvolh, zlist[index])

    for i, m in enumerate(xlf):
        ind = np.where((SEDs_dust[4,band,:] < m + dm/2) & (SEDs_dust[4,band,:] >= m - dm/2))
        med_ldisk = np.median(SEDs_dust[3,band,ind])
        med_lbd = np.median(SEDs_dust[0,band,ind])
        med_lbs = np.median(SEDs_dust[1,band,ind])
        if(med_ldisk != -999):
            med_ldisk = 10**((med_ldisk + 48.6)/(-2.5))
        else:
            med_ldisk = 0
        if(med_lbd != -999):
            med_lbd = 10**((med_lbd + 48.6)/(-2.5))
        else:
            med_lbd = 0
        if(med_lbs != -999):
            med_lbs = 10**((med_lbs + 48.6)/(-2.5))
        else:
            med_lbs = 0
        fracs[index,i,0] = med_ldisk / (med_ldisk + med_lbd + med_lbs)
        fracs[index,i,1] = med_lbd / (med_ldisk + med_lbd + med_lbs)
        fracs[index,i,2] = med_lbs / (med_ldisk + med_lbd + med_lbs)
 
    print("#UV mag fractional contribution of channels at redshift", zlist[index])
    for a,b,c,d in zip(xlf, fracs[index,:,0], fracs[index,:,1],fracs[index,:,2]):
        print(a,b,c,d,zlist[index])
    ind = np.where((SEDs_nodust[4, band,:] < -5) & (SEDs_nodust[4, band,:] > -40))
    uvmstar[index,:] = bin_it(x=np.log10((mdisk[ind] + mbulge[ind])/h0), y=SEDs_dust[4, band, ind][0])
    uvmstar_in[index,:] = bin_it_lf(x=SEDs_dust[4, band, ind][0], y=np.log10((mdisk[ind] + mbulge[ind])/h0))

    #uvmhalo[index,:] = bin_it2(x=np.log10((mhalo[ind])/h0), y=SEDs_dust[4,0,ind][0])
    mstarhalo[index,:] = bin_it(x=np.log10((mhalo[ind])/h0), y=np.log10((mdisk[ind] + mbulge[ind])/h0))
    mstarsfr[index,:] = bin_it(x=np.log10((mdisk[ind] + mbulge[ind])/h0), y=np.log10((sfr_disk[ind]+sfr_burst[ind])/1e9/h0))
    mstarbeta[index,:] = bin_it(x=np.log10((mdisk[ind] + mbulge[ind])/h0), y=const_sn[ind])

    ind = np.where((mdisk + mbulge > 0) & (SEDs_nodust[4,band,:] < -5) & (SEDs_nodust[4,band,:] > -40) & (typeg == 0))
    uvmhalo[index,:] = bin_it2(x=np.log10((mhalo[ind])/h0), y=SEDs_dust[4,band,ind][0])
    uvmhalo_in[index,:] = bin_it_lf(x=SEDs_dust[4,band,ind][0], y=np.log10((mhalo[ind])/h0))

    print("#UV mag vs stellar mass relation at redshift", zlist[index], " binned by UV mag")
    for a,b,c,d in zip(xlf, uvmstar_in[index,0,:], uvmstar_in[index,1,:], uvmstar_in[index,2,:]):
        print(a,b,c,d, zlist[index])
    print("#UV mag vs halo mass relation at redshift",  zlist[index], " binned by UV mag")
    for a,b,c,d in zip(xlf, uvmhalo_in[index,0,:],  uvmhalo_in[index,1,:], uvmhalo_in[index,2,:]):
        print(a,b,c,d, zlist[index])

def main(model_dir, outdir, redshift_table, subvols, obsdir):

    # Loop over redshift and subvolumes
    plt = common.load_matplotlib()

    Variable_Ext = True
    multiple_batches = False

    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mvir_hosthalo',
                           'mvir_subhalo', 'type', 'mean_stellar_age',
                           'sfr_disk', 'sfr_burst', 'id_galaxy', 'vvir_subhalo')}

    fields_sed = {'SED/ab_dust': ('bulge_d','bulge_m','bulge_t','disk','total'),}
    fields_sed_nod = {'SED/ab_nodust': ('bulge_d','bulge_m','bulge_t','disk','total')}

    #z = (6, 6, 6.0, 8.0, 10.0, 13.0, 15.0, 17.0)  #(8.0, 9, 10.5, 12.5) #, 1.0, 1.5, 2.0)
    #z = (6.0, 6.0, 6.0, 8.0, 10.0, 12.7, 15.0, 17.0)  #(8.0, 9, 10.5, 12.5) #, 1.0, 1.5, 2.0)
    #z = (17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0)  #(8.0, 9, 10.5, 12.5) #, 1.0, 1.5, 2.0)
    #z = [17.0]
    z = [17]
    snapshots = redshift_table[z]
    uvmhalo = np.zeros(shape = (len(z), 3, len(xmf2)))
    uvmstar = np.zeros(shape = (len(z), 3, len(xmf)))
    uvmhalo_in = np.zeros(shape = (len(z), 3, len(xlf)))
    uvmstar_in = np.zeros(shape = (len(z), 3, len(xlf)))

    mstarhalo = np.zeros(shape = (len(z), 3, len(xmf)))
    mstarsfr = np.zeros(shape = (len(z), 3, len(xmf)))
    mstarbeta = np.zeros(shape = (len(z), 3, len(xmf)))
    fracs = np.zeros(shape = (len(z), len(xlf), 3))

    file_hdf5_sed = "Shark-SED-JWST-eagle-rr14.hdf5" 
    # Create histogram
    for index, snapshot in enumerate(snapshots):
        if(multiple_batches == False):
            hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        else:
            hdf5_data = common.read_data_multiple_batches(model_dir, snapshot, fields)


        if(Variable_Ext == False):
           seds = common.read_photometry_data(model_dir, snapshot, fields_sed, subvols)
           seds_nod = common.read_photometry_data(model_dir, snapshot, fields_sed_nod, subvols)
        else:
           if(multiple_batches == False):
               seds = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed, subvols, file_hdf5_sed)
               seds_nod = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed_nod, subvols, file_hdf5_sed)
           else:
               seds = common.read_photometry_data_variable_tau_screen_multiple_batches(model_dir, snapshot, fields_sed, file_hdf5_sed)
               seds_nod = common.read_photometry_data_variable_tau_screen_multiple_batches(model_dir, snapshot, fields_sed_nod, file_hdf5_sed)

        nbands = len(seds[0]) 

        if(index == 0):
            LFs_dust     = np.zeros(shape = (len(z), 5, nbands, len(mbins)))
            LFs_nodust   = np.zeros(shape = (len(z), 5, nbands, len(mbins)))

        prepare_data(hdf5_data, seds, seds_nod, LFs_dust, LFs_nodust, index, nbands, uvmhalo, uvmstar, uvmhalo_in, uvmstar_in, mstarhalo, mstarsfr, mstarbeta, fracs, z)

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
       outdir = os.path.join(outdir, 'eagle-rr14')

    volcorr = 3.0*np.log10(h0)

    band = 17
    for i,zin in enumerate(z):
       print("#UV LF at redshift:", zin)
       print("#rest_Frame_mag num_Density_dust num_density_nodust num_Density_dust_diskins num_Density_dust_mergers num_Density_dust_disks")
       for a,b,c,d,e,f in zip(xlf, LFs_dust[i,4,band,:], LFs_nodust[i,4,band,:],  LFs_dust[i,0,band,:],  LFs_dust[i,1,band,:],  LFs_dust[i,3,band,:]):
           print(a,b+volcorr-np.log10(dm),c+volcorr-np.log10(dm), d+volcorr-np.log10(dm), e+volcorr-np.log10(dm), f+volcorr-np.log10(dm))
    plot_uv_lf_evo(plt, outdir, obsdir, h0, LFs_dust, LFs_nodust, fracs, nbands)
    plot_uv_mass_evo(plt, outdir, obsdir, uvmhalo, uvmstar, mstarhalo, mstarsfr, mstarbeta)

if __name__ == '__main__':
    main(*common.parse_args())
