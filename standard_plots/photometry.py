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

import common


##################################

# Constants
GyrToYr = 1e9
Zsun = 0.0127
XH = 0.72
PI = 3.141592654
MpcToKpc = 1e3

# Mass function initialization

mlow = -30 + 5.0 * np.log10(0.677)
mupp = -10 + 5.0 * np.log10(0.677)
dm = 0.5
mbins = np.arange(mlow,mupp,dm)
xlf   = mbins + dm/2.0

# colour distribution initialization
clow  = -0.1
cupp  = 3.5
dc    = 0.075
cbins = np.arange(clow,cupp,dc)
xc    = cbins + dc/2.0

magbins = [-17.13,-17.88,-18.63,-19.38,-20.13,-20.88,-21.63]

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
    labels= ('GALEX FUV', 'GALEX NUV', 'SDSS u', 'SDSS g', 'SDSS r', 'SDSS i', 'VISTA Z')
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


    xtit="$\\rm mag-5log(h) (AB)$"
    ytit="$\\rm log_{10}(\Phi/dex^{-1} h^3 {\\rm Mpc}^{-3})$"

    xmin, xmax, ymin, ymax = -30, -16, -5, -1
    xleg = xmin + 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)

    fig = plt.figure(figsize=(12,12))

    subplots = (331, 332, 333, 334, 335, 336, 337, 338)
    idx = (0, 1, 2, 3, 4, 5, 6, 7)
    bands = (12, 13, 14, 15, 16, 17, 18, 19)
    labels= ('WISE 2', 'WISE 3', 'WISE 4', 'P100', 'P160', 'S250', 'S350', 'S500')
    
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


        if idx == 4:
           file = obsdir+'/lf/L160microns.dat'
           lm,p,dpn,dpu = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
           dml       = 0.15
           xm        = [25.0,25.0]
           xm[1]     = xm[1] + dml
           tenpctocm = 3.086e19
           corrpc    = np.log10(tenpctocm)*2.0
           xobs      = -2.5*(lm[0:14]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
           xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
           dmm       = np.abs(xm[0] - xm[1])
           yobs      = np.log10(pow(10.0,p[0:14]) * dml/dmm)-3.0*np.log10(0.677)
           ydn       = np.log10(pow(10.0,dpn[0:14]) * dml/dmm)-3.0*np.log10(0.677)
           yup       = np.log10(pow(10.0,dpu[0:14]) * dml/dmm)-3.0*np.log10(0.677)
          
           ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*')
          
           dml       = 0.4
           xm        = [25.0,25.0]
           xm[1]     = xm[1] + dml
           tenpctocm = 3.086e19
           corrpc    = np.log10(tenpctocm)*2.0
           xobs      = -2.5*(lm[15:21]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
           xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
           dmm       = np.abs(xm[0] - xm[1])
           yobs      = np.log10(pow(10.0,p[15:21]) * dml/dmm)-3.0*np.log10(0.677)
           ydn       = np.log10(pow(10.0,dpn[15:21]) * dml/dmm)-3.0*np.log10(0.677)
           yup       = np.log10(pow(10.0,dpu[15:21]) * dml/dmm)-3.0*np.log10(0.677)
          
           ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='p')
 

        if idx == 5:
            file = obsdir+'/lf/lf250_dye10.data'
            lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
            lx        = np.log10(lm[0:7])
            dml       = np.abs(lx[1] - lx[0])
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(np.log10(lm[0:7])-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
            dmm       = np.abs(xobs[0] - xobs[1])
            yobs      = np.log10(p[0:7] * dml/dmm)-3.0*np.log10(0.677)
            ydn       = np.log10((p[0:7]- dp[0:7]) * dml/dmm)-3.0*np.log10(0.677)
            yup       = np.log10((p[0:7]+ dp[0:7]) * dml/dmm)-3.0*np.log10(0.677)
           
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')
           
            file = obsdir+'/lf/L250microns.dat'
            lm,p,dpn,dpu = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
            dml       = 0.15
            xm        = [25.0,25.0]
            xm[1]     = xm[1] + dml
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(lm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
            pxm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
            dmm       = np.abs(xm[0] - xm[1])
            yobs      = np.log10(pow(10.0,p) * dml/dmm)-3.0*np.log10(0.677)
            ydn       = np.log10(pow(10.0,dpn) * dml/dmm)-3.0*np.log10(0.677)
            yup       = np.log10(pow(10.0,dpu) * dml/dmm)-3.0*np.log10(0.677)
           
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*')

        if idx == 6:

            file = obsdir+'/lf/L350microns.dat'
            lm,p,dpn,dpu = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
            dml       = 0.2
            xm        = [25.0,25.0]
            xm[1]     = xm[1] + dml
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(lm[0:9]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
            xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
            dmm       = np.abs(xm[0] - xm[1])
            yobs      = np.log10(pow(10.0,p[0:9]) * dml/dmm)-3.0*np.log10(0.677)
            ydn       = np.log10(pow(10.0,dpn[0:9]) * dml/dmm)-3.0*np.log10(0.677)
            yup       = np.log10(pow(10.0,dpu[0:9]) * dml/dmm)-3.0*np.log10(0.677)
            
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*')
            
            dml       = 0.3
            xm        = [25.0,25.0]
            xm[1]     = xm[1] + dml
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(lm[10:19]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
            xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
            dmm       = np.abs(xm[0] - xm[1])
            yobs      = np.log10(pow(10.0,p[10:19]) * dml/dmm)-3.0*np.log10(0.677)
            ydn       = np.log10(pow(10.0,dpn[10:19]) * dml/dmm)-3.0*np.log10(0.677)
            yup       = np.log10(pow(10.0,dpu[10:19]) * dml/dmm)-3.0*np.log10(0.677)
            
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s')

        if idx == 7:
            file = obsdir+'/lf/L500microns.dat'
            lm,p,dpn,dpu = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
            dml       = 0.2
            xm        = [25.0,25.0]
            xm[1]     = xm[1] + dml
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(lm[0:11]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
            xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
            dmm       = np.abs(xm[0] - xm[1])
            yobs      = np.log10(pow(10.0,p[0:11]) * dml/dmm)-3.0*np.log10(0.677)
            ydn       = np.log10(pow(10.0,dpn[0:11]) * dml/dmm)-3.0*np.log10(0.677)
            yup       = np.log10(pow(10.0,dpu[0:11]) * dml/dmm)-3.0*np.log10(0.677)
            
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*',label="Marchetti+2016")
            
            dml       = 0.3
            xm        = [25.0,25.0]
            xm[1]     = xm[1] + dml
            tenpctocm = 3.086e19
            corrpc    = np.log10(tenpctocm)*2.0
            xobs      = -2.5*(lm[11:19]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
            xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
            dmm       = np.abs(xm[0] - xm[1])
            yobs      = np.log10(pow(10.0,p[11:19]) * dml/dmm)-3.0*np.log10(0.677)
            ydn       = np.log10(pow(10.0,dpn[11:19]) * dml/dmm)-3.0*np.log10(0.677)
            yup       = np.log10(pow(10.0,dpu[11:19]) * dml/dmm)-3.0*np.log10(0.677)

            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s',label="Negrello+2013")

            xobs = np.zeros(shape = 2)
            xobs[0] = 1.0
            xobs[1] = 1.0
            yobs=xobs
            ydn=xobs
            yup=xobs
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='p',label="Patel+2013")
            ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Dye+2010")

        if idx == 7:
            common.prepare_legend(ax, ['k','SlateGray','Gray','k','b','r','LightSalmon','grey','grey','grey','grey'], bbox_to_anchor=[1.1,-0.1])

    common.savefig(outdir, fig, "IR_luminosity_functions.pdf")

def plot_colours(plt, outdir, obsdir, h0, colours_dist):

    #plot colour distributions
    fig = plt.figure(figsize=(9.7,11.7))
    xtit = "$\\rm (u-r)$"
    ytit = "$\\rm dp/d(u-r)$"
    xmin, xmax, ymin, ymax = 0, 3.5, 0, 3
    xleg = xmax - 0.4 * (xmax - xmin)
    yleg = ymax - 0.15 * (ymax - ymin)
   
    subplots = (321, 322, 323, 324, 325, 326)
    indeces = (0, 1, 2, 3, 4, 5)
    labels  = ("(-17.13,-17.88)","(-17.88,-18.63)","(-18.63,-19.38)","(-19.38,-20.13)","(-20.13,-20.88)","(-20.88,-21.63)")
   
    file = obsdir+'/Colours/ur_colours_SDSS.data'
    cSDSS,dpSDSS = np.loadtxt(file,usecols=[0,1],unpack=True)
    columns = [0,38,77,112,151,190,229]
   
    for subplot, idx in zip(subplots, indeces):
   
        ax = fig.add_subplot(subplot)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
        ax.text(xleg, yleg, labels[idx])
        
        if(idx == 0):
           c1 = columns[idx]
        if(idx > 0):
           c1 = columns[idx]+1
        # Observed CDF
        yzeros = np.zeros(shape = len(cSDSS[c1:columns[idx+1]]))
        ax.plot(cSDSS[c1:columns[idx+1]], dpSDSS[c1:columns[idx+1]],'grey', label='SDSS')
        ax.fill_between(cSDSS[c1:columns[idx+1]], dpSDSS[c1:columns[idx+1]], yzeros, facecolor='grey', alpha=1,interpolate=True)
        
        # Predicted CDF
        y = colours_dist[0,idx,0,:]
        ind = np.where(y > 0.)
        ax.plot(xc[ind],y[ind],'r', label='Shark')
        
        common.prepare_legend(ax, ['grey','r'], loc = 2)

    common.savefig(outdir, fig, 'ur_colour_z0.pdf')

    fig = plt.figure(figsize=(9.7,11.7))
    xtit = "$\\rm (g-r)$"
    ytit = "$\\rm dp/d(g-r)$"
    xmin, xmax, ymin, ymax = 0, 1.1, 0, 6.5
    xleg = xmax - 0.4 * (xmax - xmin)
    yleg = ymax - 0.15 * (ymax - ymin)
    
    file = obsdir+'/Colours/gr_colours_SDSS.data'
    cSDSS,dpSDSS = np.loadtxt(file,usecols=[0,1],unpack=True)
    columns = [0,21,43,65,87,109,131]
    
    for subplot, idx in zip(subplots, indeces):
    
        ax = fig.add_subplot(subplot)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
        ax.text(xleg, yleg, labels[idx])
    
        if(idx == 0):
           c1 = columns[idx]
        if(idx > 0):
           c1 = columns[idx]+1
        # Observed CDF
        yzeros = np.zeros(shape = len(cSDSS[c1:columns[idx+1]]))
        ax.plot(cSDSS[c1:columns[idx+1]], dpSDSS[c1:columns[idx+1]],'grey', label='SDSS')
        ax.fill_between(cSDSS[c1:columns[idx+1]], dpSDSS[c1:columns[idx+1]], yzeros, facecolor='grey', alpha=1,interpolate=True)
    
        # Predicted CDF
        y = colours_dist[0,idx,1,:]
        ind = np.where(y > 0.)
        ax.plot(xc[ind],y[ind],'r', label='Shark')
    
        common.prepare_legend(ax, ['grey','r'], loc = 2)
    
    common.savefig(outdir, fig, 'gr_colour_z0.pdf')
    
def prepare_data(hdf5_data, phot_data, LFs_dust, LFs_nodust, colours_dist, index, nbands):
    
    (h0, _, mdisk, mbulge, mhalo, mshalo, typeg, age, 
     sfr_disk, sfr_burst, id_gal) = hdf5_data
    
    mass   = np.zeros(shape = len(mhalo))
    masssh = np.zeros(shape = len(mhalo))

    ind = np.where((typeg <= 0) & (mdisk+mbulge > 1e5))
    mass[ind] = np.log10(mhalo[ind]) - np.log10(float(h0))
    ind = np.where((typeg <= 1) & (mdisk+mbulge > 1e5))
    masssh[ind] = np.log10(mshalo[ind]) - np.log10(float(h0))

    #components:
    #(len(my_data), 2, 2, 5, nbands)
    #0: disk instability bulge
    #1: galaxy merger bulge
    #2: total bulge
    #3: disk
    #4: total
    SEDs_dust = phot_data[:,1,0,:,:]
    SEDs_nodust = phot_data[:,0,0,:,:]
    
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


    uband = 2
    gband = 3
    rband = 4
    ubandl = SEDs_dust[:,4,uband]
    gbandl = SEDs_dust[:,4,gband]
    rbandl = SEDs_dust[:,4,rband]

    for mag in range(0,len(magbins)-1):
        ind = np.where((ubandl < -1) & (gbandl < -1) & (rbandl < magbins[mag]) & (rbandl >= magbins[mag+1]))
        H, bins_edges  = np.histogram(ubandl[ind] - rbandl[ind],bins=np.append(cbins,cupp))
        colours_dist[index,mag,0,:] = colours_dist[index,mag,0,:] + H
        colours_dist[index,mag,0,:] = colours_dist[index,mag,0,:] / (len(ubandl[ind]) * dc)
        H, bins_edges  = np.histogram(gbandl[ind] - rbandl[ind],bins=np.append(cbins,cupp))
        colours_dist[index,mag,1,:] = colours_dist[index,mag,1,:] + H
        colours_dist[index,mag,1,:] = colours_dist[index,mag,1,:] / (len(gbandl[ind]) * dc)

def main(model_dir, outdir, redshift_table, subvols, obsdir):

    # Loop over redshift and subvolumes
    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mvir_hosthalo',
                           'mvir_subhalo', 'type', 'mean_stellar_age', 
                           'sfr_disk', 'sfr_burst', 'id_galaxy')}

    z = (0, 0.25, 0.5)
    snapshots = redshift_table[z]

    # Create histogram
    for index, snapshot in enumerate(snapshots):

        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        seds, ids, nbands = common.read_photometry_data(model_dir, snapshot, subvols)
        
        if(index == 0):
            LFs_dust     = np.zeros(shape = (len(z), 5, nbands, len(mbins)))
            LFs_nodust   = np.zeros(shape = (len(z), 5, nbands, len(mbins)))
            colours_dist = np.zeros(shape = (len(z), len(magbins)-1, 2, len(cbins)))

        prepare_data(hdf5_data, seds, LFs_dust, LFs_nodust, colours_dist, index, nbands)

        h0, volh = hdf5_data[0], hdf5_data[1]
        if(volh > 0.):
            LFs_dust[index,:]   = LFs_dust[index,:]/volh
            LFs_nodust[index,:] = LFs_nodust[index,:]/volh

    # Take logs
    ind = np.where(LFs_dust > 0.)
    LFs_dust[ind] = np.log10(LFs_dust[ind])

    ind = np.where(LFs_nodust > 0.)
    LFs_nodust[ind] = np.log10(LFs_nodust[ind])

    plot_lfs(plt, outdir, obsdir, h0, LFs_dust, LFs_nodust)
    plot_colours(plt, outdir, obsdir, h0, colours_dist)
   
if __name__ == '__main__':
    main(*common.parse_args())
