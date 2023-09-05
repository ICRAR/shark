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

import functools

import numpy as np
import os
import h5py

import common
import utilities_statistics as us

##################################
# Constants
mlow = 0.0
mupp = 12.5
dm = 0.2
mbins = np.arange(mlow, mupp, dm)
xmf = mbins + dm/2.0

tlow = 0.0
tupp = 5
dt = 0.15
tbins = np.arange(tlow, tupp, dt)
xtf = tbins + dt/2.0


zsun = 0.02 #189
MpctoKpc = 1e3

#model of Mattson et al. (2014) for the dependence of the dust-to-metal mass ratio and metallicity X/H.
corrfactor_dm = 2.0
polyfit_dm = [ 0.00544948, 0.00356938, -0.07893235,  0.05204814,  0.49353238] 

#choose dust model between mm14, rr14 and constdust
m14 = False
rr14 = True
constdust = False
rr14xcoc = False

#read EAGLE tables
sdust_eaglet, taumed_eagle, taulow_eagle, tauhigh_eagle = common.load_observation('../data', 'Models/EAGLE/Tau5500-Trayford-EAGLE.dat', [0,1,2,3])
sdust_eaglem, mmed_eagle, mlow_eagle, mhigh_eagle = common.load_observation('../data/','Models/EAGLE/CFPowerLaw-Trayford-EAGLE.dat', [0,1,2,3])

def interp (sdust_eagle, med_eagle, low_eagle, high_eagle):

    med = np.zeros(shape = (2,len(sdust_eagle)-1))
    low = np.zeros(shape = (2,len(sdust_eagle)-1))
    hig = np.zeros(shape = (2,len(sdust_eagle)-1))
   
    for i in range(0,len(sdust_eagle)-1):
        delta_dust = sdust_eagle[i+1] - sdust_eagle[i]
        med[0,i] = (med_eagle[i+1] - med_eagle[i] ) / delta_dust
        med[1,i] =  med_eagle[i+1]- med[0,i] * sdust_eagle[i+1]
        low[0,i] = (low_eagle[i+1] - low_eagle[i]) / delta_dust
        low[1,i] =  low_eagle[i+1]- low[0,i] * sdust_eagle[i+1]
        hig[0,i] = (high_eagle[i+1] - high_eagle[i]) / delta_dust
        hig[1,i] =  high_eagle[i+1]- hig[0,i] * sdust_eagle[i+1]

    return (med, low, hig)

(m_med, m_low, m_hig) = interp (sdust_eaglet, taumed_eagle, taulow_eagle, tauhigh_eagle)
(s_med, s_low, s_hig) = interp (sdust_eaglem, mmed_eagle, mlow_eagle, mhigh_eagle)

def dust_mass(mz, mg, h0):
    md = np.zeros(shape = len(mz))
    ind = np.where((mz > 0) & (mg > 0))
    XHd = np.log10(mz[ind]/mg[ind]/zsun)
    if(m14 == True):
        DToM = (polyfit_dm[0] * XHd**4.0 + polyfit_dm[1] * XHd**3.0 + polyfit_dm[2] * XHd**2.0 + polyfit_dm[3] * XHd + polyfit_dm[4])/corrfactor_dm
        DToM = np.clip(DToM, 1e-6, 0.5)
        md[ind] = mz[ind]/h0 * DToM
        DToM_MW = polyfit_dm[4]/corrfactor_dm
    elif(rr14 == True):
         y = np.zeros(shape = len(XHd))
         highm = np.where(XHd > -0.59)
         y[highm] = 10.0**(2.21 - XHd[highm]) #gas-to-dust mass ratio
         lowm = np.where(XHd <= -0.59)
         y[lowm] = 10.0**(0.26 - (3.1+1.33) * XHd[lowm]) #gas-to-dust mass ratio
         DToM = 1.0 / y / (mz[ind]/mg[ind])
         DToM = np.clip(DToM, 1e-6, 1)
         md[ind] = mz[ind]/h0 * DToM
         DToM_MW = 1.0 / (10.0**(2.21)) / zsun
    elif(rr14xcoc == True):
        y = np.zeros(shape = len(XHd))
        highm = np.where(XHd > -0.15999999999999998)
        y[highm] = 10.0**(2.21 - XHd[highm]) #gas-to-dust mass ratio
        lowm = np.where(XHd <= -0.15999999999999998)
        y[lowm] = 10.0**(1.66 - 4.43 * XHd[lowm]) #gas-to-dust mass ratio
        DToM = 1.0 / y / (mz[ind]/mg[ind])
        DToM = np.clip(DToM, 1e-6, 1)
        md[ind] = mz[ind]/h0 * DToM
        DToM_MW = 1.0 / (10.0**(2.21)) / zsun
    elif(constdust == True):
         md[ind] = 0.33 * mz[ind]/h0
         DToM_MW = 0.33
   
    return (md, DToM_MW)

# define tau diffuse
def tau_diff (md, rd, hd, h0):

    tau    = np.zeros(shape = len(md))
    sigma  = np.zeros(shape = len(md))

    ind = np.where((md > 0) & (rd > 0) & (hd > 0))
    sigma[ind] = np.log10(md[ind]/h0 / (2.0 * 3.1416 * rd[ind]*MpctoKpc/h0 * hd[ind]*MpctoKpc/h0)) #in Msun/kpc^2
    # cap surface density of dust to physical values based on EAGLE
    sigma[ind] = np.clip(sigma[ind],0, 12)

    #interpolate linearly in dust surface density going through the list of values in the EAGLE table
    for i in range(0,len(sdust_eaglet)-1):
        if(i == 0):
           selecinrage = np.where(sigma <= sdust_eaglet[i+1])
        elif(i == len(sdust_eaglet)-2):
           selecinrage = np.where(sigma >= sdust_eaglet[i])
        else:
           selecinrage = np.where((sigma >= sdust_eaglet[i]) & (sigma < sdust_eaglet[i+1]))

        tau[selecinrage] = m_med[0,i] * sigma[selecinrage] + m_med[1,i]
        tlow = abs(tau[selecinrage] - (m_low[0,i] * sigma[selecinrage] + m_low[1,i]))
        thigh= abs( (m_hig[0,i] * sigma[selecinrage] + m_hig[1,i]) - tau[selecinrage])
        var_gauss = (tlow + thigh) * 0.5
        pert = np.random.randn(len(var_gauss)) * np.sqrt(var_gauss)
        tau[selecinrage] = tau[selecinrage] + pert

    # cap it to maximum and minimum values in EAGLE
    tau = np.clip(tau, 1e-6, 5) 

    return (tau, sigma) 

# define clump tau
def tau_clump(mz,mg, h0, sigmag, tau_diff):
    sigmaclump = np.zeros(shape = len(mg))
    sigmaclump[:] = 85.0*1e6 #in Msun/kpc^3
    ind = np.where(sigmaclump < sigmag)
    sigmaclump[ind] = sigmag[ind]
    tau = np.zeros(shape = len(mz))
    ind = np.where((mz > 0) & (mg > 0))
    (md, DToM_MW)  = dust_mass(mz[ind],mg[ind],h0)
    norm = 85.0*1e6 * DToM_MW * zsun #dust surface density of clumps
    tau[ind] = 0.5 * (sigmaclump[ind] * md/(mg[ind]/h0) / norm)
    ind = np.where(tau < tau_diff)
    tau[ind] = tau_diff[ind]
    # cap it to maximum and minimum values in EAGLE but also forcing the clump tau to be at least as high as the diffuse tau
    tau = np.clip(tau, 1e-6, 5)
    return tau

def tau_clump2(mz,mg,h0):
    tau = np.zeros(shape = len(mz))
    ind = np.where((mz > 0) & (mg > 0))
    (md, DToM_MW)  = dust_mass(mz[ind],mg[ind],h0)
    tau[ind] = 1.5 * (md/mz[ind]/DToM_MW)
    # cap it to maximum and minimum values in EAGLE but also forcing the clump tau to be at least as high as the diffuse tau
    tau = np.clip(tau, 1e-6, 5)
    return tau


# define slope diffuse
def slope_diff (md, rd, hd, h0):

    m    = np.zeros(shape = len(md))
    sigma  = np.zeros(shape = len(md))

    ind = np.where((md > 0) & (rd > 0) & (hd > 0))
    sigma[ind] = np.log10(md[ind]/h0 / (2.0 * 3.1416 * rd[ind]*MpctoKpc/h0 * hd[ind]*MpctoKpc/h0)) #in Msun/kpc^2
    # cap surface density of dust to physical values based on EAGLE
    sigma[ind] = np.clip(sigma[ind],0, 12)

    #interpolate linearly in dust surface density going through the list of values in the EAGLE table
    for i in range(0,len(sdust_eaglem)-1):
        if(i == 0):
           selecinrage = np.where(sigma <= sdust_eaglem[i+1])
        elif(i == len(sdust_eaglet)-2):
           selecinrage = np.where(sigma >= sdust_eaglem[i])
        else:
           selecinrage = np.where((sigma >= sdust_eaglem[i]) & (sigma < sdust_eaglem[i+1]))

        m[selecinrage] = s_med[0,i] * sigma[selecinrage] + s_med[1,i]
        tlow = abs(m[selecinrage] - (s_low[0,i] * sigma[selecinrage] + s_low[1,i]))
        thigh= abs((s_hig[0,i] * sigma[selecinrage] + s_hig[1,i]) - m[selecinrage])
        var_gauss = (tlow + thigh) * 0.5
        pert = np.random.randn(len(var_gauss)) * np.sqrt(var_gauss)
        m[selecinrage] = m[selecinrage] + pert

    # cap it to maximum and minimum values in EAGLE
    m = np.clip(m, -3, -0.001)
    return m 

def prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit):
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)
    xleg = xmax - 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)
    #ax.text(xleg, yleg, 'z=0')

def prepare_data(hdf5_data, index, tdiff, tcloud, sigmad_diff, sigmag_diff, sfr_rat, met_evo,  m_diff, hist_sigmad, tau_comp, 
                 model_dir, snapshot, subvol, writeon):

    bin_it = functools.partial(us.wmedians, xbins=xmf)
    bin_it_tau = functools.partial(us.wmedians, xbins=xtf)

    # Unpack data
    (h0, volh, typeg, rgasd, rgasb, mHId, mH2d, mgasd, mHIb, mH2b, mgasb, mzd, mzb, mdisk, mbulge, sfrd, sfrb, idgal) = hdf5_data
    XH = 0.72
    h0log = np.log10(float(h0))

    sigma_g_d = mgasd/h0/(2.0 * 3.1416 * (rgasd/h0*1e3)**2.0)
    sigma_g_b = mgasb/h0/(2.0 * 3.1416 * (rgasb/h0*1e3)**2.0)

    (mdustd, DToM_MW) = dust_mass(mzd, mgasd, h0)
    (mdustb, DToM_MW) = dust_mass(mzb, mgasb, h0)

    #random numbers from 0 to 1
    sin_inclination = np.random.rand(len(mdustd)) 
    inclination = np.arcsin(sin_inclination) * 180.0/3.1416 #degrees
    bd  = sin_inclination*(rgasd - rgasd/7.3)  +  rgasd/7.3 #scaleheight at r50

    (tau_dust_bulge, sigmab) = tau_diff(mdustb, rgasb, rgasb, h0)
    (tau_dust_disk, sigmad) = tau_diff(mdustd, rgasd, bd, h0)
    
    tau_clump_bulge = tau_clump(mzb, mgasb, h0, sigma_g_b, tau_dust_bulge)
    tau_clump_disk  = tau_clump(mzd, mgasd, h0, sigma_g_d, tau_dust_disk)
    slope_dust_bulge = slope_diff(mdustb, rgasb, rgasb, h0) 
    slope_dust_disk  = slope_diff(mdustd, rgasd, bd, h0)

    tau_clump_bulge2 = tau_clump2(mzb, mgasb, h0)
    tau_clump_disk2  = tau_clump2(mzd, mgasd, h0)

    ind = np.where((tau_clump_disk > 0) & (tau_clump_disk2 > 0))
    tau_comp[index,0,:] = bin_it_tau(x=tau_clump_disk2[ind], y=tau_clump_disk[ind])
    ind = np.where((tau_clump_bulge > 0) & (tau_clump_bulge2 > 0))
    tau_comp[index,1,:] = bin_it_tau(x=tau_clump_bulge2[ind], y=tau_clump_bulge[ind])

    mass = np.log10((mdisk +  mbulge)/h0)
    #ignore these low mass satellites because they are not converged
    lowmasssat = np.where((mass <= 8.3) & (sfrd+sfrb > 0) & (typeg > 0))
    sfrd[lowmasssat] = 0
    sfrb[lowmasssat] = 0
 
    ind = np.where((mass >= 6) & (sfrd > 0))
    tdiff[index,0,:]  = bin_it(x=mass[ind], y=tau_dust_disk[ind])
    tcloud[index,0,:] = bin_it(x=mass[ind], y=tau_clump_disk[ind])
    m_diff[index,0,:]  = bin_it(x=mass[ind], y=slope_dust_disk[ind])
    sigmad_diff[index,0,:]  = bin_it(x=mass[ind], y=sigmad[ind])
    sigmag_diff[index,0,:]  = bin_it(x=mass[ind], y=np.log10(sigma_g_d[ind]))
    H, _ = np.histogram(sigmad[ind],bins=np.append(mbins,mupp))
    hist_sigmad[0,index,:] = hist_sigmad[0,index,:] + H
    met_evo[index,0,:] = bin_it(x=mass[ind], y=np.log10((mzd[ind])/(mgasd[ind])/zsun))

    ind = np.where((mass >= 6) & (sfrb > 0))
    tdiff[index,1,:]  = bin_it(x=mass[ind], y=tau_dust_bulge[ind])
    tcloud[index,1,:] = bin_it(x=mass[ind], y=tau_clump_bulge[ind])
    sigmad_diff[index,1,:]  = bin_it(x=mass[ind], y=sigmab[ind])
    m_diff[index,1,:]  = bin_it(x=mass[ind], y=slope_dust_bulge[ind])
    sigmag_diff[index,1,:]  = bin_it(x=mass[ind], y=np.log10(sigma_g_b[ind]))
    H, _ = np.histogram(sigmad[ind],bins=np.append(mbins,mupp))
    hist_sigmad[1,index,:] = hist_sigmad[1,index,:] + H
    #divide by volume and interval
    hist_sigmad[0,index,:] = hist_sigmad[0,index,:] / volh / dm
    hist_sigmad[1,index,:] = hist_sigmad[1,index,:] / volh / dm
    ind = np.where((mass >= 6) & (mgasb > 0))
    met_evo[index,1,:] = bin_it(x=mass[ind], y=np.log10((mzb[ind])/(mgasb[ind])/zsun))

    ind = np.where((mass >= 6) & (sfrd + sfrb > 0))
    sfr_rat[index,:] = bin_it(x=mass[ind], y=sfrb[ind]/(sfrd[ind]+sfrb[ind]))

def plot_taus(plt, output_dir, tdiff, tcloud, sigmad_diff, sigmag_diff, sfr_rat, met_evo, m_diff, 
              hist_sigmad, tau_comp, zlist):

    #tau diffuse medium
    fig = plt.figure(figsize=(14,4.5))

    xmin, xmax, ymin, ymax = 6.0, 12.0, 0.0, 2.5
    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\tau$"

    xleg = xmax - 0.6 * (xmax - xmin)
    yleg = ymax + 0.02 * (ymax - ymin)

    colors  = ('DarkRed','Salmon','Orange','YellowGreen','MediumSeaGreen','DarkTurquoise','MediumBlue','Purple')
    subplots = (141, 142, 143, 144)
    labels = ('disk ISM','bulge ISM','disk BC','bulge BC')

    for s in range(0,len(subplots)):
        ax = fig.add_subplot(subplots[s])
        prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)
        ax.text(xleg, yleg, labels[s])

        if(s <= 1):
           j = s 
           for i in range(0,len(tcloud[:,j,0,0])):
               # Predicted relation
               ind = np.where(tdiff[i,j,0,:] > 0)
               xplot = xmf[ind]
               yplot = tdiff[i,j,0,ind]
               errdn = tdiff[i,j,1,ind]
               errup = tdiff[i,j,2,ind]
               ax.plot(xplot,yplot[0],color=colors[i])
               ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=colors[i], alpha=0.45,interpolate=True)
               ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=colors[i], alpha=0.45,interpolate=True)
        else:
           j = s - 2
           for i in range(0,len(tcloud[:,j,0,0])):
               # Predicted relation
               ind = np.where(tcloud[i,j,0,:] > 0)
               xplot = xmf[ind]
               yplot = tcloud[i,j,0,ind]
               errdn = tcloud[i,j,1,ind]
               errup = tcloud[i,j,2,ind]
               ax.plot(xplot,yplot[0],color=colors[i],label='z=%s' % str(zlist[i]))
               ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=colors[i], alpha=0.45,interpolate=True)
               ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=colors[i], alpha=0.45,interpolate=True)
        if(s <= 1):
           x=[6,12]
           y=[0.3,0.3]
           ax.plot(x,y,linestyle='solid',color='k')
        else:
           x=[6,12]
           y=[1,1]
           ax.plot(x,y,linestyle='solid',color='k')
        if(s == 2):
           common.prepare_legend(ax, colors, loc='upper left')
    common.savefig(output_dir, fig, "extinction_EAGLE_predictions.pdf")

    #tau comparison 
    fig = plt.figure(figsize=(7,4.5))

    xmin, xmax, ymin, ymax = 0,5, 0.0, 5
    xtit="$\\tau_{\\rm clump,simple}$"
    ytit="$\\tau_{\\rm clump, complex}$"

    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (121, 122)
    labels = ('disk','bulge')

    for j in range(0,len(subplots)):
        ax = fig.add_subplot(subplots[j])
        prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)
        ax.text(xleg, yleg, labels[j])

        for i in range(0,len(tau_comp[:,j,0,0])):
            # Predicted relation
            ind = np.where(tau_comp[i,j,0,:] > 0)
            xplot = xtf[ind]
            yplot = tau_comp[i,j,0,ind]
            errdn = tau_comp[i,j,1,ind]
            errup = tau_comp[i,j,2,ind]
            ax.plot(xplot,yplot[0],color=colors[i],label='z=%s' % str(zlist[i]))
            #ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=colors[i], alpha=0.5,interpolate=True)
            #ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=colors[i], alpha=0.5,interpolate=True)

        x=[0,5]
        y=[0,5]
        ax.plot(x,y,linestyle='solid',color='k')
    common.prepare_legend(ax, colors, loc='upper left')
    common.savefig(output_dir, fig, "tau_clumps_comparison.pdf")

    #tau clumps 
    fig = plt.figure(figsize=(7,4.5))

    xmin, xmax, ymin, ymax = 6.0, 12.0, 0.0, 2.5
    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\tau_{\\rm clump}$"

    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (121, 122)
    labels = ('disk','bulge')

    for j in range(0,len(subplots)):
        ax = fig.add_subplot(subplots[j])
        prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)
        ax.text(xleg, yleg, labels[j])

        for i in range(0,len(tcloud[:,j,0,0])):
            # Predicted relation
            ind = np.where(tcloud[i,j,0,:] > 0)
            xplot = xmf[ind]
            yplot = tcloud[i,j,0,ind]
            errdn = tcloud[i,j,1,ind]
            errup = tcloud[i,j,2,ind]
            ax.plot(xplot,yplot[0],color=colors[i],label='z=%s' % str(zlist[i]))
            ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=colors[i], alpha=0.5,interpolate=True)
            ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=colors[i], alpha=0.5,interpolate=True)

        x=[6,12]
        y=[1.5,1.5]
        ax.plot(x,y,linestyle='solid',color='k')
    common.prepare_legend(ax, colors, loc='upper left')
    common.savefig(output_dir, fig, "extinction_clumps_EAGLE_predictions.pdf")

    #slope diffuse 
    fig = plt.figure(figsize=(8,5.5))

    xmin, xmax, ymin, ymax = 6.0, 12.0, -3.2, 0
    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\eta_{\\rm ISM}$"

    xleg = xmax - 0.55 * (xmax - xmin)
    yleg = ymax + 0.02 * (ymax - ymin)

    subplots = (121, 122)
    labels = ('disk','bulge')

    for j in range(0,len(subplots)):
        ax = fig.add_subplot(subplots[j])
        prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)
        ax.text(xleg, yleg, labels[j],fontsize=12)

        for i in range(0,len(m_diff[:,j,0,0])):
            # Predicted relation
            ind = np.where(m_diff[i,j,0,:] != 0)
            xplot = xmf[ind]
            yplot = m_diff[i,j,0,ind]
            errdn = m_diff[i,j,1,ind]
            errup = m_diff[i,j,2,ind]
            ax.plot(xplot,yplot[0],color=colors[i],label='z=%s' % str(zlist[i]))
            ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=colors[i], alpha=0.45,interpolate=True)
            ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=colors[i], alpha=0.45,interpolate=True)
        if j == 0:
           common.prepare_legend(ax, colors, loc='lower right')

        x=[6,12]
        y=[-0.7,-0.7]
        ax.plot(x,y,linestyle='solid',color='k')
    common.savefig(output_dir, fig, "slope_extinction_diffuse_EAGLE_predictions.pdf")

    #distribution of surface densities of dust
    fig = plt.figure(figsize=(7,4.5))

    xmin, xmax, ymin, ymax = 0.3, 10.5, -6, -1
    xtit="$\\rm log_{10} (\\Sigma_{\\rm dust}/M_{\odot} kpc^{-2})$"
    ytit="$\\rm log_{10}(\Phi/dlog_{10}({\\Sigma_{\\rm dust}})/{\\rm Mpc}^{-3} )$"

    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (211, 212)
    labels = ('disk','bulge')

    for j in range(0,len(subplots)):
        ax = fig.add_subplot(subplots[j])
        if(j == 1):
           ytitle = ''
        else:
           ytitle = ytit
        prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytitle)

        ax.text(xleg, yleg, labels[j])
        for i in range(0,len(zlist)):
            # Predicted relation
            ind = np.where(hist_sigmad[j,i,:] != 0)
            xplot = xmf[ind]
            yplot = hist_sigmad[j,i,ind]
            errdn = hist_sigmad[j,i,ind]
            errup = hist_sigmad[j,i,ind]
            ax.plot(xplot,yplot[0],color=colors[i],label='z=%s' % str(zlist[i]))

    #common.prepare_legend(ax, colors, loc='upper left')
    common.savefig(output_dir, fig, "sigma_dust_predictions_distribution.pdf")

    #surface densities of dust
    fig = plt.figure(figsize=(7,4.5))

    xmin, xmax, ymin, ymax = 6, 12.0, 0, 9
    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm log_{10} (\\Sigma_{\\rm dust}/M_{\odot} kpc^{-2})$"

    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymin + 0.1 * (ymax - ymin)

    subplots = (121, 122)
    labels = ('disk','bulge')

    for j in range(0,len(subplots)):
        ax = fig.add_subplot(subplots[j])
        if(j == 1):
           ytitle = ''
        else:
           ytitle=ytit
        prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytitle)
        ax.text(xleg, yleg, labels[j])

        for i in range(0,len(sigmad_diff[:,j,0,0])):
            # Predicted relation
            ind = np.where(sigmad_diff[i,j,0,:] > 0)
            xplot = xmf[ind]
            yplot = sigmad_diff[i,j,0,ind]
            errdn = sigmad_diff[i,j,1,ind]
            errup = sigmad_diff[i,j,2,ind]
            ax.plot(xplot,yplot[0],color=colors[i],label='z=%s' % str(zlist[i]))
            ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=colors[i], alpha=0.5,interpolate=True)
            ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=colors[i], alpha=0.5,interpolate=True)
        if(j == 0):
           common.prepare_legend(ax, colors, loc='upper left')
    common.savefig(output_dir, fig, "sigma_dust_predictions.pdf")

    #surface densities of gas
    fig = plt.figure(figsize=(7,4.5))

    xmin, xmax, ymin, ymax = 6.0, 12.0, 4, 11.1
    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm log_{10} (\\Sigma_{\\rm gas}/M_{\odot} kpc^{-2})$"

    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (121, 122)
    labels = ('disk','bulge')

    for j in range(0,len(subplots)):
        ax = fig.add_subplot(subplots[j])
        prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)
        ax.text(xleg, yleg, labels[j])

        for i in range(0,len(sigmag_diff[:,j,0,0])):
            # Predicted relation
            ind = np.where(sigmag_diff[i,j,0,:] > 0)
            xplot = xmf[ind]
            yplot = sigmag_diff[i,j,0,ind]
            errdn = sigmag_diff[i,j,1,ind]
            errup = sigmag_diff[i,j,2,ind]
            ax.plot(xplot,yplot[0],color=colors[i],label='z=%s' % str(zlist[i]))
            ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=colors[i], alpha=0.5,interpolate=True)
            ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=colors[i], alpha=0.5,interpolate=True)

    common.prepare_legend(ax, colors, loc='upper left')
    common.savefig(output_dir, fig, "sigma_gas_predictions.pdf")

    #SFR ratio
    fig = plt.figure(figsize=(4.5,4.5))

    xmin, xmax, ymin, ymax = 6.0, 12.0, 0, 1
    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm SFR_{\\rm burst}/SFR_{\\rm tot}$"

    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(111)
    prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)

    for i in range(0,len(sfr_rat[:,0,0])):
        # Predicted relation
        ind = np.where(sfr_rat[i,0,:] > 0)
        xplot = xmf[ind]
        yplot = sfr_rat[i,0,ind]
        errdn = sfr_rat[i,1,ind]
        errup = sfr_rat[i,2,ind]
        ax.plot(xplot,yplot[0],color=colors[i],label='z=%s' % str(zlist[i]))
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=colors[i], alpha=0.5,interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=colors[i], alpha=0.5,interpolate=True)

    common.prepare_legend(ax, colors, loc='upper left')
    common.savefig(output_dir, fig, "sfr_ratio_predictions.pdf")

    #metallicity evolution
    fig = plt.figure(figsize=(7,4.5))

    xmin, xmax, ymin, ymax = 6.0, 12.0, -2, 0.8
    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm log_{10} (Z_{\\rm gas}/Z_{\\odot})$"

    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (121, 122)
    labels = ('disk','bulge')

    for j in range(0,len(subplots)):
        ax = fig.add_subplot(subplots[j])
        prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)
        ax.text(xleg, yleg, labels[j])

        for i in range(0,len(met_evo[:,j,0,0])):
            # Predicted relation
            ind = np.where(met_evo[i,j,0,:] != 0)
            xplot = xmf[ind]
            yplot = met_evo[i,j,0,ind]
            errdn = met_evo[i,j,1,ind]
            errup = met_evo[i,j,2,ind]
            ax.plot(xplot,yplot[0],color=colors[i],label='z=%s' % str(zlist[i]))
            ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=colors[i], alpha=0.5,interpolate=True)
            ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=colors[i], alpha=0.5,interpolate=True)

    common.prepare_legend(ax, colors, loc='upper left')
    common.savefig(output_dir, fig, "metallicity_evo_predictions.pdf")

def main(model_dir, output_dir, redshift_table, subvols, obs_dir):

    #zlist = (0, 0.25, 0.5, 1, 2, 3, 4, 6, 8)
    zlist = (0, 0.5, 1, 2, 3, 4, 6, 8)

    plt = common.load_matplotlib()
    fields = {'galaxies': ('type', 'rgas_disk', 'rgas_bulge', 'matom_disk', 'mmol_disk', 'mgas_disk',
                           'matom_bulge', 'mmol_bulge', 'mgas_bulge', 'mgas_metals_disk', 
                           'mgas_metals_bulge', 'mstars_disk', 'mstars_bulge','sfr_disk','sfr_burst','id_galaxy')}

    tau_diff = np.zeros(shape = (len(zlist), 2, 3, len(xmf)))
    tau_cloud = np.zeros(shape = (len(zlist), 2, 3, len(xmf)))
    sigmad_diff = np.zeros(shape = (len(zlist), 2, 3, len(xmf)))
    sigmag_diff = np.zeros(shape = (len(zlist), 2, 3, len(xmf)))
    hist_sigmad       = np.zeros(shape = (2, len(zlist), len(xmf)))

    m_diff = np.zeros(shape = (len(zlist), 2, 3, len(xmf)))
    tau_comp = np.zeros(shape = (len(zlist), 2, 3, len(xtf)))

    sfr_rat = np.zeros(shape = (len(zlist), 3, len(xmf)))
    met_evo = np.zeros(shape = (len(zlist), 2, 3, len(xmf)))
    
    writeon = False

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        prepare_data(hdf5_data, index, tau_diff, tau_cloud, sigmad_diff, sigmag_diff, sfr_rat, 
                     met_evo, m_diff, hist_sigmad, tau_comp, model_dir, snapshot, subvols, writeon)

    ind = np.where(hist_sigmad > 0.)
    hist_sigmad[ind] = np.log10(hist_sigmad[ind])

    output_dir = os.path.join(output_dir, 'eagle-rr14')

    plot_taus(plt, output_dir, tau_diff, tau_cloud, sigmad_diff, sigmag_diff, sfr_rat, met_evo, 
              m_diff, hist_sigmad, tau_comp, zlist)

if __name__ == '__main__':
    main(*common.parse_args())
