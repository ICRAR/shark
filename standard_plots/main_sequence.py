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
MpcToKpc = 1e3

##################################
# Mass function initialization
mlow = 5
mupp = 14
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0
imf   = 'cha'

mlowlr  = 9.0
mupplr  = 11.4
dmlr    = 0.4
mbinslr = np.arange(mlowlr,mupplr,dmlr)
xmflr   = mbinslr + dmlr/2.0

mlow2 = 5.0
mupp2 = 15.0
dm2 = 0.4
mbins2 = np.arange(mlow2,mupp2,dm2)
xmf2 = mbins2 + dm2/2.0

ssfrlow = -6
ssfrupp = 4
dssfr = 0.2
ssfrbins = np.arange(ssfrlow,ssfrupp,dssfr)
xssfr    = ssfrbins + dssfr/2.0

def plot_sfr_mstars_z0(plt, outdir, obsdir, h0, sfr_seq, mainseqsf, sigmamainseqsf, slope_ms_z0, offset_ms_z0):

    fig = plt.figure(figsize=(5,5))
    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm log_{10}(\\rm SFR/M_{\odot} yr^{-1})$"

    xmin, xmax, ymin, ymax = 8, 12, -5, 3
    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #predicted relation
    ind = np.where((sfr_seq[0,:] > 0) & (sfr_seq[1,:] != 0) )
    xdata = sfr_seq[0,ind]
    ydata = sfr_seq[1,ind]
    us.density_contour(ax, xdata[0], ydata[0], 30, 30) #, **contour_kwargs)

    ind = np.where(mainseqsf[0,0,:] != 0)
    xplot = xmf[ind]
    yplot = mainseqsf[0,0,ind]
    ax.plot(xplot,yplot[0],color='k',linestyle='solid', linewidth = 1, label="Shark all galaxies")

    #SFR relation z=0
    lm, SFR = common.load_observation(obsdir, 'SFR/Brinchmann04.dat', (0, 1))
    hobs = 0.7
    #add cosmology correction plus IMF correction that goes into the stellar mass.
    corr_cos = np.log10(pow(hobs,2)/pow(h0,2)) - 0.09
    # apply correction to both stellar mass and SFRs.
    ax.plot(lm[0:35] + corr_cos, SFR[0:35] + corr_cos, color='PaleVioletRed', linewidth = 3, linestyle='dashed', label='Brinchmann+04')
    ax.plot(lm[36:70] + corr_cos, SFR[36:70] + corr_cos, color='PaleVioletRed',linewidth = 5, linestyle='dotted')
    ax.plot(lm[71:len(SFR)] + corr_cos, SFR[71:len(SFR)] + corr_cos, color='PaleVioletRed',linewidth = 5, linestyle='dotted')

    xdataD16 = [9.3, 10.6]
    ydataD16 = [-0.39, 0.477]
    ax.plot(xdataD16,ydataD16, color='b',linestyle='dashdot',linewidth = 4, label='Davies+16')

    ax.plot(xmf, (slope_ms_z0 * xmf + offset_ms_z0) + xmf - 9.0, linewidth = 3, linestyle='solid', color='grey',label='MS fit Shark')
    # Legend
    common.prepare_legend(ax, ['k','PaleVioletRed', 'b', 'grey'], loc=2)
    common.savefig(outdir, fig, 'SFR_Mstars_z0_MSShark.pdf')


    #plot scatter of the MS for all different definitions.
    fig = plt.figure(figsize=(5,5))
    ytit="$\\rm \\sigma_{\\rm SSFR}$"
    xtit="$\\rm log_{10}(\\rm M_{\\star}/M_{\odot})$"

    xmin, xmax, ymin, ymax = 7, 11, 0 , 1.4
    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(7.2,1.3,'$\\rm SSFR>\\alpha MS(M_{\\star},fit)$')

    #predicted relation
    inputs = [0,1,2,3,4,5]
    labels = ['$\\rm log_{10}(SSFR/Gyr^{-1})>-2.5$','$\\alpha=0.05$','$\\alpha=0.1$','$\\alpha=0.12$','$\\alpha=0.158$','$\\alpha=0.25$']
    colors = ['k','b','LightSteelBlue','Orange','Red','DarkSalmon']
    lines  = ['solid','dashed','dotted','dashdot','solid','dotted']
    for j in zip(inputs[:]):
        i = j[0]
        ind = np.where(sigmamainseqsf[0,i,:] != 0)
        yplot = sigmamainseqsf[0,i,ind]
        ax.plot(xmf[ind],yplot[0],color=colors[i],linestyle=lines[i], linewidth = 1, label=labels[i])

    common.prepare_legend(ax, colors,  bbox_to_anchor=(0.1, 0.45))
    common.savefig(outdir, fig, 'Scatter_mainsquence_z0.pdf')


    #plot scatter of the MS across redshifts.
    fig = plt.figure(figsize=(5,5))
    ytit="$\\rm \\sigma_{\\rm SSFR}$"
    xtit="$\\rm log_{10}(\\rm M_{\\star}/M_{\odot})$"

    xmin, xmax, ymin, ymax = 7, 11, 0 , 1.4
    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #predicted relation
    inputs = [0,1,2,3,4,5]
    labels = ['z=0','z=0.5','z=1','z=2','z=3','z=4']
    colors = ['k','b','DarkTurquoise','YellowGreen','Gold','red']
    lines  = ['solid','dashed','dotted','dashdot','solid','dotted']
    for j in zip(inputs[:]):
        i = j[0]
        ind = np.where(sigmamainseqsf[i,4,:] != 0)
        yplot = sigmamainseqsf[i,4,ind]
        ax.plot(xmf[ind],yplot[0],color=colors[i],linestyle=lines[i], linewidth = 1, label=labels[i])

    common.prepare_legend(ax, colors,  bbox_to_anchor=(0.1, 0.45))
    common.savefig(outdir, fig, 'Scatter_mainsquence_z_evolution.pdf')

def plot_passive_fraction(plt, outdir, obsdir, passive_fractions, hist_ssfr, passive_fractions_cens_sats):

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


    
    #plot passive fraction for satellites/centrals in bins of stellar mass and halo mass
    fig = plt.figure(figsize=(12,9))
    xtit="$\\rm log_{10} (\\rm M_{\\rm halo}/M_{\odot})$"
    ytit="$\\rm passive\,fraction$"

    subplots= (231, 232, 233, 234, 235, 236)
    titles=('$\\rm 9<log_{10} (\\rm M_{\\star}/M_{\odot})<9.4$','$\\rm 9.4<log_{10} (\\rm M_{\\star}/M_{\odot})<9.8$',
            '$\\rm 9.8<log_{10} (\\rm M_{\\star}/M_{\odot})<10.2$', '$\\rm 10.2<log_{10} (\\rm M_{\\star}/M_{\odot})<10.6$',
            '$\\rm 10.6<log_{10} (\\rm M_{\\star}/M_{\odot})<11$','$\\rm 11<log_{10} (\\rm M_{\\star}/M_{\odot})<11.4$')
    xmin, xmax, ymin, ymax = 11, 14.5, 0, 1

    for j in range(0,6):
        ax = fig.add_subplot(subplots[j])
        plt.subplots_adjust(bottom=0.15, left=0.15)

        xtitin = xtit
        if(j < 3):
           xtitin = ""
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitin, ytit, locators=(0.1, 1, 0.1, 1))
        ax.text(11.1, 1.03, titles[j])
        
        ind   = np.where(passive_fractions_cens_sats[0,0,j,:] != 0)
        xplot = xmf2[ind] 
        yplot = passive_fractions_cens_sats[0,0,j,ind]
        ax.plot(xplot,yplot[0],color='b',linestyle='solid', linewidth = 1, label="centrals")

        ind   = np.where(passive_fractions_cens_sats[0,1,j,:] != 0)
        xplot = xmf2[ind] 
        yplot = passive_fractions_cens_sats[0,1,j,ind]
        ax.plot(xplot,yplot[0],color='r',linestyle='dashed', linewidth = 1, label="satellites")

        if(j == 0):
           # Legend
           common.prepare_legend(ax, ['b','r'], loc=2)

    common.savefig(outdir, fig, 'passive_fraction_z0_mhalo-mstars.pdf')


def prepare_data(hdf5_data, index, redshift, mainseqsf, passive_fractions, hist_ssfr, sigmamainseqsf, 
                 passive_fractions_cens_sats):

    (h0, volh, sfr_disk, sfr_burst, mdisk, mbulge, rstar_disk, mBH, mHI, mH2, 
     mgas_disk, mHI_bulge, mH2_bulge, mgas_bulge, mgas_metals_disk, mgas_metals_bulge, 
     mstars_metals_disk, mstars_metals_bulge, typeg, mvir_hosthalo, rstar_bulge) = hdf5_data

    bin_it = functools.partial(us.wmedians, xbins=xmf)
    bin_it_2sigma = functools.partial(us.wmedians_2sigma, xbins=xmf)

    mgas = mgas_disk+mgas_bulge
    mgas_metals = mgas_metals_disk+mgas_metals_bulge

    mass          = np.zeros(shape = len(mdisk))
    mhalo         = np.zeros(shape = len(mdisk))
    sfr           = np.zeros(shape = len(mdisk))
    active_flag   = np.zeros(shape = len(mdisk))

    ind = np.where(mdisk+mbulge > 0)
    mass[ind] = np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0))
    mhalo[ind] = np.log10(mvir_hosthalo[ind]) - np.log10(float(h0))

    ind = np.where(sfr_disk+sfr_burst > 0)
    sfr[ind] = np.log10((sfr_disk[ind]+sfr_burst[ind])/h0/GyrToYr)

    ind = np.where((sfr_disk+sfr_burst > 0) & (mdisk+mbulge > 0))
    passive_fractions[index,0,:] = us.fractions(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)), y = np.log10((sfr_disk[ind]+sfr_burst[ind])/(mdisk[ind]+mbulge[ind])), xbins=xmf2, ythresh=-2.2)
    passive_fractions[index,0,:] = 1.0 - passive_fractions[index,0,:] 
    H, _ = np.histogram(np.log10((sfr_disk[ind]+sfr_burst[ind])/(mdisk[ind]+mbulge[ind])),bins=np.append(ssfrbins,ssfrupp))
    hist_ssfr[index,:] = hist_ssfr[index,:] + H

    ind = np.where((sfr_disk+sfr_burst > 0) & (mdisk+mbulge > 0) & (mvir_hosthalo < 1e11))
    passive_fractions[index,1,:] = us.fractions(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)), y = np.log10((sfr_disk[ind]+sfr_burst[ind])/(mdisk[ind]+mbulge[ind])), xbins=xmf2, ythresh=-2.2)
    passive_fractions[index,1,:] = 1.0 - passive_fractions[index,1,:] 

    ind = np.where((sfr_disk+sfr_burst > 0) & (mdisk+mbulge > 0) & (mvir_hosthalo >= 1e11))
    passive_fractions[index,2,:] = us.fractions(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)), y = np.log10((sfr_disk[ind]+sfr_burst[ind])/(mdisk[ind]+mbulge[ind])), xbins=xmf2, ythresh=-2.2)
    passive_fractions[index,2,:] = 1.0 - passive_fractions[index,2,:] 

    ind = np.where((sfr_disk+sfr_burst > 0) & (mdisk+mbulge > 0) & ((sfr_disk+sfr_burst)/(mdisk+mbulge) > 1e-3))
    mainseqsf[index,:] = bin_it_2sigma(x=mass[ind], y=np.log10((sfr_disk[ind]+sfr_burst[ind])/h0/GyrToYr))

    # calculate main sequence:
    ms = np.zeros(shape = len(xmf))
    for j in range(0,len(xmf)):
        ind = np.where((mass > xmf[j]-dm/2.0) & (mass <= xmf[j]+dm/2.0) & (sfr - mass + 9.0 > -5 + 0.5*redshift) & (typeg == 0))
        if(len(mass[ind] > 0)):
               ms[j] = np.median(sfr[ind] - mass[ind] + 9.0)
    ind = np.where((ms != 0) & (xmf > 8) & (xmf < 9.8))
    (ms_fit_slope, ms_fit_offs) = np.polyfit(xmf[ind],ms[ind],1)

    # calculate main sequence. sfr, m need to be given in log10 units.
    def calculate_sigma_sfr_fromfixssfr(sfr, m, sigma, mscut, offms, flag):
        ms = np.zeros(shape = len(xmf))
        for j in range(0,len(xmf)):
 	    ind = np.where((m > xmf[j]-dm/2.0) & (m <= xmf[j]+dm/2.0) & (sfr - m + 9.0 > mscut + np.log10(offms)))
            if(len(m[ind] > 0)):
                   sigma[j] = np.std(sfr[ind] - m[ind] + 9.0)
                   flag[ind] = 1

    def calculate_sigma_sfr_fromms(sfr, m, sigma, mscut, offms, flag):
        ms = np.zeros(shape = len(xmf))
        for j in range(0,len(xmf)):
            ind = np.where((m > xmf[j]-dm/2.0) & (m <= xmf[j]+dm/2.0) & (sfr - m + 9.0 > mscut))
            if(len(m[ind] > 0)):
                   ms[j] = np.median(sfr[ind] - m[ind] + 9.0)
 	    ind = np.where((m > xmf[j]-dm/2.0) & (m <= xmf[j]+dm/2.0) & (sfr - m + 9.0 > ms[j]+np.log10(offms)))
            if(len(m[ind] > 0)):
                   sigma[j] = np.std(sfr[ind] - m[ind] + 9.0)
                   flag[ind] = 1

    def calculate_sigma_sfr_frommsfit(sfr, m, sigma, ms_fit_slope, ms_fit_offs, offms, flag):
        ms = np.zeros(shape = len(xmf))
        for j in range(0,len(xmf)):
 	    ind = np.where((m > xmf[j]-dm/2.0) & (m <= xmf[j]+dm/2.0) & (sfr - m + 9.0 > ms_fit_slope * m + ms_fit_offs + np.log10(offms)))
            if(len(m[ind] > 0)):
                   sigma[j]  = np.std(sfr[ind] - m[ind] + 9.0)
                   flag[ind] = 1

    ind = np.where((mass > 0) & (sfr != 0))
    calculate_sigma_sfr_fromfixssfr(sfr[ind], mass[ind], sigmamainseqsf[index,0,:], -2.0, 1.0, active_flag[ind])
    calculate_sigma_sfr_frommsfit(sfr[ind], mass[ind], sigmamainseqsf[index,1,:], ms_fit_slope,ms_fit_offs, 1.0/20.0, active_flag[ind])
    calculate_sigma_sfr_frommsfit(sfr[ind], mass[ind], sigmamainseqsf[index,2,:], ms_fit_slope,ms_fit_offs, 1.0/10.0, active_flag[ind])
    calculate_sigma_sfr_frommsfit(sfr[ind], mass[ind], sigmamainseqsf[index,3,:], ms_fit_slope,ms_fit_offs, 1.0/8.0, active_flag[ind])
    calculate_sigma_sfr_frommsfit(sfr[ind], mass[ind], sigmamainseqsf[index,4,:], ms_fit_slope,ms_fit_offs, 1.0/6.0, active_flag[ind])
    calculate_sigma_sfr_frommsfit(sfr[ind], mass[ind], sigmamainseqsf[index,5,:], ms_fit_slope,ms_fit_offs, 1.0/4.0, active_flag[ind])

    calculate_sigma_sfr_frommsfit(sfr, mass, sigmamainseqsf[index,6,:], ms_fit_slope, ms_fit_offs, 0.15848931924, active_flag)


    for i in range(0,len(xmflr)):
        for j in range(0,len(xmf2)):
            #centrals
            ind = np.where((mass > xmflr[i]-dmlr/2.0) & (mass <= xmflr[i]+dmlr/2.0) & (mhalo > xmf2[j] - dm2/2.0) & (mhalo <= xmf2[j] + dm2/2.0) & (typeg <= 0))
            if(len(mass[ind]) > 9):
                   totnumber      = len(mass[ind])
                   onms           = active_flag[ind]
                   moffs          = np.where(onms <= 0)
                   passivenumber  = len(onms[moffs])
                   passive_fractions_cens_sats[index,0,i,j] = (passivenumber + 0.0)/(totnumber + 0.0)
            ind = np.where((mass > xmflr[i]-dmlr/2.0) & (mass <= xmflr[i]+dmlr/2.0) & (mhalo > xmf2[j] - dm2/2.0) & (mhalo <= xmf2[j] + dm2/2.0) & (typeg > 0))
            if(len(mass[ind]) > 9):
                   totnumber      = len(mass[ind])
                   onms           = active_flag[ind]
                   moffs          = np.where(onms <= 0)
                   passivenumber  = len(onms[moffs])
                   passive_fractions_cens_sats[index,1,i,j] = (passivenumber + 0.0)/(totnumber + 0.0)

    return (mass, ms_fit_slope, ms_fit_offs)

def main(modeldir, outdir, subvols, obsdir):

    zlist = ["199","174", "156", "131", "113", "99"]
    redshifts = [0, 0.5, 1.0, 2.0, 3.0, 4.0]
    plt = common.load_matplotlib()

    mainseqsf          = np.zeros(shape = (len(zlist), 3, len(xmf)))
    sigmamainseqsf     = np.zeros(shape = (len(zlist), 7, len(xmf)))

    passive_fractions = np.zeros(shape = (len(zlist), 3, len(xmf2)))
    passive_fractions_cens_sats = np.zeros(shape = (len(zlist), 2, len(xmflr), len(xmf2)))

    hist_ssfr = np.zeros(shape = (len(zlist), len(ssfrbins)))

    fields = {'galaxies': ('sfr_disk', 'sfr_burst', 'mstars_disk', 'mstars_bulge',
                           'rstar_disk', 'm_bh', 'matom_disk', 'mmol_disk', 'mgas_disk',
                           'matom_bulge', 'mmol_bulge', 'mgas_bulge',
                           'mgas_metals_disk', 'mgas_metals_bulge',
                           'mstars_metals_disk', 'mstars_metals_bulge', 'type', 
			   'mvir_hosthalo', 'rstar_bulge')}

    for index in range(0,len(zlist)):
        hdf5_data = common.read_data(modeldir, zlist[index], fields, subvols)
        (mass, slope, offset) = prepare_data(hdf5_data, index, redshifts[index], mainseqsf, passive_fractions, 
                            hist_ssfr, sigmamainseqsf, passive_fractions_cens_sats)

        h0 = hdf5_data[0]
        if index == 0:
            (sfr_disk, sfr_burst, mdisk, mbulge) = hdf5_data[2:6]
            sfr_seq = np.zeros(shape = (2, len(mdisk)))
            ind  = np.where((sfr_disk + sfr_burst > 0) & (mdisk + mbulge > 0))
            sfr_seq[0,ind] = mass[ind]
            sfr_seq[1,ind] = np.log10((sfr_disk[ind] + sfr_burst[ind]) / h0 / GyrToYr)
            slope_ms_z0  = slope
            offset_ms_z0 = offset
            print 'scatter MS'
            for m,a,b,c,d,e,f,g in zip(xmf[:], sigmamainseqsf[index,0,:], sigmamainseqsf[index,1,:], sigmamainseqsf[index,2,:], sigmamainseqsf[index,3,:], sigmamainseqsf[index,4,:], sigmamainseqsf[index,5,:], sigmamainseqsf[index,6,:]):
                print m,a,b,c,d,e,f,g

            print 'passive fractions centrals'
            for m,a,b,c,d,e,f in zip(xmf2[:], passive_fractions_cens_sats[0,0,0,:], passive_fractions_cens_sats[0,0,1,:], passive_fractions_cens_sats[0,0,2,:], passive_fractions_cens_sats[0,0,3,:], passive_fractions_cens_sats[0,0,4,:], passive_fractions_cens_sats[0,0,5,:],):
                print m,a,b,c,d,e,f
            print 'passive fractions satellites'
            for m,a,b,c,d,e,f in zip(xmf2[:], passive_fractions_cens_sats[0,1,0,:], passive_fractions_cens_sats[0,1,1,:], passive_fractions_cens_sats[0,1,2,:], passive_fractions_cens_sats[0,1,3,:], passive_fractions_cens_sats[0,1,4,:], passive_fractions_cens_sats[0,1,5,:],):
                print m,a,b,c,d,e,f


    # This should be the same in all HDF5 files
    plot_sfr_mstars_z0(plt, outdir, obsdir, h0, sfr_seq, mainseqsf, sigmamainseqsf, slope_ms_z0, offset_ms_z0)
    plot_passive_fraction(plt, outdir, obsdir, passive_fractions, hist_ssfr, 
                          passive_fractions_cens_sats) 


if __name__ == '__main__':
    main(*common.parse_args(requires_snapshot=False))
