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
mlow = 7
mupp = 15
dm = 0.125
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

mlow2 = 8
mupp2 = 12
dm2 = 0.15
mbins2 = np.arange(mlow2,mupp2,dm2)
xmf2 = mbins2 + dm2/2.0

def load_smf_passive_observations(obsdir, h0):

    # Weaver et al. (2022; COSMOS2020)
    def add_weaver22_data(zappend, file_read='0.2z0.5'):
        lm, pD, dn, du = np.loadtxt(obsdir+'/mf/SMF/COSMOS2020/SMF_Farmer_v2.1_' + file_read + '_qg.txt', delimiter=' ', skiprows=0, usecols = [0,2,3,4], unpack = True)
        hobs = 0.7
        ind = np.where(dn < 0)
        dn[ind] = 1e-9
        pDlog = np.log10(pD) +  3.0 * np.log10(hobs/h0)
        dnlog = np.log10(pD) - np.log10(dn)
        dulog = np.log10(du) - np.log10(pD)
        lm = lm -  2.0 * np.log10(hobs/h0)
        zappend.append((observation("Weaver+2023", lm, pDlog, abs(dulog), abs(dnlog), err_absolute=False), 'D'))

    # z0.5 obs
    z05obs = []
    add_weaver22_data(z05obs, file_read='0.2z0.5')
    add_weaver22_data(z05obs, file_read='0.5z0.8')

    # z1 obs
    z1obs = []
    add_weaver22_data(z1obs, file_read='0.8z1.1')
    add_weaver22_data(z1obs, file_read='1.1z1.5')


    #z2 obs
    z2obs = []
    add_weaver22_data(z2obs, file_read='2.0z2.5')
    add_weaver22_data(z2obs, file_read='1.5z2.0')

    # z3 obs
    z3obs = []
    add_weaver22_data(z3obs, file_read='3.0z3.5')
    add_weaver22_data(z3obs, file_read='2.5z3.0')


    # z4 obs
    z4obs = []
    add_weaver22_data(z4obs, file_read='3.5z4.5')

    # z5 obs
    z5obs = []
    add_weaver22_data(z5obs,file_read='4.5z5.5')


    return (z05obs, z1obs, z2obs, z3obs, z4obs, z5obs)

def plot_SMHM_z(plt, outdir, obsdir, zlist, halo_mass_rel):

    def plot_ssfr_based_selection(ax, z, label=True):
        mh, med, percenlow, percenup, zin, flagp = common.load_observation(obsdir, 'Models/SharkVariations/SMHM_Galaxies_L23.dat', [0,1,2,3,4,5])
        ind = np.where((zin == z) & (flagp == 1))
        ax.fill_between(mh[ind], percenlow[ind], percenup[ind], facecolor='salmon', alpha = 0.5, interpolate=True)
        ax.plot(mh[ind], med[ind], linestyle='dashed', linewidth=3, color='salmon',label='P (sSFR)' if label == True else None)

        ind = np.where((zin == z) & (flagp == 0))
        ax.fill_between(mh[ind], percenlow[ind], percenup[ind], facecolor='DeepSkyBlue', alpha = 0.5, interpolate=True)
        ax.plot(mh[ind], med[ind], linestyle='dashed', linewidth=3, color='DeepSkyBlue',label='SF (sSFR)' if label == True else None)



    def plot_moster13(ax, z, labels, label):
        #Moster et al. (2013) abundance matching SMHM relation
        M10 = 11.590
        M11 = 1.195
        N10 = 0.0351
        N11 = -0.0247
        beta10 = 1.376
        beta11 = -0.826
        gamma10 = 0.608
        gamma11 = 0.329
        M1 = pow(10.0, M10 + M11 * z/(z+1))
        N = N10 + N11 * z/(z+1)
        beta = beta10 + beta11 * z/(z+1)
        gamma = gamma10 + gamma11 * z/(z+1)

        mh = pow(10.0,xmf)
        m = mh * 2*N * pow (( pow(mh/M1, -beta ) + pow(mh/M1, gamma)), -1)

        if not labels:
            ax.plot(xmf,np.log10(m),'magenta', linestyle='dashed', linewidth=3)
        else:
            ax.plot(xmf,np.log10(m),'magenta', linestyle='dashed', linewidth=3, label=label)

    def plot_berhoozi13(ax, z, labels, label):
        a = 1.0/(1.0+z)
        nu = np.exp(-4*a*a)
        log_epsilon = -1.777 + (-0.006*(a-1)) * nu
        M1= 11.514 + ( - 1.793 * (a-1) - 0.251 * z) * nu
        alpha = -1.412 + 0.731 * nu * (a-1)
        delta = 3.508 + (2.608*(a-1)-0.043 * z) * nu
        gamma = 0.316 + (1.319*(a-1)+0.279 * z) * nu
        Min = xmf-M1
        fx = -np.log10(pow(10,alpha*Min)+1.0)+ delta * pow(np.log10(1+np.exp(Min)),gamma) / (1+np.exp(pow(10,-Min)))
        f = -0.3+ delta * pow(np.log10(2.0),gamma) / (1+np.exp(1))

        m = log_epsilon + M1 + fx - f

        if not labels:
            ax.plot(xmf,m, 'b', linestyle='dashdot',linewidth=3)
        else:
            ax.plot(xmf,m, 'b', linestyle='dashdot',linewidth=3, label=label)



    fig = plt.figure(figsize=(5,9.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm halo}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    xmin, xmax, ymin, ymax = 11, 14, 9, 12
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymin + 0.1 * (ymax - ymin)


    subplots = (311, 312, 313)
    #zlist = (0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0)
    indices = [2,3,4]
    for i,idx in enumerate(indices):
        
        # z=0 ##################################
        ax = fig.add_subplot(subplots[i])
        if(i == 2):
            xtitle = xtit
        else:
            xtitle = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitle, ytit, locators=(0.1, 1, 0.1))
        ax.text(xleg, yleg, 'z=%s' % str(zlist[idx]))

        ax.tick_params(labelsize=13)

        if i == 0:
            plot_ssfr_based_selection(ax, zlist[idx], label=True)
        else:
            plot_ssfr_based_selection(ax, zlist[idx], label=False)

        #Predicted SMHM
        ind = np.where(halo_mass_rel[idx,0,0,:] != 0)
        xplot = xmf[ind]
        yplot = halo_mass_rel[idx,0,0,ind]
        errdn = halo_mass_rel[idx,0,1,ind]
        errup = halo_mass_rel[idx,0,2,ind]
   
        ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='red', alpha = 0.5, interpolate=True)
        ax.errorbar(xplot, yplot[0], color='r', label='P (colour)')

        ind = np.where(halo_mass_rel[idx,2,0,:] != 0)
        xplot = xmf[ind]
        yplot = halo_mass_rel[idx,2,0,ind]
        errdn = halo_mass_rel[idx,2,1,ind]
        errup = halo_mass_rel[idx,2,2,ind]
   
        ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='RoyalBlue', alpha = 0.5, interpolate=True)
        ax.errorbar(xplot, yplot[0], color='navy', label = 'SF (colour)')
 

        #plot_moster13(ax, zlist[idx], True, 'Moster+13')
        #plot_berhoozi13(ax, zlist[idx], True, 'Behroozi+13')

        if(i == 0):
           common.prepare_legend(ax, ['Salmon','DeepSkyBlue','r','navy','magenta','b'], loc=2)

    plt.tight_layout()
    common.savefig(outdir, fig, 'SMHM_z_passivegals_colorbased.pdf')


def plot_stellarmf_passive_z(plt, outdir, obsdir, h0, hist_smf, hist_smf_err, hist_smf_pass_cen, hist_smf_pass_sat):

    (z05obs, z1obs, z2obs, z3obs, z4obs, z5obs) = load_smf_passive_observations(obsdir, h0)
    PlotLagos23 = True
    def plot_lagos23_smf(ax, z, label=True):
        sm, z0, z0p5, z1, z2, z3, z4, z5 = common.load_observation(obsdir, 'Models/SharkVariations/SMF_Passive_Centrals_Lagos23.dat', [0,1,2,3,4,5,6,7])
        y = z0
        if z == 1:
           y = z0p5
        elif z == 2:
           y = z1
        elif z == 3:
           y = z2
        elif z == 4:
           y = z3
        elif z == 5:
           y = z4
        elif z == 6:
           y = z5
        ax.plot(sm, y, linestyle='dashed', color='DarkOrange',label='centrals (sSFR)' if label == True else None)

        sm, z0, z0p5, z1, z2, z3, z4, z5 = common.load_observation(obsdir, 'Models/SharkVariations/SMF_Passive_Satellites_Lagos23.dat', [0,1,2,3,4,5,6,7])
        y = z0
        if z == 1:
           y = z0p5
        elif z == 2:
           y = z1
        elif z == 3:
           y = z2
        elif z == 4:
           y = z3
        elif z == 5:
           y = z4
        elif z == 6:
           y = z5
        ax.plot(sm, y, linestyle='dashed', color='YellowGreen',label='satellites (sSFR)' if label == True else None)
  
    fig = plt.figure(figsize=(11.7,8))
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\Phi/dlog_{10}{\\rm M_{\\star}}/{\\rm Mpc}^{-3} )$"
    xmin, xmax, ymin, ymax = 8, 13, -6.5, -1
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (231, 232, 233, 234, 235, 236)
    indeces = (1, 2, 3, 4, 5, 6)
    zs = (0.5, 1, 2, 3, 4, 5)
    observations = (z05obs, z1obs, z2obs, z3obs, z4obs, z5obs)
    limits_obs = [7, 8.4, 9.4, 10, 10.1, 10.2, 10.45]

    for subplot, idx, z, obs_and_markers in zip(subplots, indeces, zs, observations):
   
        ax = fig.add_subplot(subplot)
        if(idx ==1 or idx == 4):
            ytitle = ytit
        else:
            ytitle = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytitle, locators=(0.1, 1, 0.1))
        ax.text(xleg, yleg, 'z=%s' % str(z))
        xlim = limits_obs[idx]
        ax.plot([xlim, xlim], [-6.5,-1], linestyle='dotted', color='grey')
        # Observations
        for obs, marker in obs_and_markers:
            if idx == 6:
                common.errorbars(ax, obs.x, obs.y, obs.yerrdn, obs.yerrup, 'grey',
                                marker, err_absolute=obs.err_absolute, label=obs.label, markersize=6)
            else:
                common.errorbars(ax, obs.x, obs.y, obs.yerrdn, obs.yerrup, 'grey',
                                marker, err_absolute=obs.err_absolute, markersize=6)


        # Predicted SMF
        #y = hist_smf[idx,:]
        #ind = np.where(y < 0.)
        #ax.plot(xmf[ind],y[ind],'r', label='Shark v2.0 (color)' if idx == 1 else None)
        #y = hist_smf_err[idx,:]
        #ind = np.where(y < 0.)
        #ax.plot(xmf[ind],y[ind],'r', linestyle='solid', linewidth=2, label ='Shark v2.0 (color)' if idx == 1 else None)
        y = hist_smf_pass_cen[idx,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'DarkOrange', linestyle='solid', label='centrals (color)' if idx == 1 else None)
        y = hist_smf_pass_sat[idx,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'YellowGreen', label='satellites (color)' if idx == 1 else None)

        if(PlotLagos23):
            if(idx == 1):
               plot_lagos23_smf(ax, idx, label=True)
            else:
               plot_lagos23_smf(ax, idx, label=False)

        colors = []
        if idx == 1:
            colors = ['DarkOrange','YellowGreen','DarkOrange','YellowGreen']
        if idx == 6:
            colors=['grey','grey']
        colors += ['grey', 'grey','grey']

        if idx ==1:
           common.prepare_legend(ax, colors, loc = 3)
        elif idx == 6:
           common.prepare_legend(ax, colors, loc = 2)


    plt.tight_layout()
    common.savefig(outdir, fig, 'stellarmf_passive_z_colorbased.pdf')


def prepare_data(hdf5_data, seds, index, hist_smf_pass, hist_smf_pass_err, hist_smf_pass_cen, hist_smf_pass_sat, zlist, halo_mass_rel):

    (h0, volh, mdisk, mbulge, sfrd, sfrb, typeg, rstar_disk, rstar_bulge, mvir) = hdf5_data

    bin_it = functools.partial(us.wmedians, xbins=xmf, nmin=10)

   #(0): "FUV_GALEX", "NUV_GALEX", "u_SDSS", "g_SDSS", "r_SDSS", "i_SDSS",
   #(6): "z_SDSS", "Y_VISTA", "J_VISTA", "H_VISTA", "K_VISTA", "W1_WISE",
   #(12): "I1_Spitzer", "I2_Spitzer", "W2_WISE", "I3_Spitzer", "I4_Spitzer",
   #(17): "W3_WISE", "W4_WISE", "P70_Herschel", "P100_Herschel",
   #(21): "P160_Herschel", "S250_Herschel", "S350_Herschel", "S450_JCMT",
   #(25): "S500_Herschel", "S850_JCMT", "Band_ionising_photons", "FUV_Nathan",
   #(29): "Band9_ALMA", "Band8_ALMA", "Band7_ALMA", "Band6_ALMA",
   #(33): "Band4_ALMA", "Band3_ALMA", "BandX_VLA", "BandC_VLA", "BandS_VLA",
   #(38): "BandL_VLA", "Band_610MHz", "Band_325MHz", "Band_150MHz"

    mag_total = seds[0]
    col1 = mag_total[1,:] - mag_total[4,:]
    col2 = mag_total[4,:] -  mag_total[8,:]


    ind = np.where((mdisk+mbulge) > 0.0)
    mvir = mvir[ind]
    typeg = typeg[ind]
    pass_flag = np.zeros(shape = len(mdisk[ind]))
    passive_col = np.where((col1 > 3*col2+1) & (col1 > 3.1))
    pass_flag[passive_col] = 1.0
    mass = np.log10(mdisk[ind] + mbulge[ind]) - np.log10(float(h0))
    ssfr = (sfrd[ind] + sfrb[ind]) / 1e9 / (mdisk[ind] + mbulge[ind])
    lowsfr = np.where(ssfr <= 1e-13)
    ssfr[lowsfr] = 1e-13
    ssfr = np.log10(ssfr)

    ran_err = np.random.normal(0.0, 0.3, len(mass))
    mass_err = mass + ran_err

    ind = np.where((mass > 0) & (pass_flag >= 1))
    H, _ = np.histogram(mass[ind], bins=np.append(mbins,mupp))
    hist_smf_pass[index,:] = hist_smf_pass[index,:] + H

    if index == 0:
        scatter = 0.2
    else:
        scatter = 0.3
    ran_err = np.random.normal(0.0, scatter, len(mass))
    ssfr_err = ssfr + ran_err
    ind = np.where((mass_err > 0) & (pass_flag >=1))
    H, _ = np.histogram(mass_err[ind], bins=np.append(mbins,mupp))
    hist_smf_pass_err[index,:] = hist_smf_pass_err[index,:] + H

    ind = np.where((mass_err > 0) & (pass_flag >=1) & (typeg == 0))
    H, _ = np.histogram(mass_err[ind], bins=np.append(mbins,mupp))
    hist_smf_pass_cen[index,:] = hist_smf_pass_cen[index,:] + H
    halo_mass_rel[index,1,:] = bin_it(x=np.log10(mvir[ind]), y = mass_err[ind])

    ind = np.where((mass > 0) & (pass_flag >=1) & (typeg == 0))
    halo_mass_rel[index,0,:] = bin_it(x=np.log10(mvir[ind]), y = mass[ind])

    ind = np.where((mass_err > 0) & (pass_flag >=1) & (typeg > 0))
    H, _ = np.histogram(mass_err[ind], bins=np.append(mbins,mupp))
    hist_smf_pass_sat[index,:] = hist_smf_pass_sat[index,:] + H

    if volh > 0:
        vol = volh/pow(h0,3.)  # In Mpc^3
        hist_smf_pass[index,:]  = hist_smf_pass[index,:]/vol/dm
        hist_smf_pass_err[index,:]  = hist_smf_pass_err[index,:]/vol/dm
        hist_smf_pass_cen[index,:]  = hist_smf_pass_cen[index,:]/vol/dm
        hist_smf_pass_sat[index,:]  = hist_smf_pass_sat[index,:]/vol/dm


    ind = np.where((mass > 0) & (typeg == 0))
    halo_mass_rel[index,3,:] = bin_it(x=np.log10(mvir[ind]), y = mass[ind])

    ind = np.where((mass > 0) & (pass_flag <= 0) & (typeg == 0))
    halo_mass_rel[index,2,:] = bin_it(x=np.log10(mvir[ind]), y = mass[ind])

    return mass

def main(modeldir, outdir, redshift_table, subvols, obsdir):

    zlist = (0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0) #(0.15, 0.25, 0.4, 0.625, 0.875, 1.125, 1.375, 1.625, 1.875)

    plt = common.load_matplotlib()
    file_hdf5_sed = "Shark-SED-eagle-rr14.hdf5"

    # Histograms
    hist_smf_pass   = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_pass_cen   = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_pass_sat   = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_pass_err = np.zeros(shape = (len(zlist), len(mbins)))
    halo_mass_rel  = np.zeros(shape = (len(zlist), 4, 3, len(mbins)))

    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'sfr_disk', 'sfr_burst', 'type', 'rstar_disk', 'rstar_bulge', 'mvir_hosthalo')}
    fields_sed = {'SED/ab_dust': ('total','disk')}
    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        seds = common.read_photometry_data_variable_tau_screen(modeldir, snapshot, fields_sed, subvols, file_hdf5_sed)
        mass = prepare_data(hdf5_data, seds, index, hist_smf_pass, hist_smf_pass_err, hist_smf_pass_cen, hist_smf_pass_sat, zlist, halo_mass_rel)
        h0 = hdf5_data[0]

    # Take logs
    def take_log(array):
        ind = np.where(array > 0.)
        array[ind] = np.log10(array[ind])
    take_log(hist_smf_pass)
    take_log(hist_smf_pass_err)
    take_log(hist_smf_pass_cen)
    take_log(hist_smf_pass_sat)

    plot_stellarmf_passive_z(plt, outdir, obsdir, h0, hist_smf_pass, hist_smf_pass_err, hist_smf_pass_cen, hist_smf_pass_sat)
    plot_SMHM_z(plt, outdir, obsdir, zlist, halo_mass_rel)
    #print("#SMF passive galaxies")
    #for a,b,c,d,e,f,g,h in zip(xmf, hist_smf_pass_cen[0,:], hist_smf_pass_cen[1,:], hist_smf_pass_cen[2,:], hist_smf_pass_cen[3,:], hist_smf_pass_cen[4,:], hist_smf_pass_cen[5,:], hist_smf_pass_cen[6,:]):
    #    print(a,b,c,d,e,f,g,h)

if __name__ == '__main__':
    main(*common.parse_args())
