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
import h5py
#import hickle

import common
import utilities_statistics as us

##################################
# Constants
mlow = 7.0
mupp = 13.0
dm = 0.5
mbins = np.arange(mlow, mupp, dm)
xmf = mbins + dm/2.0

GyrtoYr  = 1e9
MpcToKpc = 1e3
G        = 4.299e-9 #Gravity constant in units of (km/s)^2 * Mpc/Msun
PI       = 3.1416

def add_observations_to_plot(obsdir, fname, ax, marker, label, color='k', err_absolute=False):
    fname = '%s/Gas/%s' % (obsdir, fname)
    x, y, yerr_down, yerr_up = common.load_observation(obsdir, fname, (0, 1, 2, 3))
    common.errorbars(ax, x, y, yerr_down, yerr_up, color, marker, label=label, err_absolute=err_absolute)

def prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit):
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)
    xleg = xmax - 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)
    #ax.text(xleg, yleg, 'z=0')

def plot_num_density_passive(plt, outdir, obs_dir, zlist, num_den_passive, num_den_passive_we, num_den_passive_mass, num_den_passive_we_mass):

    def plot_observations_num_density(ax, labels=True):
    
        z, zerr, n, nerrh, nerrl = np.loadtxt(obs_dir + '/num_densities/' + 'Gould23data.txt', unpack = True, usecols = [0,1,2,3,4])
        obs = ['Straatman+14', 'Gould+23', 'Carnall+20', 'Schreiber+18', 'Shahidi+20', 'Merlin+19', 'Girelli+19', 'Weaver+22', 'Carnall+23', 'Valentino+23', 'Nanayakkara+22']
        nerrh = np.log10(n+nerrh) - np.log10(n)
        nerrl = np.log10(n) - np.log10(n-nerrl)
        n = np.log10(n)
        Nlow = [0, 1, 4, 8, 9, 11, 13, 17, 22, 24, 28]
        lenght = [1, 3, 4, 1, 2, 2, 4, 5, 2, 4, 1]
        symbols = ['o', 's', '*', 'd', 'D', '<', '>' , 'p', 'h', 'P', 'X']
        for j in range(0,len(obs)):
            if(j <= 7):
                col = 'grey'
            else:
                col = 'darkgreen'
            if(labels):
               ax.errorbar(z[Nlow[j]:Nlow[j]+lenght[j]], n[Nlow[j]:Nlow[j]+lenght[j]], xerr=zerr[Nlow[j]:Nlow[j]+lenght[j]], yerr=[nerrl[Nlow[j]:Nlow[j]+lenght[j]], nerrh[Nlow[j]:Nlow[j]+lenght[j]]], color=col, ecolor=col, marker=symbols[j], ls='None', mfc='None', label=obs[j])
            else:
                ax.errorbar(z[Nlow[j]:Nlow[j]+lenght[j]], n[Nlow[j]:Nlow[j]+lenght[j]], xerr=zerr[Nlow[j]:Nlow[j]+lenght[j]], yerr=[nerrl[Nlow[j]:Nlow[j]+lenght[j]], nerrh[Nlow[j]:Nlow[j]+lenght[j]]], color=col, ecolor=col, marker=symbols[j], ls='None', mfc='None')

    def load_lagos18_predictions(ax, low_mass= True, labels=True):
        z, n1, n2, n3, n4 = np.loadtxt(obs_dir + '/Models/SharkVariations/NumDensityPassiveGals.dat', unpack = True, usecols = [0,1,2,3,4])
        if(low_mass):
            nplot1 = np.log10(n1)
            nplot2 = np.log10(n2)
        else:
            nplot1 = np.log10(n3)
            nplot2 = np.log10(n4)

        if(labels):
           ax.plot(z,nplot1,ls='solid', color='MidnightBlue', lw=1, label = '$\\rm sSFR/yr^{-1} < 10^{-10}$ v1.1')
           ax.plot(z,nplot2,ls='solid', color='LightCoral', lw=1, label = '$\\rm sSFR/yr^{-1} < 10^{-11}$ v1.1')
        else:
           ax.plot(z,nplot1,ls='solid', color='MidnightBlue', lw=1)
           ax.plot(z,nplot2,ls='solid', color='LightCoral', lw=1)

    def load_sharkv2p0_variations(ax, labels=True):
        z, n1, n2, n3, n4 = np.loadtxt(obs_dir + '/Models/SharkVariations/NumDensityPassiveGals_Sharkv2p0_NoQSOFeedback.dat', unpack = True, usecols = [0,1,2,3,4])

        if(labels):
           ax.plot(z,np.log10(n1),ls='solid', color='LightCoral', lw=2, label = 'QSO feedback off')
           ax.plot(z,np.log10(n2),ls='dashed', color='LightCoral', lw=2, label = 'QSO feedback off with error')
        else:
           ax.plot(z,np.log10(n1),ls='solid', color='LightCoral', lw=2)
           ax.plot(z,np.log10(n2),ls='dashed', color='LightCoral', lw=2)




    fig = plt.figure(figsize=(7,10))
    ytit = "$\\rm log_{10} (\\rm \\phi/Mpc^{-3})$"
    xtit = "$\\rm redshift$"
    xmin, xmax, ymin, ymax = 0, 5, -7.5, -2.5
    xleg = xmax - 0.3 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(211)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg,yleg,'$\\rm M_{\\star}>10^{10}\,M_{\\odot}$', fontsize=13)

    cols = ['blue','red']
    lines = ['solid','dashed']
    labels = ['$\\rm sSFR/yr^{-1} < 10^{-10}$ v2.0', '$\\rm sSFR/yr^{-1} < 10^{-11}$ v2.0']
    
    for i in range(0,2):
         #Predicted relation
         ind = np.where(num_den_passive[:,i] != 0)
         yplot = np.log10(num_den_passive[ind,i])
         xplot = zlist[ind]
         ax.plot(xplot,yplot[0],ls=lines[0], color=cols[i], lw=3)
         ind = np.where(num_den_passive_we[:,i] != 0)
         yplot = np.log10(num_den_passive_we[ind,i])
         xplot = zlist[ind]
         ax.plot(xplot,yplot[0],ls=lines[1], color=cols[i])
    plot_observations_num_density(ax, labels=True)
    load_lagos18_predictions(ax, low_mass= True, labels=False)
    common.prepare_legend(ax, ['grey','grey', 'grey', 'grey', 'grey','grey', 'grey', 'grey', 'darkgreen','darkgreen', 'darkgreen'], loc=3)

    ax = fig.add_subplot(212)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg,yleg,'$\\rm M_{\\star}>10^{10.5}\,M_{\\odot}$', fontsize=13)

    for i in range(0,2):
         #Predicted relation
         ind = np.where(num_den_passive_mass[:,i] != 0)
         yplot = np.log10(num_den_passive_mass[ind,i])
         xplot = zlist[ind]
         ax.plot(xplot,yplot[0],ls=lines[0], color=cols[i], lw = 3, label=labels[i])
         ind = np.where(num_den_passive_we_mass[:,i] != 0)
         yplot = np.log10(num_den_passive_we_mass[ind,i])
         xplot = zlist[ind]
         ax.plot(xplot,yplot[0],ls=lines[1], color=cols[i], label=labels[i] + ' with error')
    for a,b,c,d,e in zip(zlist, num_den_passive_mass[:,0], num_den_passive_mass[:,1], num_den_passive_we_mass[:,0], num_den_passive_we_mass[:,1]):
        print(a,b,c,d,e)

    plot_observations_num_density(ax, labels=False)
    load_lagos18_predictions(ax, low_mass= False, labels=True)

    common.prepare_legend(ax, ['b','b', 'r', 'r', 'MidnightBlue', 'LightCoral'], loc=3)
    plt.tight_layout()
    common.savefig(outdir, fig, 'num_density_passive.pdf')


    fig = plt.figure(figsize=(7,5))
    ytit = "$\\rm log_{10} (\\rm \\phi/Mpc^{-3})$"
    xtit = "$\\rm redshift$"
    xmin, xmax, ymin, ymax = 0, 5, -7.5, -2.5
    xleg = xmax - 0.3 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg,yleg,'$\\rm M_{\\star}>10^{10.5}\,M_{\\odot}$', fontsize=13)

    for i in range(0,2):
         #Predicted relation
         ind = np.where(num_den_passive_mass[:,i] != 0)
         yplot = np.log10(num_den_passive_mass[ind,i])
         xplot = zlist[ind]
         ax.plot(xplot,yplot[0],ls=lines[0], color=cols[i], lw = 3, label=labels[i])
         #ind = np.where(num_den_passive_we_mass[:,i] != 0)
         #yplot = np.log10(num_den_passive_we_mass[ind,i])
         #xplot = zlist[ind]
         #ax.plot(xplot,yplot[0],ls=lines[1], color=cols[i], label=labels[i] + ' with error')

    plot_observations_num_density(ax, labels=False)
    load_sharkv2p0_variations(ax, labels=True)

    common.prepare_legend(ax, ['r', 'r', 'LightCoral', 'LightCoral'], loc=3)
    plt.tight_layout()
    common.savefig(outdir, fig, 'num_density_passive_modelvariations.pdf')


def prepare_data(hdf5_data, index, num_den_passive, num_den_passive_we, num_den_passive_mass, num_den_passive_we_mass, redshift, num_dens_halos):


    # Unpack data
    (h0, volh, typeg, mdisk, mbulge, sfrd, sfrb, rdisk, mBH, mHI, mH2, mgas,
     mHI_bulge, mH2_bulge, mgas_bulge, mvir) = hdf5_data
       
    #look at number densities of galaxies with sSFR<1e-10yr^-1
    ms_tot = np.log10((mdisk+mbulge)/h0)
    mvir = np.log10(mvir/h0)

    ind = np.where(ms_tot > 10.8)
    ngals = len(ms_tot[ind])
    num_dens = ngals / (volh / h0**3.0)
    print("number density of galaxies with m>10.8", num_dens, "at redshift", redshift)

    

    sfr_tot = np.log10((sfrd + sfrb)/h0/1e9)
    ssfr = sfr_tot - ms_tot
    ind = np.where((ms_tot > 10) & (ssfr < -10))
    npass = len(ms_tot[ind])
    num_den_passive[index,0] = (npass + 0.0) / (volh / h0**3)
    ind = np.where((ms_tot > 10) & (ssfr < -11))
    npass = len(ms_tot[ind])
    num_den_passive[index,1] = (npass + 0.0) / (volh / h0**3)

    ind = np.where((ms_tot > 10.5) & (ssfr < -10))
    npass = len(ms_tot[ind])
    num_den_passive_mass[index,0] = (npass + 0.0) / (volh / h0**3)
    ind = np.where((ms_tot > 10.5) & (ssfr < -11))
    npass = len(ms_tot[ind])
    num_den_passive_mass[index,1] = (npass + 0.0) / (volh / h0**3)

    #add random errors to sfr and mstar
    ms_tot_we = ms_tot + np.random.normal(0.0, 0.3, len(ms_tot)) 
    ind = np.where(ms_tot_we > 10.8)
    ngals = len(ms_tot_we[ind])
    num_dens = ngals / (volh / h0**3.0)
    print("number density of galaxies with m>10.8 (after adding errors)", num_dens, "at redshift", redshift)


    ssfr_we = sfr_tot + np.random.normal(0.0, 0.3, len(ms_tot)) - ms_tot_we
    ind = np.where((ms_tot_we > 10) & (ssfr_we < -10))
    npass = len(ms_tot[ind])
    num_den_passive_we[index,0] = (npass + 0.0) / (volh / h0**3)

    ind = np.where((ms_tot_we > 10) & (ssfr_we < -11))
    npass = len(ms_tot[ind])
    num_den_passive_we[index,1] = (npass + 0.0) / (volh / h0**3)

    ind = np.where((ms_tot_we > 10.6) & (ssfr_we < -10))
    npass = len(ms_tot[ind])
    num_den_passive_we_mass[index,0] = (npass + 0.0) / (volh / h0**3)

    ind = np.where((ms_tot_we > 10.6) & (ssfr_we < -11))
    npass = len(ms_tot[ind])
    num_den_passive_we_mass[index,1] = (npass + 0.0) / (volh / h0**3)

    ind = np.where(num_den_passive_we_mass == 0)
    num_den_passive_we_mass[ind] = 0.5 / (volh / h0**3)
    ind = np.where(num_den_passive_mass == 0)
    num_den_passive_mass[ind] = 0.5 / (volh / h0**3)
    ind = np.where(num_den_passive == 0)
    num_den_passive[ind] = 0.5 / (volh / h0**3)
    ind = np.where(num_den_passive_we == 0)
    num_den_passive_we[ind] = 0.5 / (volh / h0**3)

    mhalo_bins = [11.5, 11.6,11.7,11.8,11.9,12,12.1,12.2]
    
    for j, b in enumerate(mhalo_bins):
        ind = np.where((mvir > b) & (typeg == 0))
        ngals = len(mvir[ind])
        num_dens_halos[index, j] = (ngals +0.0) / (volh / h0**3.0)

def main(model_dir, output_dir, redshift_table, subvols, obs_dir):


    plt = common.load_matplotlib()

    zlist = np.array([0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]) #, 8.0, 9.0, 10.0, 11.0, 12.0])
    #zlist = np.array([5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0])
    fields = {'galaxies': ('type', 'mstars_disk', 'mstars_bulge', 'sfr_disk', 'sfr_burst',
                           'rstar_disk', 'm_bh', 'matom_disk', 'mmol_disk', 'mgas_disk',
                           'matom_bulge', 'mmol_bulge', 'mgas_bulge', 'mvir_hosthalo')}

    num_den_passive = np.zeros(shape = (len(zlist), 2))
    num_den_passive_we = np.zeros(shape = (len(zlist), 2))
    num_den_passive_mass = np.zeros(shape = (len(zlist), 2))
    num_den_passive_we_mass = np.zeros(shape = (len(zlist), 2))
    num_dens_halos = np.zeros(shape = (len(zlist), 8))
    for index, snapshot in enumerate(redshift_table[zlist]):

        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)

        prepare_data(hdf5_data, index, num_den_passive, num_den_passive_we, num_den_passive_mass, num_den_passive_we_mass, zlist[index], num_dens_halos)

    #for a,b,c,d,e in zip(zlist, num_den_passive[:,0], num_den_passive[:,1], num_den_passive_mass[:,0], num_den_passive_mass[:,1]):
    #    print(a,b,c,d,e)
    plot_num_density_passive(plt, output_dir, obs_dir, zlist, num_den_passive, num_den_passive_we, num_den_passive_mass, num_den_passive_we_mass)

    #for i in range(0,len(zlist)):
    #    print(zlist[i],num_dens_halos[i,0], num_dens_halos[i,1], num_dens_halos[i,2], num_dens_halos[i,3], num_dens_halos[i,4], num_dens_halos[i,5], num_dens_halos[i,6],num_dens_halos[i,7])

if __name__ == '__main__':
    main(*common.parse_args())
