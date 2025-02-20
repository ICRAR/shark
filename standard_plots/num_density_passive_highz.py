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
mupp = 15.0
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

def plot_num_density_passive(plt, outdir, obs_dir, zlist, num_densities, ssfr_thresh, mass_threshs, mass_threshs2, num_densities2, limit):

    def plot_observations_num_density(ax, labels=True):
    
        z, zerr, n, nerrh, nerrl = np.loadtxt(obs_dir + '/num_densities/' + 'Gould23data.txt', unpack = True, usecols = [0,1,2,3,4])
        obs = ['Straatman+14', 'Gould+23', 'Carnall+20', 'Schreiber+18', 'Shahidi+20', 'Merlin+19', 'Girelli+19', 'Weaver+23', 'Carnall+23', 'Valentino+23', 'Nanayakkara+24', 'Alberts+24']
        nerrh = np.log10(n+nerrh) - np.log10(n)
        nerrl = np.log10(n) - np.log10(n-nerrl)
        n = np.log10(n)
        Nlow = [0, 1, 4, 8, 9, 11, 13, 17, 22, 24, 28, 29]
        lenght = [1, 3, 4, 1, 2, 2, 4, 5, 2, 4, 1, 2]
        symbols = ['o', 's', '*', 'd', 'D', '<', '>' , 'p', 'h', 'P', 'X', 'H']
        for j in range(0,len(obs)):
            if(j <= 7):
                col = 'grey'
            else:
                col = 'darkgreen'
            if(labels):
               ax.errorbar(z[Nlow[j]:Nlow[j]+lenght[j]], n[Nlow[j]:Nlow[j]+lenght[j]], xerr=zerr[Nlow[j]:Nlow[j]+lenght[j]], yerr=[nerrl[Nlow[j]:Nlow[j]+lenght[j]], nerrh[Nlow[j]:Nlow[j]+lenght[j]]], color=col, ecolor=col, marker=symbols[j], ls='None', mfc='None', label=obs[j])
            else:
                ax.errorbar(z[Nlow[j]:Nlow[j]+lenght[j]], n[Nlow[j]:Nlow[j]+lenght[j]], xerr=zerr[Nlow[j]:Nlow[j]+lenght[j]], yerr=[nerrl[Nlow[j]:Nlow[j]+lenght[j]], nerrh[Nlow[j]:Nlow[j]+lenght[j]]], color=col, ecolor=col, marker=symbols[j], ls='None', mfc='None')

    fig = plt.figure(figsize=(7,6))
    ytit = "$\\rm log_{10} (\\rm \\phi/Mpc^{-3})$"
    xtit = "$\\rm redshift$"
    xmin, xmax, ymin, ymax = 2, 5, -7.5, -3.5
    xleg = xmax - 0.3 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)
   
    ax = fig.add_subplot(111)
    ax2 = fig.add_subplot(111, facecolor="none")
   
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    common.prepare_ax(ax2, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
   
    ax.plot([xleg-0.3,xleg-0.05], [yleg+0.15,yleg+0.15], ls='solid',color='DarkOrange', lw=4)
    ax.plot([xleg-0.3,xleg-0.05], [yleg-0.1,yleg-0.1], ls='solid',color='DarkOrange', lw=1)
   
   
    ax.text(xleg,yleg+0.15,'$\\rm M_{\\star}>10^{10}\,M_{\\odot}$', fontsize=13, color='DarkOrange')
    ax.text(xleg,yleg-0.1,'$\\rm M_{\\star}>10^{10.5}\,M_{\\odot}$', fontsize=13, color='DarkOrange')
    ax.text(3.2, -3.7,'SHARK', fontsize=16)
 
    lines = ['solid','dashed','dotted']
    lws = [4,1]
   
    for i in range(0,len(ssfr_thresh)):
        for g  in range(0,len(mass_threshs)):
            #Predicted relation
            ind = np.where(num_densities[:,i,g] != 0)
            yplot = np.log10(num_densities[ind,i,g])
            xplot = zlist[ind]
            if(g == 0):
               ax.plot(xplot,yplot[0],ls=lines[i], color='DarkOrange', lw=lws[g], label="sSFR<%s$\\rm yr^{-1}$" % str(ssfr_thresh[i]))
            else:
               ax.plot(xplot,yplot[0],ls=lines[i], color='DarkOrange', lw=lws[g], label=None)
   
    plot_observations_num_density(ax2, labels=True)
    common.prepare_legend(ax, ['DarkOrange', 'DarkOrange', 'DarkOrange'], loc=4)
    common.prepare_legend(ax2, ['grey','grey', 'grey', 'grey', 'grey','grey', 'grey', 'grey', 'darkgreen','darkgreen', 'darkgreen','darkgreen'], loc=3)

    plt.tight_layout()
    common.savefig(outdir, fig, 'num_density_passive_highzonly.pdf')

    fig = plt.figure(figsize=(7,6))
    ytit = "$\\rm log_{10} (\\rm \\phi/Mpc^{-3})$"
    xtit = "$\\rm redshift$"
    xmin, xmax, ymin, ymax = 2, 7, -8.6, -3.5
    xleg = xmax - 0.3 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)
   
    ax = fig.add_subplot(111)
    ax2 = fig.add_subplot(111, facecolor="none")
   
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    common.prepare_ax(ax2, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
   
    ax.text(3.2, -3.7,'SHARK', fontsize=16)
 
    lines = ['solid','dashed','dotted']
    lws = [4,1]
   
    i = 0
    g = 0
    #Predicted relation
    ind = np.where(num_densities[:,i,g] != 0)
    yplot = np.log10(num_densities[ind,i,g])
    xplot = zlist[ind]
    ax.plot(xplot,yplot[0],ls=lines[2], color='DarkOrange', lw=lws[g], label=None)
    for a,b in zip(xplot,yplot[0]):
        print(a,b) 
    ind = np.where((num_densities[:,i,g] != 0) & (num_densities[:,i,g] > limit))
    yplot = np.log10(num_densities[ind,i,g])
    xplot = zlist[ind]
    ax.plot(xplot,yplot[0],ls=lines[i], color='DarkOrange', lw=lws[g], label=None)
    for a,b in zip(xplot,yplot[0]):
        print(a,b)   

    plot_observations_num_density(ax2, labels=True)
    common.prepare_legend(ax, ['DarkOrange', 'DarkOrange', 'DarkOrange'], loc=4)
    common.prepare_legend(ax2, ['grey','grey', 'grey', 'grey', 'grey','grey', 'grey', 'grey', 'darkgreen','darkgreen', 'darkgreen','darkgreen'], loc=3)

    plt.tight_layout()
    common.savefig(outdir, fig, 'num_density_passive_highzonly_v2.pdf')


def plot_mvir_final(plt, output_dir, obs_dir, zlist, halo_hists):

    fig = plt.figure(figsize=(7,6))
    xtit = "$\\rm log_{10} (M_{\\rm halo}(z=0)/M_{\\odot})$"
    ytit = "$\\rm pdf$"
    xmin, xmax, ymin, ymax = 12, 15, 0, 1
    xleg = xmax - 0.3 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)
  
    subplots = [221, 222, 223, 224]
    zint = [2, 3, 4, 5]
    color= 'DarkOrange'

    for i in range(0,len(zint)):
        ax = fig.add_subplot(subplots[i])
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1,1 , 0.25, 0.25))
        ax.text(14.4, 0.9, 'z=%s' % str(zint[i]))
        #select correct redshift
        j = np.where(zlist == zint[i])
        hist_int = halo_hists[j,:][0]
        print(hist_int.shape)
        hist_all = hist_int[:,0,:] / (np.sum(hist_int[:,0,:])*dm)
        hist_pass = hist_int[:,1,:] / (np.sum(hist_int[:,1,:])*dm)
        print(hist_pass, hist_pass.shape)
        ax.plot(xmf, hist_all[0], ls='dashed',color='DarkOrange', lw=4, label = 'all massive')
        ax.plot(xmf, hist_pass[0], ls='solid',color='DarkOrange', lw=4, label = 'massive-passive')
        common.prepare_legend(ax, ['DarkOrange', 'DarkOrange', 'DarkOrange'], loc=2)

    plt.tight_layout()
    common.savefig(output_dir, fig, 'final_halo_mass_passivegals.pdf')


def prepare_data(hdf5_data, hdf5_data_halo, index, num_densities, ssfr_thresh, mass_threshs, mass_threshs2, num_densities2, halo_hists, zlist):


    # Unpack data
    (h0, volh, typeg, mdisk, mbulge, sfrd, sfrb, idhalo, mvir, mbh) = hdf5_data

    (_, _, mvirz0, idhalo_cat) = hdf5_data_halo

    vol = volh/h0**3

    #look at number densities of galaxies with sSFR<1e-10yr^-1
    ms_tot = ((mdisk+mbulge)/h0)
    sfr_tot = ((sfrd + sfrb)/h0/1e9)

    #ind = np.where(ms_tot >= 1e10)
    #allgals  = np.zeros(shape = (len(ms_tot[ind]), 5))
    #allgals[:,0] = ms_tot[ind]
    #allgals[:,1] = sfr_tot[ind]
    #allgals[:,2] = mbh[ind]/h0
    #allgals[:,3] = mvir[ind]/h0
    #allgals[:,4] = typeg[ind]
    #np.savetxt("AllGalaxiesz3_withhalo.txt", allgals)

    #ids_interest = idhalo[ind]
    #mvir_final = np.zeros(shape = len(ms_tot[ind]))
    #for g in range(0,len(ms_tot[ind])):
    #    idh = np.where(idhalo_cat == ids_interest[g])
    #    mvir_final[g] = mvirz0[idh]

    #ms_tot = ms_tot[ind]
    #sfr_tot = sfr_tot[ind]

    #H, _ = np.histogram(np.log10(mvir_final/h0),bins=np.append(mbins,mupp))
    #halo_hists[index,0,:] = halo_hists[index,0,:] + H

    #ind = np.where((sfr_tot/ms_tot <= 1e-10) & (ms_tot > 1e10))
    #print("Number density of passive galaxies with M>1e10Msun", (len(sfr_tot[ind]) + 0.0)/vol, " at ", index)
    #H, _ = np.histogram(np.log10(mvir_final[ind]/h0),bins=np.append(mbins,mupp))
    #halo_hists[index,1,:] = halo_hists[index,1,:] + H
    ind = np.where(typeg ==0)
    mvirin = mvir[ind]
    mvir_ordered = mvirin[np.argsort(1/mvirin)]/h0
    ind = np.where((sfr_tot/ms_tot <= 1e-10) & (ms_tot > 1e10))
    print("Median halo mass", np.median(np.log10(mvir[ind]/h0)), " at ", zlist[index], " while the 100 most massive halo in this redshift have a median", np.log10(np.median(mvir_ordered[0:99])))
    ind = np.where(ms_tot >= 1e10)
    print("Number density of all massive galaxies, ", np.log10(len(ms_tot[ind])/vol), " at redshift", zlist[index])

    for j in range(0,len(ssfr_thresh)):
        for g in range(0,len(mass_threshs)):
            passive = np.where(( ms_tot>=mass_threshs[g]) & (sfr_tot/ms_tot <= ssfr_thresh[j]))
            npass = len(ms_tot[passive])
            num_densities[index,j,g] = (npass + 0.0) / vol
    ind = np.where(num_densities == 0)
    num_densities[ind] = (0.99)/vol

    scatter = np.random.normal(0.0, 0.25, len(ms_tot))
    ms_tot_err = ms_tot + scatter
    for g in range(0,len(mass_threshs2)):
        massi = np.where((ms_tot >=mass_threshs2[g]) & (sfr_tot/ms_tot <= 1e-10))
        npass = len(ms_tot[massi])
        num_densities2[index,g] = (npass + 0.0) / vol

    return (0.99)/vol

def main(model_dir, output_dir, redshift_table, subvols, obs_dir):


    plt = common.load_matplotlib()

    zlist = np.array([2, 2.5, 3.0, 3.53362989, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]) #, 8.0, 9.0, 10.0, 11.0, 12.0])
    #zlist = np.array([5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0])
    fields = {'galaxies': ('type', 'mstars_disk', 'mstars_bulge', 'sfr_disk', 'sfr_burst', 'id_halo_tree', 'mvir_hosthalo', 'm_bh')}

    fields_halo = {'halo': ('final_z0_mvir', 'halo_id')}
    ssfr_thresh = [1e-10, 0.5e-10, 1e-11]
    mass_threshs = [1e10, 3e10]
    mass_threshs2 = [10**10.9, 1e11]
    mass_threshs2 = [1e10, 10**10.3, 10**10.6]

    num_densities = np.zeros(shape = (len(zlist),len(ssfr_thresh),len(mass_threshs)))
    num_densities2 = np.zeros(shape = (len(zlist),len(mass_threshs2)))
    halo_hists = np.zeros(shape = (len(zlist), 2, len(xmf)))

    for index, snapshot in enumerate(redshift_table[zlist]):

        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        hdf5_data_halo = common.read_data(model_dir, snapshot, fields_halo, subvols)

        limit = prepare_data(hdf5_data, hdf5_data_halo, index, num_densities, ssfr_thresh, mass_threshs, mass_threshs2, num_densities2, halo_hists, zlist)
    #for a,b,c,d in zip(zlist, num_densities2[:,0], num_densities2[:,1], num_densities2[:,2]):
    #    print(a,b,c,d)

    plot_num_density_passive(plt, output_dir, obs_dir, zlist, num_densities, ssfr_thresh, mass_threshs, mass_threshs2, num_densities2, limit)
    plot_mvir_final(plt, output_dir, obs_dir, zlist, halo_hists)

if __name__ == '__main__':
    main(*common.parse_args())
