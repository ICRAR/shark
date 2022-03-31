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
"""Size plots"""

import functools

import numpy as np

import os
import common
import utilities_statistics as us

# Initialize arguments
zlist = [0]

##################################
#Constants
RExp     = 1.67
MpcToKpc = 1e3
G        = 4.299e-9 #Gravity constant in units of (km/s)^2 * Mpc/Msun

mlow = 6.5
mupp = 12.5
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

def prepare_data(hdf5_data, age_bins, sfh, delta_t, LBT, index, disk_size, disk_size_age):

    (h0, _, mdisk, mbulge, rdisk, typeg, age) = hdf5_data
    (bulge_diskins_hist, bulge_mergers_hist, disk_hist) = sfh

    bin_it   = functools.partial(us.wmedians, xbins=xmf)

    BTthresh = 0.5
 
    mstar = mdisk + mbulge
    ind = np.where(mstar > 0)
    mdisk = mdisk[ind]
    mbulge = mbulge[ind]
    rdisk = rdisk[ind]
    typeg = typeg[ind]
    age = age[ind]
    mstar = mstar[ind]

    age_disk = np.zeros(shape = len(mdisk))
    for j in range(0,len(mdisk)):
          sfrtimesage = np.sum(disk_hist[j,:] * delta_t[:] * LBT[:])
          massformed = np.sum(disk_hist[j,:] * delta_t[:])
          age_disk[j] = sfrtimesage / massformed

    ind = np.where((mdisk > 0) & (mdisk / mstar >= BTthresh))
    disk_size[index,:] = bin_it(x=np.log10(mdisk[ind]) - np.log10(float(h0)),
                                y=np.log10(rdisk[ind]*MpcToKpc) - np.log10(float(h0)))

    BtoT = mbulge / mstar
    ind = np.where((mdisk > 0))
    print("#log10(mdisk/Msun) log10(rdisk/kpc) age_disk age_gal B/T")
    for a,b,c,d,e in zip(np.log10(mdisk[ind]) - np.log10(float(h0)), np.log10(rdisk[ind]*MpcToKpc) - np.log10(float(h0)), age_disk[ind], 13.7969 - age[ind], BtoT[ind]):
        print(a,b,c,d,e)

    age_disk = 13.7969 - age 
    age_meds = np.zeros(shape = len(age_bins)-1)
    for j in range(0,len(age_bins)-1):
        ind = np.where((mdisk > 0) & (age_disk >= age_bins[j]) & (age_disk < age_bins[j+1]) & (mdisk / mstar >= BTthresh))
        disk_size_age[index,j,:] = bin_it(x=np.log10(mdisk[ind]) - np.log10(float(h0)),
                                y=np.log10(rdisk[ind]*MpcToKpc) - np.log10(float(h0)))
        age_meds[j] = np.mean(age_disk[ind])

    ind = np.where((mdisk >= 5e7) & (mdisk / mstar >= BTthresh))
    return(np.log10(mdisk[ind]) - np.log10(float(h0)), np.log10(rdisk[ind]*MpcToKpc) - np.log10(float(h0)), age_disk[ind], age_meds)
 
def plot_age_disk(plt, outdir, obsdir, mdisk_z0, rdisk_z0, age_z0, disk_size, disk_size_age, age_bins, age_meds):

    fig = plt.figure(figsize=(5,4.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm disk}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm r_{\\star,disk}/kpc)$"
    xmin, xmax, ymin, ymax = 8, 12, -0.5, 2

    # LTG ##################################
    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    im = ax.hexbin(mdisk_z0, rdisk_z0, age_z0, xscale='linear', yscale='linear', gridsize=(20,20), cmap='magma', mincnt=30)
    #cbar_ax = fig.add_axes([0.86, 0.15, 0.025, 0.7])
    cbar = fig.colorbar(im)#, cax=cbar_ax)
    cbar.ax.set_ylabel('mean stellar age')

    #Predicted size-mass for disks in disk=dominated galaxies
    ind = np.where(disk_size[0,0,:] != 0)
    xplot = xmf[ind]
    yplot = disk_size[0,0,ind]
    errdn = disk_size[0,1,ind]
    errup = disk_size[0,2,ind]
    ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'k', mec='k',marker='o',label="Shark")

    def plot_robotham_size_relation_use_bins(ax, age_meds):
        for age in age_meds:
            y1 = 0.35 * (xmf - 10) - 0.06 * (age) + 1.01
            ax.plot(xmf, y1, linestyle='dotted', color='k')


    def plot_lange_relation(ax):
        #Lange et al. (2016)
        a = 5.56
        aerr = 1.745
        b = 0.274
        ind = np.where(xmf < 10.3)
        rL16 = np.log10(a*pow((pow(10.0,xmf)/1e10),b))
        ax.plot(xmf[ind],rL16[ind],'b', linestyle='solid', label ='L16 disks')
    
        m,r = common.load_observation(obsdir, 'SizesAndAM/rdisk_L16.dat', [0,1])
        ax.plot(m[0:36], r[0:36], linestyle='dotted',color='b',label="50th, 68th, 90th")
        ax.plot(m[38:83], r[38:83], linestyle='dotted',color='b')
        ax.plot(m[85:128], r[85:129], linestyle='dotted',color='b')

    plot_robotham_size_relation_use_bins(ax, np.array([3, 6, 9, 12]))
    plot_lange_relation(ax)

    common.prepare_legend(ax, ['b','b','k','r'], loc=2)
    common.savefig(outdir, fig, 'age_disks.pdf')


    fig = plt.figure(figsize=(5,4.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm disk}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm r_{\\star,disk}/kpc)$"
    xmin, xmax, ymin, ymax = 8, 12, -0.5, 2

    # LTG ##################################
    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    for j in range(0,len(age_bins)-1):
        ind = np.where(disk_size_age[0,j,0,:] != 0)
        xplot = xmf[ind]
        yplot = disk_size_age[0,j,0,ind]
        errdn = disk_size_age[0,j,1,ind]
        errup = disk_size_age[0,j,2,ind]
        ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', marker='o', label='age=[%s,%s]' % (str(age_bins[j]), str(age_bins[j+1])))

    plot_robotham_size_relation_use_bins(ax, age_meds)
    common.prepare_legend(ax, ['k','k','k','k','k','k','k'], loc=2)
    common.savefig(outdir, fig, 'age_bins_disks.pdf')


def main(modeldir, outdir, redshift_table, subvols, obsdir):

    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'rstar_disk', 'type', 
                           'mean_stellar_age')}

    sfh_fields = {'bulges_diskins': ('star_formation_rate_histories'),
                  'bulges_mergers': ('star_formation_rate_histories'),
                  'disks': ('star_formation_rate_histories')}
 
    # Loop over redshift and subvolumes

    age_bins=[0,3,6,10,14]
    disk_size = np.zeros(shape = (len(zlist), 3, len(xmf)))
    disk_size_age = np.zeros(shape = (len(zlist), len(age_bins)-1, 3, len(xmf)))

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        sfh, delta_t, LBT = common.read_sfh(modeldir, 199, sfh_fields, subvols)

        (mdisk, rdisk, age, age_meds) = prepare_data(hdf5_data, age_bins, sfh, delta_t, LBT, index, disk_size, disk_size_age)
        if(index == 0):
           mdisk_z0 = mdisk
           rdisk_z0 = rdisk
           age_z0 = age
           age_meds_z0 = age_meds

    print("Will produce plots")
    plot_age_disk(plt, outdir, obsdir, mdisk_z0, rdisk_z0, age_z0, disk_size, disk_size_age, age_bins, age_meds)

if __name__ == '__main__':
    main(*common.parse_args())
