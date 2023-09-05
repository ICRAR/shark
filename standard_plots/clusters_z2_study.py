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

import common
import utilities_statistics as us

##################################
# Constants
mlow = 7.0
mupp = 13.0
dm = 0.5
mbins = np.arange(mlow, mupp, dm)
xmf = mbins + dm/2.0

#rbins = np.array([1, 2.730045591676135, 5]) #Mpc/h
rbins = np.array([1, 2.8315841879187973, 5]) #Mpc/h
zdepth = 40.30959350543804 #Mpc/h

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

def prepare_data(hdf5_data, hdf5_data_co, index, lfco, rhoh2, rhoco, rhosfr):

    bin_it = functools.partial(us.wmedians, xbins=xmf)
    stack  = functools.partial(us.stacking, xbins=xmf)

    # Unpack data
    (h0, _, typeg, mdisk, mbulge, _, _, mHI, mH2, mgas,
     mHI_bulge, mH2_bulge, mgas_bulge, mvir, sfrd, sfrb, 
     x, y, z, vvir) = hdf5_data

    (lco_disk, lco_bulge) = hdf5_data_co
    co_total = (lco_disk[:,0] + lco_bulge[:,0]) * 3.25e7 / (115.2712)**2.0 / (4.0*PI) #in K km/s pc^2

    ind = np.where( (mdisk +  mbulge) > 0)
    mvir = mvir[ind]
    vvir = vvir[ind]
    mH2 = mH2[ind]
    mH2_bulge = mH2_bulge[ind]
    x = x[ind]
    y = y[ind]
    z = z[ind]
    mdisk = mdisk[ind]
    mbulge = mbulge[ind]
    sfrd = sfrd[ind]
    sfrb = sfrb[ind]
    typeg = typeg[ind]


    XH = 0.72
    h0log = np.log10(float(h0))

    rvir  = G * mvir / pow(vvir,2.0) / h0

    mstar_tot = (mdisk + mbulge) / h0
    sfr_tot = (sfrd + sfrb) / h0 / GyrtoYr

    #define main sequence first
    indcen = np.where((mvir/h0 > 1e12) & (typeg == 0))
    mvir_in = mvir[indcen]/h0
    sortedmass = np.argsort(1.0/mvir_in)
    mass_ordered = mvir_in[sortedmass]
    mass_tresh  = mass_ordered[9]
    print(mass_tresh)
    indcen = np.where((mvir/h0 >= mass_tresh) & (typeg == 0))
    x_cen = x[indcen]
    y_cen = y[indcen]
    z_cen = z[indcen]
    rvir_cen = rvir[indcen]

    lfco_ind_this_z = np.zeros(shape = (len(rbins), len(x_cen), len(mbins)))
    codensity_ind_this_z = np.zeros(shape = (len(rbins), len(x_cen)))
    h2density_ind_this_z = np.zeros(shape = (len(rbins), len(x_cen)))
    sfrdensity_ind_this_z = np.zeros(shape = (len(rbins), len(x_cen)))

    print ('number of clusters %d'% len(x_cen))
    #xyz distance
    for g in range(0,len(x_cen)):
        d_all = np.sqrt(pow(x - x_cen[g], 2.0) + pow(y - y_cen[g], 2.0))
        z_all = np.absolute(z - z_cen[g])
        ind   = np.where((d_all < 10) & (co_total > 0) & (z_all < zdepth))
        co_total_in = co_total[ind]
        mh2total_in = (mH2[ind] + mH2_bulge[ind]) / h0 * XH
        sfrtotal_in = sfr_tot[ind]
        d_all_in = d_all[ind]
        for r in range(0, len(rbins)):
            vol = PI * (rbins[r]/h0)**2.0 * (2.0 * zdepth/h0)
            inr = (d_all_in < rbins[r])
            H, _ = np.histogram(np.log10(co_total_in[inr]),bins=np.append(mbins,mupp))
            lfco_ind_this_z[r,g,:] = lfco_ind_this_z[r,g,:] + H 
            #divide by the volume of the cluster
            pos = np.where(lfco_ind_this_z[r,g,:] > 0)
            lfco_ind_this_z[r,g,pos] = np.log10(lfco_ind_this_z[r,g,pos]/vol/dm)
            codensity_ind_this_z[r,g] = np.log10(np.sum(co_total_in[inr])/vol)
            h2density_ind_this_z[r,g] = np.log10(np.sum(mh2total_in[inr])/vol)
            sfrdensity_ind_this_z[r,g] = np.log10(np.sum(sfrtotal_in[inr])/vol)


    for r in range(0, len(rbins)):
        rhoco[:,index,r] = us.gpercentiles(x=codensity_ind_this_z[r,:])
        rhoh2[:,index,r] = us.gpercentiles(x=h2density_ind_this_z[r,:])
        rhosfr[:,index,r] = us.gpercentiles(x=sfrdensity_ind_this_z[r,:])
 
        for b in range(0, len(mbins)):
            lfco[:,index,b,r] = us.gpercentiles(x=lfco_ind_this_z[r,:,b])

    return (lfco_ind_this_z, codensity_ind_this_z)

def plot_density_radii(plt, output_dir, zlist, rhoh2, rhoco, rhosfr):

    ###################################
    #   Plots global mass densities
    fig = plt.figure(figsize=(6,8))
    plt.subplots_adjust(bottom=0.15, left=0.15)

    xmin, xmax, ymin, ymax = 0.5, 5, 6, 13
    xtitle = [' ', '$\\rm redshift$']
    ytitle = ['$\\rm log_{10}(\\rho_{\\rm LCO(1-0)} [K\\, km/s\\, pc^{2}\\, Mpc^{-3}])$', '$\\rm log_{10}(\\rho_{\\rm H_2} [M_{\\odot}\\, Mpc^{-3}])$', '$\\rm log_{10}(\\rho_{\\rm SFR} [M_{\\odot}\\, yr^{-1}\\, Mpc^{-3}]$)']
    xleg = xmax - 0.1 * (xmax - xmin)
    yleg = ymin + 0.1 * (ymax - ymin)

    cols = ('Crimson','Green','blue')
    labels = ('$r<2.5\\, \\rm cMpc$', '$r<5\\, \\rm cMpc$','$r<7.5\\, \\rm cMpc}$')

    ax = fig.add_subplot(311)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitle[0], ytitle[0], locators=(1, 1, 1))

    #ax.text(xleg, yleg, '$\\rm M_{\\rm halo} > 10^{14} M_{\odot}$')

    print('will print CO density')
    for r in range(0, len(rbins)):
        for a,b,c,d in zip(zlist,  rhoco[0,:,r], rhoco[1,:,r], rhoco[2,:,r]):
            print(a,b,c,d)

        ax.errorbar(zlist, rhoco[0,:,r], yerr=[rhoco[1,:,r], rhoco[2,:,r]], ls='None', mfc='None', ecolor =cols[r],  mec=cols[r],  marker='o', label=labels[r])
    common.prepare_legend(ax, cols, loc = 4)

    ax = fig.add_subplot(312)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitle[0], ytitle[1], locators=(1, 1, 1))
    print('will print H2 density')
    for r in range(0, len(rbins)):
        for a,b,c,d in zip(zlist,  rhoh2[0,:,r], rhoh2[1,:,r], rhoh2[2,:,r]):
            print(a,b,c,d)

        ax.errorbar(zlist, rhoh2[0,:,r], yerr=[rhoh2[1,:,r], rhoh2[2,:,r]], ls='None', mfc='None', ecolor =cols[r],  mec=cols[r],  marker='o')

    ax = fig.add_subplot(313)
    common.prepare_ax(ax, xmin, xmax, -2, 1, xtitle[1], ytitle[2], locators=(1, 1, 1))
    print('will print SFR density')
    for r in range(0, len(rbins)):
        for a,b,c,d in zip(zlist,  rhosfr[0,:,r], rhosfr[1,:,r], rhosfr[2,:,r]):
            print(a,b,c,d)

        ax.errorbar(zlist, rhosfr[0,:,r], yerr=[rhosfr[1,:,r], rhosfr[2,:,r]], ls='None', mfc='None', ecolor =cols[r],  mec=cols[r],  marker='o')

    common.savefig(output_dir, fig, "cluster_H2CO10density.pdf")

def plot_clusters_lco(plt, output_dir, zlist, lfco):

    ###################################
    #   Plots global mass densities
    fig = plt.figure(figsize=(6,7))
    plt.subplots_adjust(bottom=0.15, left=0.15)
    xmin, xmax, ymin, ymax = 6, 12, -4, 1
    xleg = xmax - 1 * (xmax - xmin)
    yleg = ymin + 0.1 * (ymax - ymin)

    xtitle = '$\\rm log_{10}(\\rm L^{\\prime}_{\\rm CO(1-0)}[K\\, km/s\\, pc^2])$'
    ytitle = '$\\rm log_{10}(\\rm dn [Mpc^{-3}\\, dex^{-1}])$'

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitle, ytitle, locators=(1, 1, 1))
    #ax.text(xleg, yleg, '$\\rm M_{\\rm halo} > 10^{14} M_{\odot}$')

    r = 1
    cols = ['b', 'DarkTurquoise', 'Orange', 'DarkRed']
    for z in range(0, len(zlist)):
        ind = np.where(lfco[0,z,:,r] != 0)
        yplot = lfco[0,z,ind,r]
        yerrdn= lfco[1,z,ind,r]
        yerrup= lfco[2,z,ind,r]
        print('will co lf at z=%s' % str(zlist[z]))
        for a,b,c,d in zip(xmf[ind], yplot[0], yerrdn[0], yerrup[0]):
            print(a,b,c,d)
        #print(yerrdn[0])
        #common.errorbars(ax, xmf[ind], yplot[0], yerrdn[0], yerrup[0], cols[z], 'o', label='z=%s' % str(zlist[z]))
        ax.errorbar(xmf[ind], yplot[0], yerr=[yerrdn[0], yerrup[0]], ls='None', mfc='None', ecolor = cols[z], mec=cols[z], marker='o', label='z=%s' % str(zlist[z]))
    common.prepare_legend(ax, cols, loc = 3)
    common.savefig(output_dir, fig, "co_lf_clusters.pdf")

def main(model_dir, output_dir, redshift_table, subvols, obs_dir):


    plt = common.load_matplotlib()

    zlist = (0.5, 1, 1.5, 2, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0)

    fields = {'galaxies': ('type', 'mstars_disk', 'mstars_bulge',
                           'rstar_disk', 'm_bh', 'matom_disk', 'mmol_disk', 'mgas_disk',
                           'matom_bulge', 'mmol_bulge', 'mgas_bulge', 'mvir_hosthalo', 'sfr_disk',
                           'sfr_burst', 'position_x', 'position_y', 'position_z', 'vvir_hosthalo')}

    fields_co = {'galaxies': ('LCO_disk','LCO_bulge')}

    lfco = np.zeros(shape = (3, len(zlist), len(mbins), len(rbins)))
    rhoco = np.zeros(shape = (3, len(zlist), len(rbins)))
    rhoh2 = np.zeros(shape = (3, len(zlist), len(rbins)))
    rhosfr = np.zeros(shape = (3, len(zlist), len(rbins)))

    for index, snapshot in enumerate(redshift_table[zlist]):

        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        hdf5_data_co = common.read_co_data(model_dir, snapshot, fields_co, subvols)
        (lfco_ind_this_z, codensity_ind_this_z) = prepare_data(hdf5_data, hdf5_data_co, index, lfco, rhoh2, rhoco, rhosfr)

    plot_density_radii(plt, output_dir, zlist, rhoh2, rhoco, rhosfr)
    plot_clusters_lco(plt, output_dir, zlist, lfco)

if __name__ == '__main__':
    main(*common.parse_args())
