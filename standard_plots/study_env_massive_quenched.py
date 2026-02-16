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


slow = -3
supp = 3
dsfr = 0.125
sfrbins = np.arange(slow,supp,dsfr)
xsfr = sfrbins + dsfr/2.0

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


def plot_mass_sfrs_neighbours(plt, outdir, obsdir, npass, mstars_neighbours, sfr_neighbours, deltav_neighbours):

    fig = plt.figure(figsize=(9,8))
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm <N>$"
    xmin, xmax, ymin, ymax = 7, 11, 0.9, 150

    ax = fig.add_subplot(221)
    ax.set_yscale('log')
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.5, 0.5, 100, 100))

    histm = np.histogram(np.log10(mstars_neighbours),bins=np.append(mbins,mupp))[0]/npass
    ax.plot(xmf,histm, linestyle='solid', color='black', label='all')

    ind = np.where(abs(deltav_neighbours) < 1000)
    histm = np.histogram(np.log10(mstars_neighbours[ind]),bins=np.append(mbins,mupp))[0]/npass
    ax.plot(xmf,histm, linestyle='solid', color='red', label='$\\Delta \\rm v<1000\\rm km\\, s^{-1}$')

    ind = np.where(abs(deltav_neighbours) < 500)
    histm = np.histogram(np.log10(mstars_neighbours[ind]),bins=np.append(mbins,mupp))[0]/npass
    ax.plot(xmf,histm, linestyle='solid', color='crimson', label='$\\Delta \\rm v<500\\rm km\\, s^{-1}$')

    common.prepare_legend(ax, ['k', 'red', 'crimson'], loc=1)

    xtit = "$\\rm log_{10} (\\rm SFR/M_{\odot}\\, yr^{-1})$"
    xmin, xmax, ymin, ymax = -3, 3, 0.9, 150

    ax = fig.add_subplot(222)
    ax.set_yscale('log')
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 100, 100))

    ind = np.where(sfr_neighbours == 0)
    sfr_neighbours[ind] = 1e-3

    histm = np.histogram(np.log10(sfr_neighbours),bins=np.append(sfrbins,supp))[0]/npass
    ax.plot(xsfr,histm, linestyle='solid', color='black')
    ind = np.where(abs(deltav_neighbours) < 1000)
    histm = np.histogram(np.log10(sfr_neighbours[ind]),bins=np.append(sfrbins,supp))[0]/npass
    ax.plot(xsfr,histm, linestyle='solid', color='red')
    ind = np.where(abs(deltav_neighbours) < 500)
    histm = np.histogram(np.log10(sfr_neighbours[ind]),bins=np.append(sfrbins,supp))[0]/npass
    ax.plot(xsfr,histm, linestyle='solid', color='crimson')
    ind = np.where((abs(deltav_neighbours) < 1000) & (mstars_neighbours > 10**8.5))
    histm = np.histogram(np.log10(sfr_neighbours[ind]),bins=np.append(sfrbins,supp))[0]/npass
    ax.plot(xsfr,histm, linestyle='solid', color='ForestGreen', label = '$\\Delta \\rm v<1000\\rm km\\, s^{-1}$ & $M_{\\star}>10^{8.5}\\rm M_{\odot}$')
    common.prepare_legend(ax, ['ForestGreen'], loc=1)

    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm <N_{\\rm cum}>$"
    xmin, xmax, ymin, ymax = 7, 11, 0.9, 500
    ax = fig.add_subplot(223)
    ax.set_yscale('log')
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.5, 0.5, 100, 100))

    def cum_calc(hist):
        cum = np.zeros(shape = len(hist))
        for j in range(0,len(cum)):
            cum[j] = np.sum(hist[j:len(hist)])
        return cum

    histm = np.histogram(np.log10(mstars_neighbours),bins=np.append(mbins,mupp))[0]/npass
    ax.plot(xmf, cum_calc(histm), linestyle='solid', color='black', label='all')

    ind = np.where(abs(deltav_neighbours) < 1000)
    histm = np.histogram(np.log10(mstars_neighbours[ind]),bins=np.append(mbins,mupp))[0]/npass
    ax.plot(xmf,cum_calc(histm), linestyle='solid', color='red', label='$\\Delta \\rm v<1000\\rm km\\, s^{-1}$')

    ind = np.where(abs(deltav_neighbours) < 500)
    histm = np.histogram(np.log10(mstars_neighbours[ind]),bins=np.append(mbins,mupp))[0]/npass
    ax.plot(xmf,cum_calc(histm), linestyle='solid', color='crimson', label='$\\Delta \\rm v<500\\rm km\\, s^{-1}$')

    def plot_lines_interest(ax):
        ax.plot([xmin, xmax], [100,100], linestyle='dotted', color='gray')
        ax.text(xmax-0.1*(xmax-xmin), 100+5, '100')
        ax.plot([xmin, xmax], [50,50], linestyle='dotted', color='gray')
        ax.text(xmax-0.1*(xmax-xmin), 50+2, '50')
        ax.plot([xmin, xmax], [30,30], linestyle='dotted', color='gray')
        ax.text(xmax-0.1*(xmax-xmin), 30+1, '30')
        ax.plot([xmin, xmax], [10,10], linestyle='dotted', color='gray')
        ax.text(xmax-0.1*(xmax-xmin), 10+0.2, '10')

    plot_lines_interest(ax)

    xtit = "$\\rm log_{10} (\\rm SFR/M_{\odot}\\, yr^{-1})$"
    xmin, xmax, ymin, ymax = -3, 3, 0.9, 500

    ax = fig.add_subplot(224)
    ax.set_yscale('log')
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 100, 100))

    ind = np.where(sfr_neighbours == 0)
    sfr_neighbours[ind] = 1e-3

    histm = np.histogram(np.log10(sfr_neighbours),bins=np.append(sfrbins,supp))[0]/npass
    ax.plot(xsfr,cum_calc(histm), linestyle='solid', color='black')
    ind = np.where(abs(deltav_neighbours) < 1000)
    histm = np.histogram(np.log10(sfr_neighbours[ind]),bins=np.append(sfrbins,supp))[0]/npass
    ax.plot(xsfr,cum_calc(histm), linestyle='solid', color='red')
    ind = np.where(abs(deltav_neighbours) < 500)
    histm = np.histogram(np.log10(sfr_neighbours[ind]),bins=np.append(sfrbins,supp))[0]/npass
    ax.plot(xsfr,cum_calc(histm), linestyle='solid', color='crimson')
    ind = np.where((abs(deltav_neighbours) < 1000) & (mstars_neighbours > 10**8.5))
    histm = np.histogram(np.log10(sfr_neighbours[ind]),bins=np.append(sfrbins,supp))[0]/npass
    ax.plot(xsfr,cum_calc(histm), linestyle='solid', color='ForestGreen', label = '$\\Delta \\rm v<1000\\rm km\\, s^{-1}$ & $M_{\\star}>10^{8.5}\\rm M_{\odot}$')
    plot_lines_interest(ax)
    plt.tight_layout()
    common.savefig(outdir, fig, 'histogram_neighbours_quenchedz5galaxies.pdf')


def prepare_data(hdf5_data, hdf5_data_halo, index, redshift):


    # Unpack data
    (h0, volh, typeg, mdisk, mbulge, sfrd, sfrb, idhalo, mvir, mbh, xg, yg, zg, vxg, vyg, vzg) = hdf5_data

    (_, _, mvirz0, idhalo_cat) = hdf5_data_halo

    vol = volh/h0**3
    xg = xg / h0 / (1 + redshift)
    yg = yg / h0 / (1 + redshift)
    zg = zg / h0 / (1 + redshift)

    #look at number densities of galaxies with sSFR<1e-10yr^-1
    ms_tot = ((mdisk+mbulge)/h0)
    sfr_tot = ((sfrd + sfrb)/h0/1e9)

    passive = np.where((ms_tot>=10**10.4) & (sfr_tot/ms_tot <= 2e-10))
    npass = len(ms_tot[passive])
    xgp = xg[passive]
    ygp = yg[passive]
    zgp = zg[passive]
    vzgp = vzg[passive]
    scatter = np.random.normal(0.0, 0.25, len(ms_tot))

    dproj = 0.4 #pMpc
    mstars_neighbours = np.array([])
    sfr_neighbours = np.array([])
    deltav_neighbours = np.array([])

    for p in range(0,npass):
        dist = np.sqrt((xg - xgp[p])**2 + (yg - ygp[p])**2)
        deltav = (vzg - vzgp[p]) + abs(zg - zgp[p]) * us.hubble_constant(redshift) #define velocity separation, including Hubble expansion

        close = np.where((dist <= dproj) & (dproj > 0) & (ms_tot > 0))
        mstars_neighbours = np.append(mstars_neighbours, ms_tot[close])
        sfr_neighbours = np.append(sfr_neighbours, sfr_tot[close])
        deltav_neighbours = np.append(deltav_neighbours, deltav[close])

    return(npass, mstars_neighbours, sfr_neighbours, deltav_neighbours)

def main(model_dir, output_dir, redshift_table, subvols, obs_dir):


    plt = common.load_matplotlib()

    zlist = np.array([4.9]) #, 8.0, 9.0, 10.0, 11.0, 12.0])
    #zlist = np.array([5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0])
    fields = {'galaxies': ('type', 'mstars_disk', 'mstars_bulge', 'sfr_disk', 'sfr_burst', 'id_halo_tree', 'mvir_hosthalo', 'm_bh', 'position_x', 'position_y', 'position_z', 'velocity_x', 'velocity_y', 'velocity_z')}

    fields_halo = {'halo': ('final_z0_mvir', 'halo_id')}

    for index, snapshot in enumerate(redshift_table[zlist]):

        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        hdf5_data_halo = common.read_data(model_dir, snapshot, fields_halo, subvols)

        npass, mstars_neighbours, sfr_neighbours, deltav_neighbours = prepare_data(hdf5_data, hdf5_data_halo, index, zlist[index])

    plot_mass_sfrs_neighbours(plt, output_dir, obs_dir, npass, mstars_neighbours, sfr_neighbours, deltav_neighbours)

if __name__ == '__main__':
    main(*common.parse_args())
