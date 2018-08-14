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

import common
import utilities_statistics as us

##################################
# Constants
mlow = 10.0
mupp = 15.0
dm = 0.3
mbins = np.arange(mlow, mupp, dm)
xmf = mbins + dm/2.0


def add_observations_to_plot(obsdir, fname, ax, marker, label, color='k', err_absolute=False):
    fname = '%s/Gas/%s' % (obsdir, fname)
    x, y, yerr_down, yerr_up = common.load_observation(obsdir, fname, (0, 1, 2, 3))
    common.errorbars(ax, x, y, yerr_down, yerr_up, color, marker, label=label, err_absolute=err_absolute)

def prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit):
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)
    xleg = xmax - 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)
    #ax.text(xleg, yleg, 'z=0')

def prepare_data(hdf5_data):

    bin_it = functools.partial(us.wmedians, xbins=xmf)

    # Unpack data
    (h0, _, typeg, mdisk, mbulge, _, _, mHI, mH2, mgas,
     mHI_bulge, mH2_bulge, mgas_bulge, mhalo, id_halo, mhot) = hdf5_data

    XH = 0.72
    h0log = np.log10(float(h0))

    n_typeg = len(typeg)
    morpho_type = np.zeros(shape = (n_typeg))
    morpho_type_stellar = np.zeros(shape = (n_typeg))

    mHI_halos_stacking = np.zeros(shape = (len(xmf))) 

    idmax = max(id_halo)
    print('number of halos: %d' % (idmax))
 
    mHI_halo   = np.zeros(shape = idmax)
    mmass_halo = np.zeros(shape = idmax)

    print("will create vectors of halo mass and HI mass of individual groups")
    #create vector with halo masses and total HI masses in halos.
    for i in range(0, idmax):
	#select galaxies that belong to this halo
        ind = np.where(id_halo == i)
        if(len(mhalo[ind]) > 0):
                total_bar_mass = sum(mdisk[ind]) + sum(mbulge[ind]) + sum(mgas[ind]) + sum(mgas_bulge[ind]) + sum(mhot[ind])
                mhalo_all      = mhalo[ind]
	        mmass_halo[i]  = mhalo_all[0] + total_bar_mass
	        mHI_halo[i]    = (sum(mHI[ind]) * XH) #only HI 

    print("will calculate total HI mass in groups")
    for i in range(0,len(xmf)):
	mlow_r  = xmf[i] - dm/2.0
        mhigh_r = xmf[i] + dm/2.0
	#select halos in the mass range above
        ind = np.where((np.log10(mmass_halo) >= mlow_r) & (np.log10(mmass_halo) < mhigh_r))
        if(len(mmass_halo[ind]) > 0):
		mHI_halos_stacking[i] = np.log10(np.mean(mHI_halo[ind]))

    return (mHI_halos_stacking)

def plot_HI_gas_fraction_groups(plt, output_dir, obs_dir, mHI_halos_stacking):

    ###################################
    #   Plots global mass densities
    fig = plt.figure(figsize=(5,4.5))

    xtit="$\\rm log_{10} (\\rm M_{\\rm halo}/M_{\odot}\,h^{-1})$"
    ytit="$\\rm log_{10}(M_{\\rm HI}/M_{\\rm halo})$"

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    prepare_ax(ax, 10, 15, -5, -0.5, xtit, ytit)

    #Predicted SMHM
    ind = np.where((mHI_halos_stacking > 0) & (xmf > 10.3))
    xplot = xmf[ind]
    yplot = mHI_halos_stacking[ind]-xmf[ind]
    #for i,j in zip (xplot,yplot):
    #	print i,j

    ax.plot(xplot,yplot, color='k', linestyle='solid', label='Shark')


    common.prepare_legend(ax, ['k'])
    common.savefig(output_dir, fig, "HI_groups_stacking.pdf")


def main(model_dir, output_dir, subvols, obs_dir, snapshot):

    plt = common.load_matplotlib()
    fields = {'galaxies': ('type', 'mstars_disk', 'mstars_bulge',
                           'rstar_disk', 'm_bh', 'matom_disk', 'mmol_disk', 'mgas_disk',
                           'matom_bulge', 'mmol_bulge', 'mgas_bulge', 'mvir_hosthalo',
                           'id_halo', 'mhot')}
    hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)

    (mHI_halos_stacking) = prepare_data(hdf5_data)

    plot_HI_gas_fraction_groups(plt, output_dir, obs_dir, mHI_halos_stacking)

if __name__ == '__main__':
    main(*common.parse_args())
