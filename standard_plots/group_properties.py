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

G        = 4.299e-9 #Gravity constant in units of (km/s)^2 * Mpc/Msun



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
     mHI_bulge, mH2_bulge, mgas_bulge, mhalo, id_halo, 
     mhot, x, y, z, vx, vy, vz, vvir) = hdf5_data

    XH = 0.72
    h0log = np.log10(float(h0))

    n_typeg = len(typeg)
    morpho_type = np.zeros(shape = (n_typeg))
    morpho_type_stellar = np.zeros(shape = (n_typeg))

    mHI_halos_stacking = np.zeros(shape = (3,len(xmf)))

    (unique_elements, counts_elements) = np.unique(id_halo, return_counts=True)

    idmax = max(unique_elements)
    print('number of halos: %d', len(unique_elements))
    ind = np.where(typeg == 0)
    print('number of central galaxies: %d', len(typeg[ind]))

  
    mHI_halo   = np.zeros(shape = len(counts_elements))
    mHI_halo_all   = np.zeros(shape = len(counts_elements))
    mmass_halo = np.zeros(shape = len(counts_elements))
    rvir_halo  = np.zeros(shape = len(counts_elements))
    vvir_halo  = np.zeros(shape = len(counts_elements))
    id_unique_halo = np.zeros(shape = len(counts_elements))
    xyz_halo       = np.zeros(shape = (3,len(counts_elements)))
    v_xyz_halo     = np.zeros(shape = (3,len(counts_elements)))

    print("will create vectors of halo mass and HI mass of individual groups")
    #create vector with halo masses and total HI masses in halos.
    for i in range(0, len(counts_elements)):
	#select galaxies that belong to this halo
        ind = np.where(id_halo == unique_elements[i])
        indc = np.where((typeg == 0) & (id_halo == unique_elements[i]))
        xyz_halo[0,i] = x[indc]
        xyz_halo[1,i] = y[indc]
        xyz_halo[2,i] = z[indc]
        mHIallgals_halo = mHI[ind] + mHI_bulge[ind]
        vvir_halo[i]    = vvir[indc]
        rvir_halo[i]    = G * mhalo[indc] / pow(vvir_halo[i], 2.0)
        dist_to_cen = np.sqrt((x[ind] - xyz_halo[0,i])**2.0 + (y[ind] - xyz_halo[1,i])**2.0)
        inr = np.where(dist_to_cen/rvir_halo[i] < 1.0)

        if(len(mhalo[ind]) > 0):
                total_bar_mass = sum(mdisk[ind]) + sum(mbulge[ind]) + sum(mgas[ind]) + sum(mgas_bulge[ind]) + sum(mhot[ind])
                mhalo_all      = mhalo[ind]
                vvir_all       = vvir[ind]
	        mmass_halo[i]  = mhalo_all[0]# + total_bar_mass
	        mHI_halo[i]    = (sum(mHIallgals_halo[inr]) * XH) #only HI
                mHI_halo_all[i] = (sum(mHI[ind] + mHI_bulge[ind]) * XH)
                vvir_halo[i]   = vvir[0]
                rvir_halo[i]   = G * mmass_halo[i] / pow(vvir_halo[i], 2.0)
                id_unique_halo[i] = unique_elements[i]
        #select central galaxy of this halo to assign positions and velocities to this halo.
        ind = np.where((typeg == 0) & (id_halo == unique_elements[i]))
        if(len(mhalo[ind]) > 0):
           xyz_halo[0,i] = x[ind]
           xyz_halo[1,i] = y[ind]
           xyz_halo[2,i] = z[ind]
           v_xyz_halo[0,i] = vx[ind]
           v_xyz_halo[1,i] = vy[ind]
           v_xyz_halo[2,i] = vz[ind]
 
    print("will calculate total HI mass in groups")
    for i in range(0,len(xmf)):
	mlow_r  = xmf[i] - dm/2.0
        mhigh_r = xmf[i] + dm/2.0
	#select halos in the mass range above
        ind = np.where((np.log10(mmass_halo) >= mlow_r) & (np.log10(mmass_halo) < mhigh_r))
        if(len(mmass_halo[ind]) > 0):
                mHI_halos_stacking[0,i] = np.median(np.log10(mmass_halo[ind]))
		mHI_halos_stacking[1,i] = np.log10(np.mean(mHI_halo[ind]))
		mHI_halos_stacking[2,i] = np.log10(np.mean(mHI_halo_all[ind]))


    #select all satellite galaxies in halos with masses > 10^13.
    ind = np.where((mhalo > 1e12) & (typeg > 0))
    sats_type = typeg[ind]
    sats_x = x[ind]
    sats_y = y[ind]
    sats_z = z[ind]
    sats_vx = vx[ind]
    sats_vy = vy[ind]
    sats_vz = vz[ind]
    sats_halo_id = id_halo[ind]
    sats_vproj = np.zeros(shape = len(sats_x))
    sats_rproj = np.zeros(shape = len(sats_x))

    (unique_elements, counts_elements) = np.unique(sats_halo_id, return_counts=True)
    for i in range(0, len(counts_elements)):
        #select galaxies that belong to this halo
        ind_gal = np.where(sats_halo_id == unique_elements[i])
        ind_halo = np.where(id_unique_halo == unique_elements[i])
        if(len(sats_vx[ind_gal]) > 0):
           rthis_halo = rvir_halo[ind_halo]
           vthis_halo = vvir_halo[ind_halo]
           sats_vproj[ind_gal] = (sats_vx[ind_gal] - v_xyz_halo[0,ind_halo])/vthis_halo
           sats_rproj[ind_gal] = (np.sqrt(pow(sats_x[ind_gal]-xyz_halo[0,ind_halo],2.0) + pow(sats_z[ind_gal]-xyz_halo[2,ind_halo],2.0) + pow(sats_y[ind_gal]-xyz_halo[1,ind_halo],2.0)))/rthis_halo

    return (mHI_halos_stacking, sats_vproj, sats_rproj, sats_type)

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
    ind = np.where((mHI_halos_stacking[1,:] > 0) & (xmf > 10.3))
    xplot = xmf[ind]
    yplot = mHI_halos_stacking[1,ind]-mHI_halos_stacking[0,ind]
    for i,j,x,y in zip (mHI_halos_stacking[0,ind],yplot,mHI_halos_stacking[1,ind],mHI_halos_stacking[2,ind]):
    	print i, j, x + np.log10(0.6751), y+ np.log10(0.6751)

    ax.plot(xplot,yplot[0], color='k', linestyle='solid', label='Shark')


    common.prepare_legend(ax, ['k'])
    common.savefig(output_dir, fig, "HI_groups_stacking.pdf")

def plot_caustic_halos(plt, outdir, sats_vproj, sats_rproj, sats_type):

    fig = plt.figure(figsize=(9,9))
    xtit = "$r/r_{\\rm vir}$"
    ytit = "$v_{\\rm r}/v_{\\rm vir}$"
    xmin, xmax, ymin, ymax = 0, 1.2, -5, 5
    xleg = xmin + 0.02 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #ind = np.where(sats_type == 2)
    #xdata = sats_rproj[ind]
    #ydata = sats_vproj[ind]
    #us.density_contour(ax, xdata, ydata, 30, 30) #, **contour_kwargs)

    ind = np.where(sats_type == 2)
    xdata = sats_rproj[ind]
    ydata = sats_vproj[ind]
    ax.plot(xdata,ydata,'ro',markersize=0.7) #, **contour_kwargs)

    ind = np.where(sats_type == 1)
    xdata = sats_rproj[ind]
    ydata = sats_vproj[ind]
    ax.plot(xdata,ydata,'ko',markersize=0.7) #, **contour_kwargs)

    common.savefig(outdir, fig, 'caustic_groups.pdf')

def main(model_dir, output_dir, redshift_table, subvols, obs_dir):

    plt = common.load_matplotlib()
    fields = {'galaxies': ('type', 'mstars_disk', 'mstars_bulge',
                           'rstar_disk', 'm_bh', 'matom_disk', 'mmol_disk', 'mgas_disk',
                           'matom_bulge', 'mmol_bulge', 'mgas_bulge', 'mvir_hosthalo',
                           'id_halo_tree', 'mhot', 'position_x', 'position_y', 'position_z', 
                           'velocity_x', 'velocity_y', 'velocity_z', 'vvir_hosthalo')}
    hdf5_data = common.read_data(model_dir, redshift_table[0], fields, subvols)

    (mHI_halos_stacking, sats_vproj, sats_rproj, sats_type) = prepare_data(hdf5_data)

    plot_HI_gas_fraction_groups(plt, output_dir, obs_dir, mHI_halos_stacking)
    #plot_caustic_halos(plt, output_dir, sats_vproj, sats_rproj, sats_type)

if __name__ == '__main__':
    main(*common.parse_args())
