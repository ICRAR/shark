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

def prepare_data(hdf5_data, hdf5_data_halo, index, pairs, zlist, ids_pairs, xg_pairs, yg_pairs, zg_pairs, vxg_pairs, vyg_pairs, vzg_pairs, redshift_pairs, mvirz0_pairs, ms_pairs, sfr_pairs):


    # Unpack data
    (h0, volh, typeg, mdisk, mbulge, sfrd, sfrb, idhalo, mvir, xg, yg, zg, vxg, vyg, vzg) = hdf5_data

    (_, _, mvirz0, idhalo_cat) = hdf5_data_halo

    vol = volh/h0**3

    #look at number densities of galaxies with sSFR<1e-10yr^-1
    ms_tot = ((mdisk+mbulge)/h0)
    sfr_tot = ((sfrd + sfrb)/h0/1e9)

    rp_CV = 4.5 * 1e-3 #in pMpc
    dv_CV = 670.0 #km/s

    rp_search = 1.1 * rp_CV
    dv_search = 1.1 * dv_CV

    err_m = np.random.normal(0.0, 0.25, len(ms_tot)) 
    err_sfr = np.random.normal(0.0, 0.25, len(sfr_tot))
    ind = np.where((ms_tot* 10**err_m >= 10**(10.65)) & ((sfr_tot * 10**err_sfr)/(ms_tot* 10**err_m) <= 1.5e-10))
    ngals = len(xg[ind])
    xgp = xg[ind] / h0 / (1 + zlist[index])
    ygp = yg[ind] / h0 / (1 + zlist[index])
    zgp = zg[ind] / h0 / (1 + zlist[index])
    vxgp = vxg[ind] 
    vygp = vyg[ind] 
    vzgp = vzg[ind] 
    msp = np.log10(ms_tot[ind]* 10**err_m[ind])
    sfrp = (sfr_tot[ind] * 10**err_sfr[ind])
    ids_gal = zlist[index] * 1e3 + np.arange(1,ngals+1,1)
    negsfr = np.where(sfrp <0)
    sfrp[negsfr] = 0
    ids_interest = idhalo[ind]
    mvir_final = np.zeros(shape = len(ms_tot[ind]))
    for g in range(0,len(ms_tot[ind])):
        idh = np.where(idhalo_cat == ids_interest[g])
        mvir_final[g] = mvirz0[idh] 

    start_point = np.where(ids_pairs != 0)
    p = len(ids_pairs[start_point])
    #loop through all passive galaxies and find pairs in projection that satisfy criterion
    for g in range(0,ngals):
        dproj1 = np.sqrt((xgp[g] - xgp)**2 + (ygp[g] - ygp)**2)
        dvproj1 = abs(vzgp[g] - vzgp)
        dproj2 = np.sqrt((xgp[g] - xgp)**2 + (zgp[g] - zgp)**2)
        dvproj2 = abs(vygp[g] - vygp)
        dproj3 = np.sqrt((zgp[g] - zgp)**2 + (ygp[g] - ygp)**2)
        dvproj3 = abs(vxgp[g] - vxgp)

        #find pairs in all projections
        p1 = np.where((dproj1 <= rp_search) & (dvproj1 <= dv_CV))
        if(len(dproj1[p1]) > 1):
           ids_in = ids_gal[p1]
           repeated_pairs = np.in1d(ids_in, ids_pairs)
           new_pairs = np.where(repeated_pairs == False)
           pairs[index] = pairs[index] + len(repeated_pairs[new_pairs])
           ids_pairs[p:p+len(repeated_pairs[new_pairs])] = ids_in[new_pairs]
           p1 = p1[0]
           xg_pairs[p:p+len(repeated_pairs[new_pairs])]  = xgp[p1[new_pairs]]
           yg_pairs[p:p+len(repeated_pairs[new_pairs])]  = ygp[p1[new_pairs]]
           zg_pairs[p:p+len(repeated_pairs[new_pairs])]  = zgp[p1[new_pairs]]
           vxg_pairs[p:p+len(repeated_pairs[new_pairs])] = vxgp[p1[new_pairs]]
           vyg_pairs[p:p+len(repeated_pairs[new_pairs])] = vygp[p1[new_pairs]]
           vzg_pairs[p:p+len(repeated_pairs[new_pairs])] = vzgp[p1[new_pairs]]
           ms_pairs[p:p+len(repeated_pairs[new_pairs])] = msp[p1[new_pairs]]
           sfr_pairs[p:p+len(repeated_pairs[new_pairs])] = sfrp[p1[new_pairs]]
           mvirz0_pairs[p:p+len(repeated_pairs[new_pairs])] = mvir_final[p1[new_pairs]]
           redshift_pairs[p:p+len(repeated_pairs[new_pairs])] = zlist[index]
           p = p + len(repeated_pairs[new_pairs])

        #proj 2
        p2 = np.where((dproj2 <= rp_search) & (dvproj2 <= dv_CV))
        if(len(dproj2[p2]) > 1):
           ids_in = ids_gal[p2]
           repeated_pairs = np.in1d(ids_in, ids_pairs)
           new_pairs = np.where(repeated_pairs == False)
           pairs[index] = pairs[index] + len(repeated_pairs[new_pairs])
           ids_pairs[p:p+len(repeated_pairs[new_pairs])] = ids_in[new_pairs]
           p2 = p2[0]
           xg_pairs[p:p+len(repeated_pairs[new_pairs])]  = xgp[p2[new_pairs]]
           yg_pairs[p:p+len(repeated_pairs[new_pairs])]  = ygp[p2[new_pairs]]
           zg_pairs[p:p+len(repeated_pairs[new_pairs])]  = zgp[p2[new_pairs]]
           vxg_pairs[p:p+len(repeated_pairs[new_pairs])] = vxgp[p2[new_pairs]]
           vyg_pairs[p:p+len(repeated_pairs[new_pairs])] = vygp[p2[new_pairs]]
           vzg_pairs[p:p+len(repeated_pairs[new_pairs])] = vzgp[p2[new_pairs]]
           ms_pairs[p:p+len(repeated_pairs[new_pairs])] = msp[p2[new_pairs]]
           sfr_pairs[p:p+len(repeated_pairs[new_pairs])] = sfrp[p2[new_pairs]]
           mvirz0_pairs[p:p+len(repeated_pairs[new_pairs])] = mvir_final[p2[new_pairs]]
           redshift_pairs[p:p+len(repeated_pairs[new_pairs])] = zlist[index]
           p = p + len(repeated_pairs[new_pairs])

        #proj 3
        p3 = np.where((dproj3 <= rp_search) & (dvproj3 <= dv_CV))
        if(len(dproj3[p3]) > 1):
           ids_in = ids_gal[p3]
           repeated_pairs = np.in1d(ids_in, ids_pairs)
           new_pairs = np.where(repeated_pairs == False)
           pairs[index] = pairs[index] + len(repeated_pairs[new_pairs])
           ids_pairs[p:p+len(repeated_pairs[new_pairs])] = ids_in[new_pairs]
           p3 = p3[0]
           xg_pairs[p:p+len(repeated_pairs[new_pairs])]  = xgp[p3[new_pairs]]
           yg_pairs[p:p+len(repeated_pairs[new_pairs])]  = ygp[p3[new_pairs]]
           zg_pairs[p:p+len(repeated_pairs[new_pairs])]  = zgp[p3[new_pairs]]
           vxg_pairs[p:p+len(repeated_pairs[new_pairs])] = vxgp[p3[new_pairs]]
           vyg_pairs[p:p+len(repeated_pairs[new_pairs])] = vygp[p3[new_pairs]]
           vzg_pairs[p:p+len(repeated_pairs[new_pairs])] = vzgp[p3[new_pairs]]
           ms_pairs[p:p+len(repeated_pairs[new_pairs])] = msp[p3[new_pairs]]
           sfr_pairs[p:p+len(repeated_pairs[new_pairs])] = sfrp[p3[new_pairs]]
           mvirz0_pairs[p:p+len(repeated_pairs[new_pairs])] = mvir_final[p3[new_pairs]]
           redshift_pairs[p:p+len(repeated_pairs[new_pairs])] = zlist[index]
           p = p + len(repeated_pairs[new_pairs])

    print("Number of passive galaxy pairs", pairs[index], " at redshift",  zlist[index], " and number density", pairs[index] / vol)
    if(zlist[index] == max(zlist)):
        pos = np.where(ids_pairs > 0)
        for a,b,c,d,e,f,g,h,i,j,k in zip(ids_pairs[pos], ms_pairs[pos], sfr_pairs[pos], xg_pairs[pos], yg_pairs[pos], zg_pairs[pos], vxg_pairs[pos], vyg_pairs[pos], vzg_pairs[pos], mvirz0_pairs[pos], redshift_pairs[pos]):
            print(a,b,c,d,e,f,g,h,i,j,k)

def main(model_dir, output_dir, redshift_table, subvols, obs_dir):

    ids_pairs = np.zeros(shape = 500)
    xg_pairs = np.zeros(shape = 500)
    yg_pairs = np.zeros(shape = 500)
    zg_pairs = np.zeros(shape = 500)
    vxg_pairs = np.zeros(shape = 500)
    vyg_pairs = np.zeros(shape = 500)
    vzg_pairs = np.zeros(shape = 500)
    redshift_pairs = np.zeros(shape = 500)
    mvirz0_pairs = np.zeros(shape = 500)
    ms_pairs = np.zeros(shape = 500) 
    sfr_pairs = np.zeros(shape = 500) 

    plt = common.load_matplotlib()

    #zlist = np.array([2, 2.1, 2.25, 2.5, 2.75, 3.0, 3.25, 3.53362989, 3.75, 4.0]) #, 8.0, 9.0, 10.0, 11.0, 12.0])
    zlist = np.array([3.53362989, 3.75, 4.0])
    #zlist = np.array([5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0])
    fields = {'galaxies': ('type', 'mstars_disk', 'mstars_bulge', 'sfr_disk', 'sfr_burst', 'id_halo_tree', 'mvir_hosthalo', 'position_x', 'position_y', 'position_z', 'velocity_x', 'velocity_y', 'velocity_z')}
    fields_halo = {'halo': ('final_z0_mvir', 'halo_id')}

    pairs = np.zeros(shape = (len(zlist)))

    for index, snapshot in enumerate(redshift_table[zlist]):

        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        hdf5_data_halo = common.read_data(model_dir, snapshot, fields_halo, subvols)
        limit = prepare_data(hdf5_data, hdf5_data_halo, index, pairs, zlist, ids_pairs, xg_pairs, yg_pairs, zg_pairs, vxg_pairs, vyg_pairs, vzg_pairs, redshift_pairs, mvirz0_pairs, ms_pairs, sfr_pairs)

    #plot_num_density_passive(plt, output_dir, obs_dir, zlist, num_densities, ssfr_thresh, mass_threshs, mass_threshs2, num_densities2, limit)
    #plot_mvir_final(plt, output_dir, obs_dir, zlist, halo_hists)

if __name__ == '__main__':
    main(*common.parse_args())
