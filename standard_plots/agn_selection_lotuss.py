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
import os
import h5py

import common
import utilities_statistics as us


def prepare_data(hdf5_data, galidx, bhidx, lmechh, lbolh, delta_t, LBT, index, model_dir, zlist):


    # Unpack data
    (h0, volh, mdisk, mbulge, sfrd, sfrb, idgal, mbh, macc_hh, macc_sb,spin, 
     qjet, lbol, typeg, mvir_hh, mvir_sh, vvir_hh, cnfw_sh, lambda_sh, mhot) = hdf5_data

    lmechh = lmechh[0]
    lbolh = lbolh[0]
    ind = np.where(((mdisk + mbulge)/h0 > 1e10) & (typeg ==0))
    gal_props = np.zeros(shape = (len(mdisk[ind]), 13))
    
    gal_props[:,0] = np.log10((mdisk[ind] + mbulge[ind])/h0)
    gal_props[:,1] = (sfrd[ind] + sfrb[ind])/h0/1e9
    gal_props[:,2] = np.log10(mbh[ind]/h0)
    gal_props[:,3] = (macc_hh[ind] + macc_sb[ind])/h0/1e9
    gal_props[:,4] = spin[ind]
    gal_props[:,5] = qjet[ind]
    gal_props[:,6] = lbol[ind]
    gal_props[:,7] = np.log10(mvir_hh[ind]/h0)
    gal_props[:,8] = np.log10(mhot[ind]/h0)
    gal_props[:,9] = vvir_hh[ind]
    gal_props[:,10] = cnfw_sh[ind]
    gal_props[:,11] = lambda_sh[ind]
    gal_props[:,12] = galidx[ind]

    np.savetxt("LoTSS_snapshot_catalogues_z" +  str(zlist[index]) + "_Lagos24_HBTTrees.txt", gal_props)

    lmech_all = np.zeros(shape = (len(mdisk[ind]), len(LBT)))
    lbol_all = np.zeros(shape = (len(mdisk[ind]), len(LBT)))
    idsin = galidx[ind]
    for j in range(0,len(mdisk[ind])):
        match = np.where(bhidx == idsin[j])
        if(len(bhidx[match]) == 1):
            lmech_all[j,:] = lmechh[match,:]
            lbol_all[j,:] = lbolh[match,:]
    np.savetxt("LoTSS_snapshot_mechanical_luminosity_history_z" +  str(zlist[index]) + "_Lagos24_HBTTrees.txt", lmech_all)
    np.savetxt("LoTSS_snapshot_bolometric_luminosity_history_z" +  str(zlist[index]) + "_Lagos24_HBTTrees.txt", lbol_all)

def main(model_dir, output_dir, redshift_table, subvols, obs_dir):


    zlist = [0, 1.0]

    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge','sfr_disk','sfr_burst', 'id_galaxy',
                           'm_bh', 'bh_accretion_rate_hh', 'bh_accretion_rate_sb', 'bh_spin', 'mechanical_power_agn', 
                           'bolometric_luminosity_agn', 'type', 'mvir_hosthalo', 'mvir_subhalo', 'vvir_hosthalo', 
                           'cnfw_subhalo', 'lambda_subhalo', 'mhot')}
    bhh_idx = {'galaxies': ('id_galaxy')}
    bhh_lmech = {'galaxies': ('mechanical_power_agn')}
    bhh_lbol = {'galaxies': ('bolometric_luminosity_agn')}
    gal_idx = {'galaxies': ('id_galaxy', 'id_subhalo_tree')}

    for index, snapshot in enumerate(redshift_table[zlist]):
        for j, ivol in enumerate(subvols):
            bh_id, delta_t, LBT = common.read_bhh(model_dir, snapshot, bhh_idx, [ivol])
            galid, subid = common.read_data(model_dir, snapshot, gal_idx, [ivol], include_h0_volh = False)
            subid = int((ivol+1) * 1e8)
            bh_id = bh_id[0] + subid
            galid = galid + subid
            if(j == 0):
               bh_ids_all = bh_id
               gal_ids_all = galid
            else:
               gal_ids_all = np.concatenate([gal_ids_all, galid])
               bh_ids_all = np.concatenate([bh_ids_all, bh_id]) 

        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        bhid, delta_t, LBT = common.read_bhh(model_dir, snapshot, bhh_idx, subvols)
        lmechh, delta_t, LBT = common.read_bhh(model_dir, snapshot, bhh_lmech, subvols)
        lbolh, delta_t, LBT = common.read_bhh(model_dir, snapshot, bhh_lbol, subvols)

        prepare_data(hdf5_data, gal_ids_all, bh_ids_all, lmechh, lbolh, delta_t, LBT, index, model_dir, zlist)
           
if __name__ == '__main__':
    main(*common.parse_args())
