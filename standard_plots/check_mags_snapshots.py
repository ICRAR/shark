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
"""HMF plots"""

import numpy as np
import os

import common


##################################

# Constants
GyrToYr = 1e9
Zsun = 0.0127
XH = 0.72
PI = 3.141592654
MpcToKpc = 1e3
c_light = 299792458.0 #m/s

def prepare_data(hdf5_data, phot_data, phot_data_nod, index, nbands, snapshot):
   
    #star_formation_histories and SharkSED have the same number of galaxies in the same order, and so we can safely assume that to be the case.
    #to select the same galaxies in galaxies.hdf5 we need to ask for all of those that have a stellar mass > 0, and then assume that they are in the same order.

    (h0, volh, mdisk, mbulge, mhalo, mshalo, typeg, age,
     sfr_disk, sfr_burst, id_gal) = hdf5_data

    #components:
    #(len(my_data), 2, 2, 5, nbands)
    #0: disk instability bulge
    #1: galaxy merger bulge
    #2: total bulge
    #3: disk
    #4: total
    seds_bulge_dummy = phot_data_nod[0]
    ngals = len(seds_bulge_dummy[0,:])
    SEDs_dust = np.zeros(shape = (5,nbands,ngals))
    SEDs_dust[0,:] = phot_data[0]
    SEDs_dust[1,:] = phot_data[1]
    SEDs_dust[2,:] = phot_data[2]
    SEDs_dust[3,:] = phot_data[3]
    SEDs_dust[4,:] = phot_data[4]
    SEDs_apdust = np.zeros(shape = (5,nbands,ngals))
    SEDs_apdust[0,:] = phot_data_nod[0]
    SEDs_apdust[1,:] = phot_data_nod[1]
    SEDs_apdust[2,:] = phot_data_nod[2]
    SEDs_apdust[3,:] = phot_data_nod[3]
    SEDs_apdust[4,:] = phot_data_nod[4]

    print("median absolute-apparent for the SDSS u band total", np.median(SEDs_dust[4,2,:] - SEDs_apdust[4,2,:])," at snapshot", snapshot)
    print("median absolute-apparent for the SDSS u band disk", np.median(SEDs_dust[3,2,:] - SEDs_apdust[3,2,:])," at snapshot", snapshot)
    #print("median absolute and apparent for the SDSS u band total", np.median(SEDs_dust[4,2,:]), np.median(SEDs_apdust[4,2,:])," at snapshot", snapshot)
    #print("median absolute and apparent for the SDSS u band disk", np.median(SEDs_dust[3,2,:]), np.median(SEDs_apdust[3,2,:])," at snapshot", snapshot)

    #print("median absolute-apparent for the SDSS u band manual total calculation", np.median(SEDs_dust[4,2,:] - new_tot[2,:])," at snapshot", snapshot)

def main(model_dir, outdir, redshift_table, subvols, obsdir):

    # Loop over redshift and subvolumes
    plt = common.load_matplotlib()

    Variable_Ext = True

    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mvir_hosthalo',
                           'mvir_subhalo', 'type', 'mean_stellar_age',
                           'sfr_disk', 'sfr_burst', 'id_galaxy')}

    fields_sed = {'SED/ab_dust': ('bulge_d','bulge_m','bulge_t','disk','total'),}
    fields_sed_nod = {'SED/ap_dust': ('bulge_d','bulge_m','bulge_t','disk','total')}

    fields_halo =  {'haloTrees': ('snapshotNumber', 'redshift')}


    #tree_dir = '/scratch/pawsey0119/clagos/medi-SURFS/1536/halos/trees/dhalo/hbt_trees_199/'
    #halo_data = common.read_halo_data(tree_dir, fields_halo, subvols)
    # 
    #(snaps_halo, z_halo) = halo_data
    #print(min(snaps_halo))
    ##z = (3.0, 4.0, 6.0, 8.0, 10.0)  #(8.0, 9, 10.5, 12.5) #, 1.0, 1.5, 2.0)
    snapshots = range(165,169,1) #redshift_table[z]

    file_hdf5_sed = "Shark-SED-eagle-rr14.hdf5" 
    # Create histogram
    for index, snapshot in enumerate(snapshots):
        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
    
        if(Variable_Ext == False):
           seds = common.read_photometry_data(model_dir, snapshot, fields_sed, subvols)
           seds_nod = common.read_photometry_data(model_dir, snapshot, fields_sed_nod, subvols)
        else:
           seds = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed, subvols, file_hdf5_sed)
           seds_nod = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed_nod, subvols, file_hdf5_sed)
        nbands = len(seds[0]) 
        prepare_data(hdf5_data, seds, seds_nod, index, nbands, snapshot)


if __name__ == '__main__':
    main(*common.parse_args())
