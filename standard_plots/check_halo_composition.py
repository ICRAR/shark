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


def prepare_data(hdf5_data, index, model_dir, snapshot, subvol):


    # Unpack data
    (h0, _, typeg, idgal, idhalo, mvir) = hdf5_data

    haloids = np.unique(idhalo)
    nhalos = len(haloids)
    ind = np.where(typeg == 0)
    ncens = len(typeg[ind])
    ids_withcentral = np.unique(idhalo[ind])
    mvir_in = mvir[ind]
    halos_missing_cens = haloids[np.in1d(haloids,ids_withcentral,invert = True)]
    mvir_nocens = np.log10(mvir[np.in1d(idhalo,halos_missing_cens)])
        
    num_halos_no_central = nhalos - ncens

    if(num_halos_no_central > 0):
       ind = np.where((typeg == 0) & (np.log10(mvir) > np.median(mvir_nocens) - 0.15) & (np.log10(mvir) < np.median(mvir_nocens) + 0.15))
       ngals_in_massbins = len(typeg[ind])
       print("Number of halos with no central:", num_halos_no_central, " in snapshot", snapshot, "(% in mass bin is)", (num_halos_no_central + 0.0)/ (ngals_in_massbins + num_halos_no_central + 0.0) * 100, " median Mvir", np.median(mvir_nocens), np.median(np.log10(mvir)), max(np.log10(mvir)))
    else:
       ngals_in_massbins = 0


def main(model_dir, output_dir, redshift_table, subvols, obs_dir):

    #zlist = np.arange(2,10,0.25)
    #zlist = (0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 0, 0.25, 0.5, 1, 2, 3, 4, 6, 8, 9, 10)

    zlist_given = True
    if(zlist_given):
        zlist = [0, 0.5, 1.0, 2.0, 3.0] #0.381963715160695, 1.77053590476006]
    else:
        snap_list = range(250,269,1)

    plt = common.load_matplotlib()
    fields = {'galaxies': ('type', 'id_galaxy', 'id_halo_tree', 'mvir_hosthalo')}

    if(zlist_given): 
       for index, snapshot in enumerate(redshift_table[zlist]):
           hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
           prepare_data(hdf5_data, index, model_dir, snapshot, subvols)
    else:
        for index, snapshot in enumerate(snap_list):
           hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
           prepare_data(hdf5_data, index, model_dir, snapshot, subvols)

if __name__ == '__main__':
    main(*common.parse_args())
