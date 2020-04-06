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
from scipy import interpolate
import common
import os
import math 
import h5py

##################################

# Constants
GyrToYr = 1e9
Zsun = 0.0127
XH = 0.72
PI = 3.141592654
MpcToKpc = 1e3
c_light = 299792458.0 #m/s
Sigma0=1e-10 #yr^{-1}*zsun
FX0=1e-2 #erg/s/cm^2
Lsun = 3.839e-7 #in 1e40 erg/s
sigma_gas = 20.0 #km/s for CO

thresh_thin_disk = 0.01
thresh_super_edd = 1.0

#the parameter below considers the difference between 
#the peak flux and the maximum in a box-shaped emission line
boost_box_profile = 1.5

gammaSFR = 1.0
alphaFx  = 1.0 
Av = 4

zsun = 0.0189

def prepare_data(hdf5_data, index, model_dir, snapshot, subvol, obsdir):

    (h0, _, mdisk, mbulge, mburst_mergers, mburst_diskins, mstars_bulge_mergers_assembly, mstars_bulge_diskins_assembly,
     mBH, rdisk, rbulge, typeg, specific_angular_momentum_disk_star, specific_angular_momentum_bulge_star,
     specific_angular_momentum_disk_gas, specific_angular_momentum_bulge_gas, specific_angular_momentum_disk_gas_atom,
     specific_angular_momentum_disk_gas_mol, lambda_sub, mvir_s, mgas_disk, mgas_bulge, matom_disk, mmol_disk, matom_bulge,
     mmol_bulge, mbh_acc_hh, mbh_acc_sb, rgas_disk, rgas_bulge) = hdf5_data

    vdisk = specific_angular_momentum_disk_gas_mol / rgas_disk / 2.0  #in km/s
    vbulge = specific_angular_momentum_bulge_gas_mol / rgas_bulge / 2.0 #in km/s


def main(model_dir, output_dir, redshift_table, subvols, obs_dir):

    plt = common.load_matplotlib()

    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mstars_burst_mergers', 'mstars_burst_diskinstabilities',
                           'mstars_bulge_mergers_assembly', 'mstars_bulge_diskins_assembly', 'm_bh', 'rstar_disk', 'rstar_bulge', 'type',
                           'specific_angular_momentum_disk_star', 'specific_angular_momentum_bulge_star',
                           'specific_angular_momentum_disk_gas', 'specific_angular_momentum_bulge_gas',
                           'specific_angular_momentum_disk_gas_atom', 'specific_angular_momentum_disk_gas_mol',
                           'lambda_subhalo', 'mvir_subhalo', 'mgas_disk', 'mgas_bulge','matom_disk', 'mmol_disk',
                           'matom_bulge', 'mmol_bulge', 'bh_accretion_rate_hh', 'bh_accretion_rate_sb','rgas_disk', 'rgas_bulge')}
    fields_co = {'galaxies': ('LCO_disk','LCO_bulge')}

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)

        prepare_data(hdf5_data, index, model_dir, snapshot, subv, obs_dir)

if __name__ == '__main__':
    main(*common.parse_args())
