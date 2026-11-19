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

import collections
import functools
import logging
import math

import numpy as np

import common
import utilities_statistics as us


observation = collections.namedtuple('observation', 'label x y yerrup yerrdn err_absolute')

logger = logging.getLogger(__name__)

##################################
# Constants
GyrToYr = 1e9
Zsun = 0.0127
XH = 0.72
MpcToKpc = 1e3

def prepare_data(hdf5_data, seds, seds_ap, index):

    (h0, volh, mdisk, mbulge, sfrd, sfrb, typeg, rstar_disk, rstar_bulge, mvir, mzd, mzb, mgd, mgb, mburst_mer, mburst_di, mass_mer, mass_di, xg, yg, zg, vx, vy, vz) = hdf5_data

    mab_dust = seds[1]
    map_dust = seds_ap[1]
    nbands = len(mab_dust[:,0])
    nbands_JWST = 17
    volnoh = volh/h0**3

   #(0): "F070W_JWST", "F090W_JWST", "F115W_JWST", "F150W_JWST", "F200W_JWST",
   #(5): "F277W_JWST", "F356W_JWST", "F444W_JWST", "F560W_JWST", "F770W_JWST",
   #(10): "F1000W_JWST", "F1130W_JWST", "F1280W_JWST", "F1500W_JWST",
   #(14): "F1800W_JWST", "F2100W_JWST", "F2550W_JWST",
   #(17): "Band_ionising_photons", "FUV_Nathan", "Band9_ALMA", "Band8_ALMA",
   #(21): "Band7_ALMA", "Band6_ALMA", "Band4_ALMA", "Band3_ALMA", "BandX_VLA",
   #(26): "BandC_VLA", "BandS_VLA", "BandL_VLA", "Band_610MHz", "Band_325MHz",
   #(31): "Band_150MHz"

    print("volume in cMpc^3:", volnoh)
    Zgas = ((mzd + mzb) / (mgd + mgb) / Zsun)

    ssfr = (sfrd + sfrb) / 1e9 / (mdisk + mbulge)
    sfr = (sfrd + sfrb) / 1e9 / h0
    ind = np.where(ssfr <= 1e-13)
    ssfr[ind] = 1e-13
    ssfr = np.log10(ssfr)

    ind   = np.where((mdisk+mbulge) > 0.0)
    mass  = np.log10(mdisk[ind] + mbulge[ind]) - np.log10(float(h0))
    sfr   = sfr[ind]
    ssfr  = ssfr[ind]
    Zgas  = Zgas[ind]
    mvir  = np.log10(mvir[ind]/h0)
    typeg = typeg[ind]
    xg = xg[ind]/h0
    yg = yg[ind]/h0
    zg = zg[ind]/h0
    vxg = vx[ind]
    vyg = vy[ind]
    vzg = vz[ind]

    ind = np.where(mass >=8)
    props = np.zeros(shape=(len(mass[ind]), 11 + nbands_JWST * 2))
    props[:,0] = mass[ind]
    props[:,1] = sfr[ind]
    props[:,2] = typeg[ind]
    props[:,3] = mvir[ind]
    props[:,4] = np.log10(Zgas[ind])
    props[:,5] = xg[ind]
    props[:,6] = yg[ind]
    props[:,7] = zg[ind]
    props[:,8] = vxg[ind]
    props[:,9] = vyg[ind]
    props[:,10] = vzg[ind]

    mab_dust_in = mab_dust[:,ind]
    map_dust_in = map_dust[:,ind]
    mab_dust_in = mab_dust_in[:,0,:]
    map_dust_in = map_dust_in[:,0,:]

    for j in range(0,nbands_JWST):
        props[:,j+11] = mab_dust_in[j,:]
        props[:,j+11+nbands_JWST] = map_dust_in[j,:]


    np.savetxt("AllGalaxiesMstarGT1e8_z3_Lagos24_mediSURFS.txt", props)

def main(model_dir, outdir, redshift_table, subvols, obsdir):

    zlist = [3]
    plt = common.load_matplotlib()

    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'sfr_disk', 'sfr_burst', 'type', 'rstar_disk', 'rstar_bulge', 'mvir_hosthalo',
                           'mgas_metals_disk', 'mgas_metals_bulge', 'mgas_disk', 'mgas_bulge', 'mstars_burst_mergers',
                           'mstars_burst_diskinstabilities', 'mstars_bulge_mergers_assembly', 'mstars_bulge_diskins_assembly',  'position_x', 
                           'position_y', 'position_z', 'velocity_x', 'velocity_y', 'velocity_z')}
    fields_sed = {'SED/ab_dust': ('disk', 'total')}
    fields_sed_ap = {'SED/ap_dust': ('disk', 'total')}

    file_hdf5_sed = "Shark-SED-JWST-eagle-rr14.hdf5"

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        seds = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed, subvols, file_hdf5_sed)
        seds_ap = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed_ap, subvols, file_hdf5_sed)
        prepare_data(hdf5_data, seds, seds_ap, index)

if __name__ == '__main__':
    main(*common.parse_args())
