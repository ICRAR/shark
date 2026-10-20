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
"""SMF plots"""

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


def prepare_data(hdf5_data, index, zlist):

    (h0, volh, mdisk, mbulge, sfrd, sfrb, typeg, rstar_disk, rstar_bulge, mvir, mzd, mzb, mgd, mgb, mburst_mer, mburst_di, mass_mer, mass_di, mmd, mmb) = hdf5_data

    volnoh = volh/h0**3

    ssfr = (sfrd + sfrb) / 1e9 / (mdisk + mbulge)
    sfr = (sfrd + sfrb) / 1e9 / h0
    ind = np.where(ssfr <= 1e-13)
    ssfr[ind] = 1e-13
    ssfr = np.log10(ssfr)

    mass          = np.zeros(shape = len(mdisk))
    ind = np.where((mdisk+mbulge) > 0.0)
    mass[ind] = np.log10(mdisk[ind] + mbulge[ind]) - np.log10(float(h0))

    ssfr_thresh = -11 + 0.5 * zlist[index]
    if(zlist[index] > 2):
        ssfr_thresh = -10
    ind = np.where((mass >=7) & (ssfr > ssfr_thresh))
    props = np.zeros(shape=(len(mass[ind]), 4))
    props[:,0] = mass[ind]
    props[:,1] = sfr[ind]
    props[:,2] = np.log10((mgd[ind] + mgb[ind])/h0)
    props[:,3] = np.log10((mmd[ind] + mmb[ind])/h0)

    print("Maximum stellar mass,", max(mass[ind]), " at redshift", zlist[index])
    #print("Median bulge-to-total ratio:", np.median(props[:,3]))
    np.savetxt("EcoGals_snapshot_catalogues_z" +  str(zlist[index]) + "_Lagos18_mediSURFS.txt", props)


def main(modeldir, outdir, redshift_table, subvols, obsdir):

    zlist = (0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10, 11, 12, 13, 14)

    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'sfr_disk', 'sfr_burst', 'type', 'rstar_disk', 'rstar_bulge', 'mvir_hosthalo',
                           'mgas_metals_disk', 'mgas_metals_bulge', 'mgas_disk', 'mgas_bulge', 'mstars_burst_mergers',
                           'mstars_burst_diskinstabilities', 'mstars_bulge_mergers_assembly', 'mstars_bulge_diskins_assembly', 'mmol_disk',
                           'mmol_bulge')}

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        prepare_data(hdf5_data, index, zlist)

if __name__ == '__main__':
    main(*common.parse_args())
