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

import numpy as np
import collections
import logging

import common
import utilities_statistics as us


observation = collections.namedtuple('observation', 'label x y yerrup yerrdn err_absolute')

logger = logging.getLogger(__name__)

mlow = 8
mupp = 15
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

mrlow = -5
mrupp = 5
dmr = 0.2
mrbins = np.arange(mrlow,mrupp,dmr)
xmrf = mrbins + dmr/2.0

def prepare_data(hdf5_data, hdf5_data_subh, index):

    Omegab = 0.0491
    OmegaM = 0.3121

    (h0, volh, mdisk, mbulge, mhalo, typeg, id_halo, mvir_subhalo, id_subhalo,  mHId, mHIb, m_bh, x, y, z) = hdf5_data
    (h0, volh, zinfall, idsubh) = hdf5_data_subh
    mHI = mHId + mHIb

    ind = np.where(typeg == 2)
    print(x[ind], y[ind], z[ind])

    ind = np.where((typeg == 0) & (mdisk+mbulge > 0) & (m_bh <= 0) & (mvir_subhalo > 1e11))
    print(np.median(mvir_subhalo[ind]), len(mvir_subhalo[ind]), np.median(mdisk[ind]+mbulge[ind]))

    ind = np.where((mvir_subhalo > 1e11) & (mvir_subhalo < 2e11) & (typeg == 0))
    print(np.median(mHI[ind]))
    ind = np.where((mvir_subhalo > 1e11) & (mvir_subhalo < 2e11) & (typeg == 1))
    print(np.median(mHI[ind]))

    #match galaxies to their subhalo
    ind= np.where(typeg == 1)
    idsubhalos_types1 = id_subhalo[ind]
    mHI_type1 = mHI[ind]
    mvir_subhalo_type1 = mvir_subhalo[ind]
    zinfall_galaxies = np.zeros(shape = (len(idsubhalos_types1)))
    for i,g in enumerate(idsubhalos_types1):
        match = np.where(idsubh == g)
        zinfall_galaxies[i] = zinfall[match]

    print(min(zinfall_galaxies), max(zinfall_galaxies))
    ind = np.where((mvir_subhalo_type1 > 1e11) & (mvir_subhalo_type1 < 2e11) & (zinfall_galaxies < 0.1))
    print(np.median(mHI_type1[ind]))

    ind = np.where((mvir_subhalo_type1 > 1e11) & (mvir_subhalo_type1 < 2e11) & (zinfall_galaxies > 0.1) & (zinfall_galaxies < 0.3))
    print(np.median(mHI_type1[ind]))

    ind = np.where((mvir_subhalo_type1 > 1e11) & (mvir_subhalo_type1 < 2e11) & (zinfall_galaxies > 0.3) & (zinfall_galaxies < 0.6))
    print(np.median(mHI_type1[ind]))

    ind = np.where((mvir_subhalo_type1 > 1e11) & (mvir_subhalo_type1 < 2e11) & (zinfall_galaxies > 0.6) & (zinfall_galaxies < 1))
    print(np.median(mHI_type1[ind]))

    ind = np.where((mvir_subhalo_type1 > 1e11) & (mvir_subhalo_type1 < 2e11) & (zinfall_galaxies > 1) & (zinfall_galaxies < 3))
    print(np.median(mHI_type1[ind]))

    return h0

def main(modeldir, outdir, redshift_table, subvols, obsdir):

    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mvir_hosthalo',
                           'type','id_halo_tree', 
                           'mvir_subhalo', 'id_subhalo_tree', 'matom_disk', 'matom_bulge','m_bh',  'position_x', 'position_y', 'position_z')}

    fields_subh = {'subhalo': ('infall_time_subhalo', 'id')}

    zlist = [0] #, 0.5, 1, 2, 3, 4)
    snapshots = redshift_table[zlist]

    for idx, snapshot in enumerate(snapshots):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        hdf5_data_subh = common.read_data(modeldir, snapshot, fields_subh, subvols)
        h0 = prepare_data(hdf5_data, hdf5_data_subh, idx)


if __name__ == '__main__':
    main(*common.parse_args(requires_observations=True))
