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

##################################
# Mass function initialization
mlow = 7
mupp = 12
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

def prepare_data(hdf5_data, index, hist_smf):

    (h0, volh, mdisk, mbulge) = hdf5_data

    mass          = np.zeros(shape = len(mdisk))
    ind = np.where((mdisk+mbulge) > 0.0)
    mass[ind] = np.log10(mdisk[ind] + mbulge[ind]) - np.log10(float(h0))
    H, _ = np.histogram(mass,bins=np.append(mbins,mupp))
    hist_smf[index,:] = hist_smf[index,:] + H

    if volh > 0:
        vol = volh/pow(h0,3.)  # In Mpc^3
        hist_smf[index,:]  = hist_smf[index,:]/vol/dm

    return mass

def main(modeldir, outdir, redshift_table, subvols, obsdir):

    zlist = (0.15, 0.25, 0.4, 0.625, 0.875, 1.125, 1.375, 1.625, 1.875)

    plt = common.load_matplotlib()

    # Histograms
    hist_smf       = np.zeros(shape = (len(zlist), len(mbins)))

    fields = {'galaxies': ('mstars_disk', 'mstars_bulge')}

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        mass = prepare_data(hdf5_data, index, hist_smf)
        h0 = hdf5_data[0]

    # Take logs
    def take_log(array):
        ind = np.where(array > 0.)
        array[ind] = np.log10(array[ind])

    take_log(hist_smf)
    for index, snapshot in enumerate(redshift_table[zlist]):
        print("#SMF at z=%s and snapshot=%s" % (str(zlist[index]), str(redshift_table[zlist[index]])))
        for a,b in zip(xmf, hist_smf[index,:]):
            print(a,b)

if __name__ == '__main__':
    main(*common.parse_args())
