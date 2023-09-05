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

import common


##################################

# Constants
GyrToYr = 1e9
Zsun = 0.0127
XH = 0.72
PI = 3.141592654
MpcToKpc = 1e3
c_light = 299792458.0 #m/s

def main(model_dir, outdir, redshift_table, subvols, obsdir):

    # Loop over redshift and subvolumes
    plt = common.load_matplotlib()
    fields = {'galaxies': ('id_galaxy', 'descendant_id_galaxy')}

    #sfh_fields = {'bulges_diskins': ('star_formation_rate_histories'),
    #              'bulges_mergers': ('star_formation_rate_histories'),
    #              'disks': ('star_formation_rate_histories')}

    z = (1.90829,1.95572)#, 2.00392)
    snapshots = redshift_table[z]

    # Create histogram
    for index, snapshot in enumerate(snapshots):
        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
	if(index == 0):
           id_gal1 = hdf5_data[2]
        if(index == 1):
           desc_id_step0 = hdf5_data[3]
        
    matches = np.in1d(desc_id_step0,id_gal1)
    positive_matches = np.where(matches == True)
    ids_matched = desc_id_step0[positive_matches]

    print('Number of galaxies with descendantes %s' % len(desc_id_step0))
    print('Number of galaxies with matches %s' % len(ids_matched))

    if(len(ids_matched) < len(desc_id_step0)):
       positive_matches = np.where(matches == False)
       print('Galaxy with problems: %s' % desc_id_step0[positive_matches])
    

if __name__ == '__main__':
    main(*common.parse_args())
