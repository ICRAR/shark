# Test two SHARK models for equality based upon their constituent galaxies
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
import h5py
from numpy.testing import assert_equal, assert_array_equal, assert_raises

def assert_galaxies_equal(model1, model2):
    galaxies1 = model1['galaxies']
    galaxies2 = model2['galaxies']

    assert_equal(galaxies1.keys(), galaxies2.keys())
    for key in galaxies1.keys():
        assert_array_equal(galaxies1[key][:], galaxies2[key][:])

def assert_galaxies_not_equal(model1, model2):
    assert_raises(AssertionError, assert_galaxies_equal, model1, model2)

def main():
    model_base = h5py.File('mini-SURFS/my_model/156/0/galaxies.hdf5', 'r')
    model_equal_seed = h5py.File('mini-SURFS/my_model_equal_seed/156/0/galaxies.hdf5', 'r')
    model_unequal_seed = h5py.File('mini-SURFS/my_model_unequal_seed/156/0/galaxies.hdf5', 'r')

    # Test the two models return exactly the same galaxy formations.
    assert_galaxies_equal(model_base, model_equal_seed)
    assert_galaxies_not_equal(model_base, model_unequal_seed)

if __name__ == '__main__':
    main()
