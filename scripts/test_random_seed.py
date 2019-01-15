#!usr/bin/env python
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
