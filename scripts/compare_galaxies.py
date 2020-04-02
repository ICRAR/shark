# Test two shark galaxies for equality
#
# Copyright by Kai Striega
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
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#

import argparse

import h5py
import numpy as np


def _check(name, a, b, equality_condition):
    if not all(equality_condition):
        not_equal = np.logical_not(equality_condition)
        raise AssertionError(
            'Galaxies dataset %s not equal: %r / %r'
            % (name, a[not_equal], b[not_equal])
        )

def full_dataset_equality(name, a, b):
    _check(name, a, b, np.equal(a, b))


def lenient_dataset_equality(name, a, b):
    a, b = np.sort(a), np.sort(b)
    _check(name, a, b, np.isclose(a, b, equal_nan=True))


def read_args():
    """Return an argparse.Namespace object with the CLI arguments"""
    arg_parser = argparse.ArgumentParser(
        "Test whether two shark models are equivalent."
        )
    arg_parser.add_argument(
        '-e', '--expect-unequal', action='store_true',
        help='Whether the models are expected to be unequal.')
    arg_parser.add_argument(
        '-l', '--lenient', action='store_true',
        help='Use lenient comparison (sorted values, not exact equality).')
    arg_parser.add_argument(
        '-m', '--models', required=True, nargs=2,
        help='Path where each model is found.'
        )
    arg_parser.add_argument(
        '-E', '--exclude-dataset', action='append',
        help='Datasets to exclude from comparison')
    arg_parser.add_argument(
        '-i', '--include-dataset', action='append',
        help='Datasets to include in comparison')
    return arg_parser.parse_args()


def assert_galaxies_equal(check, galaxy1, galaxy2, inclusions, exclusions):
    """Raise an AssertionError if two galaxies are not equal."""
    if inclusions:
        names1 = names2 = inclusions
    else:
        names1 = set(galaxy1.keys()) - set(exclusions)
        names2 = set(galaxy2.keys()) - set(exclusions)
    if names1 != names2:
        raise AssertionError('Galaxy datasets unequal')
    for name in names1:
        ds1 = galaxy1[name]
        ds2 = galaxy2[name]
        if ds1.shape != ds2.shape:
            raise AssertionError('Dataset shapes unequal: %r / %r' % (ds1.shape, ds2.shape))
        check(name, ds1, ds2)


def assert_galaxies_not_equal(check, galaxy1, galaxy2, inclusions, exclusions):
    """Raise an AssertionError if two galaxies are equal."""
    try:
        assert_galaxies_equal(check, galaxy1, galaxy2, inclusions, exclusions)
    except AssertionError as e:
        pass
    else:
        raise AssertionError('Galaxies expected to be unequal, but are equal.')

def main():
    args = read_args()
    exclusions = args.exclude_dataset or []
    inclusions = args.include_dataset or []
    check = lenient_dataset_equality if args.lenient else full_dataset_equality
    model_one, model_two = h5py.File(args.models[0], 'r'), h5py.File(args.models[1], 'r')
    galaxies_one, galaxies_two = model_one['galaxies'], model_two['galaxies']

    if args.expect_unequal:
        assert_galaxies_not_equal(check, galaxies_one, galaxies_two, inclusions, exclusions)
    else:
        assert_galaxies_equal(check, galaxies_one, galaxies_two, inclusions, exclusions)

if __name__ == '__main__':
    main()
