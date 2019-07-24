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


def full_dataset_equality(name, a, b):
    if not np.array_equal(a, b):
        raise AssertionError('Galaxies dataset %s not equal' % name)


def lenient_dataset_equality(name, a, b):
    a, b = np.sort(a), np.sort(b)
    isclose = np.isclose(a, b, equal_nan=True)
    if not all(isclose):
        not_close = np.logical_not(isclose)
        raise AssertionError(
            "Galaxies dataset %s not equal: %r / %r"
            % (name, a[not_close], b[not_close])
        )


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
    return arg_parser.parse_args()


def assert_galaxies_equal(check, galaxy1, galaxy2, exclusions):
    """Raise an AssertionError if two galaxies are not equal."""
    names1 = set(galaxy1.keys()) - set(exclusions)
    names2 = set(galaxy2.keys()) - set(exclusions)
    if names1 != names2:
        raise AssertionError('Galaxy datasets unequal')
    for name in names1:
        check(name, galaxy1[name], galaxy2[name])


def assert_galaxies_not_equal(check, galaxy1, galaxy2, exclusions):
    """Raise an AssertionError if two galaxies are equal."""
    try:
        assert_galaxies_equal(check, galaxy1, galaxy2, exclusions)
    except AssertionError as e:
        pass
    else:
        raise AssertionError('Galaxies expected to be unequal, but are equal.')

def main():
    args = read_args()
    exclusions = args.exclude_dataset or []
    check = lenient_dataset_equality if args.lenient else full_dataset_equality
    model_one, model_two = h5py.File(args.models[0]), h5py.File(args.models[1])
    galaxies_one, galaxies_two = model_one['galaxies'], model_two['galaxies']

    if args.expect_unequal:
        assert_galaxies_not_equal(check, galaxies_one, galaxies_two, exclusions)
    else:
        assert_galaxies_equal(check, galaxies_one, galaxies_two, exclusions)

if __name__ == '__main__':
    main()
