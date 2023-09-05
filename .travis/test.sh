#!/bin/bash
#
# Travis CI test script
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2018
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
#
# Contributed by Rodrigo Tobar
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA 02111-1307  USA
#

fail() {
	echo -e "$@" 1>&2
	exit 1
}

cd ${TRAVIS_BUILD_DIR}/build

# Run unit tests first
make CTEST_OUTPUT_ON_FAILURE=1 test || fail "unit tests failed"

# Run shark as a whole using the test data and the sample config file
mkdir input || fail "failed to create input/ directory"
curl -L -o input/redshifts.txt 'https://docs.google.com/uc?export=download&id=1xvNmJB_KmoBHuQz-QzdPnY0HFs7smkUB' || fail "failed to download redshifts file"
curl -L -o input/tree_199.0.hdf5 'https://docs.google.com/uc?export=download&id=1JDK8ak13bEhzg9H9xt0uE8Fh_2LD3KpZ' || fail "failed to download test hdf5 file"

run_shark() {
	model_name=$1; shift
	./shark ../sample.cfg \
	    -o simulation.redshift_file=input/redshifts.txt \
	    -o simulation.tree_files_prefix=input/tree_199 \
	    -o execution.name_model=$model_name $@ || fail "failure during execution of shark"
}

run_shark my_model -o execution.seed=123456

# Generate the HDF5 output documentation and check it's up to date
# otherwise tell the user how to update it
check_hdf5_doc() {
   ../scripts/properties_as_list.sh mini-SURFS/my_model/$1 > props.rst
	_diff="`diff -Naur ../doc/hdf5_properties/$2 props.rst`"
	if [ -n "${_diff}" ]; then
		fail "\nThe file doc/hdf5_properties/$2 is out of date. This probably means that you added a new\n" \
		     "dataset to shark's output, but forgot to update the corresponding documentation.\n" \
		     "The full difference follows:\n\n${_diff}\n\n" \
		     "Please run the script/properties_as_lish.sh script against a `basename $1` file\n" \
		     "to re-generate its documentation, then commit your changes. For example:\n\n" \
		     "scripts/properties_as_list.sh my-output/model/199/0/`basename $1` > doc/hdf5_properties/$2"
	fi
}

check_hdf5_doc 199/0/galaxies.hdf5 galaxies.rst
check_hdf5_doc 156/0/star_formation_histories.hdf5 star_formation_histories.rst

if [ -n "$PYTHON" ]; then

	# How many CPUs do we have?
	if [ $TRAVIS_OS_NAME = osx ]; then
		n_cpus=`sysctl -n hw.ncpu`
	else
		n_cpus=`grep -c processor /proc/cpuinfo`
	fi

	compare_galaxies() {
		model_name=$1; shift
		"$PYTHON" ../scripts/compare_galaxies.py \
		    -m "mini-SURFS/my_model/199/0/galaxies.hdf5" \
		       "mini-SURFS/$model_name/199/0/galaxies.hdf5" $@
	}

	# Make sure the standard plotting scripts run correctly
	echo "backend: Agg" >> matplotlibrc
	"$PYTHON" ../standard_plots/all.py -c ../sample.cfg -z input/redshifts.txt || fail "failure during execution of python plotting scripts"

	# Make sure the seed value returns a reproducible result
	run_shark my_model_same_seed -o execution.seed=123456
	compare_galaxies my_model_same_seed || fail "Models expected to be equal, they are not."

	# Like above, but in parallel using all available CPUs
	run_shark my_model_same_seed_parallel -o execution.seed=123456 -t $n_cpus
	compare_galaxies my_model_same_seed_parallel || fail "Models expected to be equal, they are not."

	# Run using a random seed in the interval 2^32 - 1
	run_shark my_model_random_seed || fail "failure during execution of shark"
	compare_galaxies my_model_random_seed --expect-unequal || fail "Models expected to be unequal, they are not."
fi
