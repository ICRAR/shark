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
	echo $1 1>&2
	exit 1
}

cd ${TRAVIS_BUILD_DIR}/build

# Run unit tests first
make CTEST_OUTPUT_ON_FAILURE=1 test || fail "unit tests failed"

# Run shark as a whole using the test data and the sample config file
mkdir input || fail "failed to create input/ directory"
curl -L -o input/redshifts.txt 'https://docs.google.com/uc?export=download&id=1xvNmJB_KmoBHuQz-QzdPnY0HFs7smkUB' || fail "failed to download redshifts file"
curl -L -o input/tree_199.0.hdf5 'https://docs.google.com/uc?export=download&id=1JDK8ak13bEhzg9H9xt0uE8Fh_2LD3KpZ' || fail "failed to download test hdf5 file"

./shark ../sample.cfg \
    -o simulation.redshift_file=input/redshifts.txt \
    -o simulation.tree_files_prefix=input/tree_199 \
    -o execution.seed=123456 \
    -o execution.name_model=my_model || fail "failure during execution of shark"

# Make sure the standard plotting scripts run correctly
if [ -n "$PYTHON" ]; then
	echo "backend: Agg" >> matplotlibrc
	"$PYTHON" ../standard_plots/all.py -c ../sample.cfg -z input/redshifts.txt || fail "failure during execution of python plotting scripts"

	# Make sure the seed value returns a reproducible result
	./shark ../sample.cfg \
	    -o simulation.redshift_file=input/redshifts.txt \
	    -o simulation.tree_files_prefix=input/tree_199 \
	    -o execution.seed=123456 \
	    -o execution.name_model=my_model_equal_seed \
	    || fail "failure during execution of shark"

	 ./shark ../sample.cfg \
	    -o simulation.redshift_file=input/redshifts.txt \
	    -o simulation.tree_files_prefix=input/tree_199 \
	    -o execution.name_model=my_model_unequal_seed \
	    || fail "failure during execution of shark"

	"$PYTHON" ../scripts/test_random_seed.py || fail "seed value was not reproducible"
fi
