#!/bin/bash
#
# Travis CI before-install script
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

cd ${TRAVIS_BUILD_DIR}

# In MacOS we simply need to brew install some things
if [ "${TRAVIS_OS_NAME}" = "osx" ]
then

	# cxxtest pulls python@2, so we need to unlink
	# the pre-installed python first
	brew unlink python || fail "cannot unlink python"

	# Minimal dependencies for testing
	pkgs="gsl hdf5 cxxtest"
	if [ "$PYTHON" = "venv" ]
	then
		pkgs="$pkgs python3"
	fi

	brew install $pkgs || fail "cannot install packages: $pkgs"

	# PYTHON==venv means that we are going to use a virtualenv'd python
	# to run the standard plot scripts; therefore we will be testing
	# against the latest version of all python packages
	if [ "$PYTHON" = "venv" ]
	then
		python3 -mvenv "${TRAVIS_BUILD_DIR}/shark-venv" || fail "cannot create virtualenv"
		source "${TRAVIS_BUILD_DIR}/shark-venv/bin/activate"
		pip install -U pip wheel setuptools h5py matplotlib scipy || fail "cannot install python packages in virtualenv"
		deactivate
		PYTHON="${TRAVIS_BUILD_DIR}/shark-venv/bin/python3"
	fi

	return
fi

# Ubuntu Travis still comes with GSL 1.X but we need >= 2
# We cache the binary version through travis' cache, so let's
# check first if it exists
export GSL_ROOT_DIR=${TRAVIS_BUILD_DIR}/gsl/2.4
export LD_LIBRARY_PATH=${GSL_ROOT_DIR}/lib:$LD_LIBRARY_PATH
if [ ! -d "${GSL_ROOT_DIR}/lib" ]
then
	curl -O https://mirror.freedif.org/GNU/gsl/gsl-2.4.tar.gz
	tar xf gsl-2.4.tar.gz
	cd gsl-2.4
	./configure --prefix=${GSL_ROOT_DIR}
	make all -j2
	make install
	cd ..
fi
