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

	export HOMEBREW_NO_AUTO_UPDATE=1
	export HOMEBREW_NO_INSTALL_CLEANUP=1

	# cxxtest pulls python, so we need to unlink the pre-installed python@2 first
	# (only in xcode12.2 though, it doesn't come installed in the rest)
	test "${TRAVIS_OSX_IMAGE}" != xcode12.2 || brew unlink python@2 || fail "cannot unlink python@2"

	# Minimal dependencies for testing
	pkgs="gsl cxxtest libomp"
	brew install $pkgs || fail "cannot install packages: $pkgs"

	# PYTHON==venv means that we are going to use a virtualenv'd python
	# to run the standard plot scripts; therefore we will be testing
	# against the latest version of all python packages
	if [ "$PYTHON" = "venv" ]
	then
		python3 -mvenv "${TRAVIS_BUILD_DIR}/shark-venv" || fail "cannot create virtualenv"
		source "${TRAVIS_BUILD_DIR}/shark-venv/bin/activate"
		pip install -U pip wheel setuptools h5py matplotlib!=3.4 scipy || fail "cannot install python packages in virtualenv"
		deactivate
		PYTHON="${TRAVIS_BUILD_DIR}/shark-venv/bin/python3"
	fi

	return
fi
