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

cd ${TRAVIS_BUILD_DIR}

# In MacOS we simply need to brew install some things
if [ "${TRAVIS_OS_NAME}" = "osx" ]
then

	# The xcode7.3 osx image needs an update
	if [ "${XCODE}" = "7.3" ]
	then
		brew update
	fi

	# The xcode 8.1 and 9.1 images need oclint to be uninstalled
	# (see travis-ci issue #8826)
	if [ "${XCODE}" = "8.1" -o "${XCODE}" = "9.1" ]
	then
		brew cask uninstall oclint
	fi

	# cxxtest pulls python@2, so we need to unlink
	# the pre-installed python first
	brew unlink python

	# Minimal dependencies for testing
	brew install gsl hdf5 cxxtest
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
