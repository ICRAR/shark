#!/bin/bash
#
# Travis CI install script
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
mkdir build
cd build

SHARK_CMAKE_OPTIONS="-DCMAKE_CXX_COMPILER=$COMPILER -DSHARK_TEST=ON -DCMAKE_CXX_FLAGS=\"-Wall -Werror\""
if [ ${TRAVIS_OS_NAME} = osx ]; then
	SHARK_CMAKE_OPTIONS+=" -DOpenMP_CXX_LIB_NAMES=libomp"
	SHARK_CMAKE_OPTIONS+=" -DOpenMP_CXX_FLAGS='-Xpreprocessor -fopenmp -I/usr/local/include'"
	SHARK_CMAKE_OPTIONS+=" -DOpenMP_libomp_LIBRARY=/usr/local/lib/libomp.dylib"
	SHARK_CMAKE_OPTIONS+=" -DCMAKE_CXX_FLAGS=-Wno-unused-command-line-argument"
fi

# Go, go, go!
eval cmake .. ${SHARK_CMAKE_OPTIONS} || fail "cmake failed"
make all -j2 || fail "make failed"
cd ..
