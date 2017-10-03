# tests CMakeLists.txt
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2017
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
#
# Contributed by Rodrigo Tobar
#
# This file is part of SHArk.
#
# libprofit is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# libprofit is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with libprofit.  If not, see <http://www.gnu.org/licenses/>.


# Note: This file is included from the main CMakeLists.txt
#       so we skip most of the configuration here and go
#       straight to the point

include_directories(${CXXTEST_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
set(CXXTEST_TESTGEN_ARGS --error-printer --have-eh)

set(SHARK_TEST_NAMES options)

foreach(test_name ${SHARK_TEST_NAMES})
	CXXTEST_ADD_TEST(test_${test_name} test_${test_name}.cpp ${CMAKE_CURRENT_SOURCE_DIR}/test_${test_name}.h)
	target_link_libraries(test_${test_name} sharklib)
endforeach()
