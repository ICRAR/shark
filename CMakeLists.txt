# CMakeLists.txt
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2016
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
#
# Contributed by Rodrigo Tobar
#
# This file is part of SHArk.
#
# shark is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# shark is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with libprofit.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.1)

# We have the version here
file(STRINGS ${CMAKE_SOURCE_DIR}/VERSION PROFIT_VERSION)

# We have only C++ sources, but until cmake 3.4 the check_symbol_exists
# macro didn't work properly if only defining C++ as the project language
if( ${CMAKE_MAJOR_VERSION} EQUAL 3 AND ${CMAKE_MINOR_VERSION} LESS 4 )
	set(PROFIT_LANGS C CXX)
else()
	set(PROFIT_LANGS CXX)
endif()

# Our project...
project(SHArk VERSION ${SHARK_VERSION} LANGUAGES ${PROFIT_LANGS})
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMakeModules/)

# We use HDF5 here
find_package(HDF5 REQUIRED COMPONENTS CXX)
set(SHARK_LIBS ${HDF5_CXX_LIBRARIES})
include_directories(${HDF5_INCLUDE_DIRS})

# Common definitions for the shared lib and the binary
add_definitions(-DSHARK_VERSION="${SHARK_VERSION}")
include_directories(${CMAKE_SOURCE_DIR}/include)

# The shark shared library
set(SHARKLIB_SRCS
   src/options.cpp
   src/utils.cpp
   src/hdf5/reader.cpp
)
add_library(sharklib SHARED ${SHARKLIB_SRCS})
target_link_libraries(sharklib ${SHARK_LIBS})
set_target_properties(sharklib PROPERTIES LIBRARY_OUTPUT_NAME shark)

# The shark-importer executable
set(SHARK_IMPORTER_SRCS
	src/importer/descendants.cpp
	src/importer/main.cpp
	src/importer/reader.cpp
	src/importer/velociraptor.cpp
)
add_executable(shark-importer ${SHARK_IMPORTER_SRCS})
target_link_libraries(shark-importer sharklib)

# The shark executable
set(SHARK_SRCS
	src/main.cpp
	src/evolve_halos.cpp
)
add_executable(shark ${SHARK_SRCS})
target_link_libraries(shark sharklib)

# Installing stuff
install(TARGETS sharklib shark-importer
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib)
install(FILES
        include/reader.h
        DESTINATION include/shark)
