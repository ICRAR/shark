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
file(STRINGS ${CMAKE_SOURCE_DIR}/VERSION SHARK_VERSION)

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

#
# This is how we look for GSL
# Some cmake installations don't have a FindGSL, so if that fails
# we then try via pkg-config
#
macro(find_gsl)

	# Try with FindGSL.cmake first, if found
	# Otherwise try with pkg-config
	find_package(GSL)
	if( NOT GSL_FOUND )
		find_package(PkgConfig)
		if( PKG_CONFIG_FOUND )
			pkg_check_modules(GSL gsl)
		endif()
	endif()

	if( NOT GSL_FOUND )
		message(FATAL_ERROR
				  "\nNo GSL found in the system. Please install the GSL "
				  "development files and try again\n")
	elseif( GSL_VERSION LESS "2.0" )
		message(FATAL_ERROR
				  "\nGSL version ${GSL_VERSION} found, but we need at least 2.0\n")
	endif()

	include_directories(${GSL_INCLUDE_DIRS})
	set(SHARK_LIBS ${SHARK_LIBS} ${GSL_LIBRARIES})
endmacro()

#
# This is how we look for HDF5
# We are currently interested in the CXX library,
# which is a hard requirement -- if it's not there we fail misserably
#
macro(find_hdf5)
	find_package(HDF5 REQUIRED COMPONENTS CXX)
	set(SHARK_LIBS ${SHARK_LIBS} ${HDF5_CXX_LIBRARIES})
	include_directories(${HDF5_INCLUDE_DIRS})
endmacro()

#
# This is how we look for Boost
# We are currently interested in the log library
# which is a hard requirement -- if it's not there we fail misserably
#
macro(find_boost)
	find_package(Boost 1.54 REQUIRED COMPONENTS log program_options system)
	set(SHARK_LIBS ${SHARK_LIBS} ${Boost_LOG_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_SYSTEM_LIBRARY})
	include_directories(${Boost_INCLUDE_DIRS})
endmacro()

#
# Go out there and find stuff
#
find_gsl()
find_hdf5()
find_boost()

#
# Save all compile-time options here
#
configure_file(include/config.h.in include/config.h @ONLY)

# Our own include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_BINARY_DIR}/include)

# The shark shared library
set(SHARKLIB_SRCS
   src/cosmology.cpp
   src/gas_cooling.cpp
   src/options.cpp
   src/ode_solver.cpp
   src/physical_model.cpp
   src/star_formation.cpp
   src/stellar_feedback.cpp
   src/utils.cpp
   src/hdf5/iobase.cpp
   src/hdf5/reader.cpp
)
add_library(sharklib SHARED ${SHARKLIB_SRCS})
target_link_libraries(sharklib ${SHARK_LIBS})
set_target_properties(sharklib PROPERTIES LIBRARY_OUTPUT_NAME shark)

# The shark-importer executable
set(SHARK_IMPORTER_SRCS
	src/importer/descendants.cpp
	src/importer/main.cpp
	src/importer/options.cpp
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
