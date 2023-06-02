# CMakeLists.txt
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2017
# Copyright by UWA (in the framework of the ICRAR)
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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

# cmake 3.2 comes with FindGSL.cmake
cmake_minimum_required(VERSION 3.2)

# We have the version here
file(STRINGS ${CMAKE_SOURCE_DIR}/VERSION SHARK_VERSION)

# We have only C++ sources, but we need to allow compilation of C sources
# so that finding the HDF5 API works correctly
set(SHARK_LANGS C CXX)

# Our project...
project(SHArk VERSION ${SHARK_VERSION} LANGUAGES ${SHARK_LANGS})
set(CMAKE_CXX_STANDARD 11)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# If not specified, always build in Release mode
if ("${CMAKE_BUILD_TYPE}" STREQUAL "")
	set(CMAKE_BUILD_TYPE "Release")
endif()

# Options users can give on the command line
option(SHARK_TEST       "Include test compilation in the build" OFF)
option(SHARK_NO_OPENMP  "Don't attempt to include OpenMP support in shark" OFF)

#
# Make sure we have thread support
# We don't directly use this (yet), but if we don't manually
# add this library (most likely pthreads) to the set of libraries
# that shark needs to link against we may end up in compilation
# errors due to our dependencies using pthreads and not having it
# listed at compile time.
# This whole situation looks quite buggy to me, but I'd rather add
# the pthread dependency here than spend a life trying to find out
# how to "properly" solve the problem
#
macro(find_threads)
	find_package(Threads)
	set(SHARK_LIBS ${SHARK_LIBS} ${CMAKE_THREAD_LIBS_INIT})
endmacro()

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
	elseif (GSL_VERSION VERSION_LESS "2.0")
		message(FATAL_ERROR
				  "\nGSL version ${GSL_VERSION} found, but shark needs at least 2.0\n")
	endif()

	include_directories(${GSL_INCLUDE_DIRS})
	set(SHARK_LIBS ${SHARK_LIBS} ${GSL_LIBRARIES})
endmacro()

#
# This is how we look for HDF5
# We are currently interested in the C library,
# which is a hard requirement -- if it's not there we fail miserably
#
macro(find_hdf5)

	# The CONFIG mode is needed at least when using vcpkg's hdf5
	if (SHARK_HDF5_USE_CONFIG)
		find_package(HDF5 REQUIRED COMPONENTS C CONFIG)
	else()
		find_package(HDF5 REQUIRED COMPONENTS C)
	endif()

	#
	# HDF5_VERSION is defined only in cmake 3.3+
	# Even in these newer cmakes, there are cases in which the version is still
	# not properly calculated though, so let us check simply for an empty string
	# and compile a small program that prints the version of the HDF5 found
	# by cmake.
	#
	if (NOT HDF5_VERSION)
		message("-- HDF5 version not found automatically, taking it into our hands now...")
		set(scratch_directory ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/hdf5)
		set(test_file ${scratch_directory}/cmake_version.cxx)
		file(WRITE ${test_file}
		     "#include <stdio.h>\n"
		     "#include <hdf5.h>\n"
		     "int main(int argc, char **argv) {\n"
		     "  printf(\"INFO:\" H5_VERSION);\n"
		     "  return 0;\n"
		     "}")
		try_compile(_success ${scratch_directory} ${test_file}
		            CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${HDF5_INCLUDE_DIRS}"
		            COPY_FILE ${scratch_directory}/hdf5_version)
		if (NOT _success)
			message(FATAL_ERROR "\nFailed to compile test program to find HDF5 version")
		endif()
		execute_process(COMMAND ${scratch_directory}/hdf5_version OUTPUT_VARIABLE _VER)
		string(REGEX MATCH "^INFO:([0-9]+\\.[0-9]+\\.[0-9]+)(-patch([0-9]+))?" INFO_VER "${_VER}")
		set(HDF5_VERSION ${CMAKE_MATCH_1})
		if(CMAKE_MATCH_3)
			set(HDF5_VERSION ${HDF5_VERSION}.${CMAKE_MATCH_3})
		endif()
		mark_as_advanced(HDF5_VERSION)
		message("-- Found HDF5 version ${HDF5_VERSION}")
	endif()

	# We support 1.8+
	if (HDF5_VERSION VERSION_LESS "1.8")
		message(FATAL_ERROR
				  "\nHDF5 version ${HDF5_VERSION} found, "
				  "but shark needs versions at least 1.8.0\n")
	endif()

	# We need to pass down the information about the HDF5 version to shark
	message("-- Using HDF5 version ${HDF5_VERSION}")
	message("-- HDF5 libraries / C++ libraries:  ${HDF5_LIBRARIES} / ${HDF5_CXX_LIBRARIES}")
	string(REPLACE "." ";" VERSION_LIST ${HDF5_VERSION})
	list(GET VERSION_LIST 0 HDF5_VERSION_MAJOR)
	list(GET VERSION_LIST 1 HDF5_VERSION_MINOR)
	list(GET VERSION_LIST 2 HDF5_VERSION_PATCH)
	include_directories(${HDF5_INCLUDE_DIRS})
	add_definitions(${HDF5_DEFINITIONS})
	add_definitions("-DHDF5_VERSION_MAJOR=${HDF5_VERSION_MAJOR}"
	                "-DHDF5_VERSION_MINOR=${HDF5_VERSION_MINOR}"
	                "-DHDF5_VERSION_PATCH=${HDF5_VERSION_PATCH}")

	if (SHARK_HDF5_USE_CONFIG)
		set(SHARK_LIBS ${SHARK_LIBS} hdf5::hdf5-shared)
	else()
		set(SHARK_LIBS ${SHARK_LIBS} ${HDF5_C_LIBRARIES})
	endif()
endmacro()

#
# This is how we look for Boost
# We are currently interested in the log library
# which is a hard requirement -- if it's not there we fail misserably
#
macro(find_boost)
	# Boost 1.71's own CMake config file doesn't respect the
	# Boost_USE_MULTITHREADED setting, which is important in environments where
	# both single- and multi-threaded libraries are installed (e.g., MacOS Mojave
	# with homebrew). Let's thus not use Boost's cmake config file and stick to
	# cmake's own version of it for the time being
	set(Boost_NO_BOOST_CMAKE ON)
	find_package(Boost 1.54 COMPONENTS filesystem log program_options system REQUIRED)
	set(SHARK_LIBS ${SHARK_LIBS} ${Boost_LIBRARIES})
	include_directories(${Boost_INCLUDE_DIRS})
endmacro()

#
# Find OpenMP, if required
#
macro(find_openmp)
	find_package(OpenMP)
	if (OPENMP_FOUND)

		# We require at least OpenMP 2.0, so let's double check
		# we have that at least
		# cmake support for OpenMP has varied across the 3.x series:
		# * OPENMP_FOUND and OpenMP_<lang>_FLAGS have always been there
		# * 3.7 brings OpenMP_<lang>_SPEC_DATE variables, which *may* be set
		# * 3.9 brings OpenMP_<lang>_VERSION variables, which *may* be set,
		#   and OpenMP_<lang>_LIBRARIES/LIBRARY/LIB_NAMES variables
		# This is important for Intel compilers in particular, which require
		# additional linking
		if ((DECLARED OPENMP_VERSION) AND (NOT OPENMP_VERSION VERSION_LESS 2.0))
			set(SHARK_OPENMP ON)
		elseif((DECLARED OPENMP_CXX_SPEC_DATE) AND (NOT OPENMP_CXX_SPEC_DATE LESS 200203))
			set(SHARK_OPENMP ON)
		else()
			set(OPENMP_VERSION_CHECK_SOURCE "
#include <stdio.h>
int main(int argc, char *argv[]) {
#if _OPENMP >= 200203
	return 0;
#else
	fail to compile please
#endif
}")
			set(WORK_DIR ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/OpenMPVersionCheck)
			set(SRC_FILE ${WORK_DIR}/ompver.cpp)
			file(WRITE ${SRC_FILE} "${OPENMP_VERSION_CHECK_SOURCE}")
			try_compile(COMPILE_RESULT ${CMAKE_BINARY_DIR} ${SRC_FILE}
			            COMPILE_DEFINITIONS ${OpenMP_CXX_FLAGS}
			            LINK_LIBRARIES ${OpenMP_CXX_LIBRARIES})

			if (COMPILE_RESULT)
				set(SHARK_OPENMP ON)
			endif()
		endif()

		if (NOT SHARK_OPENMP)
			message(STATUS "OpenMP found, but is <= 2.0. Compiling without OpenMP support")
			if (CMAKE_CXX_COMPILER_ID STREQUAL Intel AND CMAKE_VERSION VERSION_LESS 3.9)
				message(WARNING
"**********************************************************************"
"\nIntel compilers actually support the OpenMP version required by shark,"
"but versions of cmake < 3.9 might fail to recognize such support."
"Try using cmake >= 3.9 to recognize this properly"
"and generate an OpenMP-enabled shark build.\n"
"**********************************************************************")
			endif()
		endif()
	endif()

	if (SHARK_OPENMP)
		add_compile_options("${OpenMP_CXX_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
		set(SHARK_LIBS ${SHARK_LIBS} ${OpenMP_CXX_LIBRARIES})
	endif()
endmacro()

#
# Go out there and find stuff
#
find_threads()
find_gsl()
find_hdf5()
find_boost()
if (NOT SHARK_NO_OPENMP)
	find_openmp()
endif()

# Windows builds need to link to ws2_32 (due to usage of gethostname)
if (WIN32)
	set(SHARK_LIBS ${SHARK_LIBS} ws2_32)
endif()

# Compiler-specific tweaks
if (MSVC)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

#
# Save all compile-time options here
#
set(config_h "${CMAKE_CURRENT_BINARY_DIR}/include/config.h")
configure_file(include/config.h.in "${config_h}" @ONLY)

#
# We use the source code directory and the installation prefix to lookup data
# files distributed with shark itself, so users don't need to worry about their
# location
#
set(data_cpp "${CMAKE_CURRENT_BINARY_DIR}/data.cpp")
configure_file(src/data.cpp.in "${data_cpp}" @ONLY)

#
# Add git revision information in a generated .cpp file, which we then include
# as part of the shark library
#
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_local_changes(GIT_HAS_LOCAL_CHANGES)
set(git_revision_cpp "${CMAKE_CURRENT_BINARY_DIR}/git_revision.cpp")
configure_file("src/git_revision.cpp.in" "${git_revision_cpp}" @ONLY)


# Our own include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_BINARY_DIR}/include)

# The shark shared library
set(SHARKLIB_SRCS
   "${config_h}"
   "${data_cpp}"
   "${git_revision_cpp}"
   include/agn_feedback.h
   include/baryon.h
   include/components.h
   include/cosmology.h
   include/dark_matter_halos.h
   include/data.h
   include/disk_instability.h
   include/environment.h
   include/evolve_halos.h
   include/exceptions.h
   include/execution.h
   include/galaxy.h
   include/galaxy_creator.h
   include/galaxy_mergers.h
   include/galaxy_writer.h
   include/gas_cooling.h
   include/git_revision.h
   include/halo.h
   include/integrator.h
   include/interpolator.h
   include/logging.h
   include/merger_tree.h
   include/merger_tree_reader.h
   include/mixins.h
   include/naming_convention.h
   include/nfw_distribution.h
   include/numerical_constants.h
   include/ode_solver.h
   include/omp_utils.h
   include/options.h
   include/physical_model.h
   include/recycling.h
   include/reincorporation.h
   include/reionisation.h
   include/root_solver.h
   include/shark.h
   include/shark_runner.h
   include/simulation.h
   include/star_formation.h
   include/stellar_feedback.h
   include/subhalo.h
   include/timer.h
   include/total_baryon.h
   include/tree_builder.h
   include/utils.h
   include/hdf5/attribute.h
   include/hdf5/data_set.h
   include/hdf5/data_space.h
   include/hdf5/data_type.h
   include/hdf5/entity.h
   include/hdf5/file.h
   include/hdf5/group.h
   include/hdf5/location.h
   include/hdf5/utils.h
   include/hdf5/io/iobase.h
   include/hdf5/io/reader.h
   include/hdf5/io/traits.h
   include/hdf5/io/writer.h
   src/agn_feedback.cpp
   src/cosmology.cpp
   src/execution.cpp
   src/dark_matter_halos.cpp
   src/disk_instability.cpp
   src/environment.cpp
   src/evolve_halos.cpp
   src/galaxy_creator.cpp
   src/galaxy_mergers.cpp
   src/galaxy_writer.cpp
   src/gas_cooling.cpp
   src/halo.cpp
   src/integrator.cpp
   src/interpolator.cpp
   src/logging.cpp
   src/merger_tree_reader.cpp
   src/naming_convention.cpp
   src/options.cpp
   src/ode_solver.cpp
   src/physical_model.cpp
   src/recycling.cpp
   src/reincorporation.cpp
   src/reionisation.cpp
   src/root_solver.cpp
   src/shark_runner.cpp
   src/simulation.cpp
   src/star_formation.cpp
   src/stellar_feedback.cpp
   src/subhalo.cpp
   src/total_baryon.cpp
   src/tree_builder.cpp
   src/utils.cpp
   src/hdf5/attribute.cpp
   src/hdf5/data_set.cpp
   src/hdf5/data_space.cpp
   src/hdf5/data_type.cpp
   src/hdf5/entity.cpp
   src/hdf5/file.cpp
   src/hdf5/group.cpp
   src/hdf5/location.cpp
   src/hdf5/io/iobase.cpp
   src/hdf5/io/reader.cpp
   src/hdf5/io/traits.cpp
   src/hdf5/io/writer.cpp
)

if (WIN32)
	set(SHARKLIB_MODE STATIC)
else()
	set(SHARKLIB_MODE SHARED)
endif()
add_library(sharklib ${SHARKLIB_MODE} ${SHARKLIB_SRCS})
target_link_libraries(sharklib ${SHARK_LIBS})
set_target_properties(sharklib PROPERTIES
    LIBRARY_OUTPUT_NAME shark
    RUNTIME_OUTPUT_NAME shark)

# The shark-importer executable
set(SHARK_IMPORTER_SRCS
	include/importer/descendants.h
	include/importer/reader.h
	include/importer/velociraptor.h
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
)
add_executable(shark ${SHARK_SRCS})
target_link_libraries(shark sharklib)

# Installing stuff: programs, scripts, static data
install(TARGETS sharklib shark shark-importer
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
install(PROGRAMS hpc/shark-run hpc/shark-submit
        DESTINATION bin)
install(DIRECTORY data
        DESTINATION share/shark)

#
# Unit tests
#
if( SHARK_TEST )
	find_package(CxxTest)

	# In cmake < 3.7 the cxxtestgen executable is always executed through an
	# interpreter, and doesn't honour the shebang line that could be present.
	# We try to rectify this situation here with a solution similar to what
	# is found in cmake >= 3.7
	if (${CMAKE_VERSION} VERSION_LESS 3.7 AND
	    NOT "${CXXTEST_PYTHON_TESTGEN_EXECUTABLE}" STREQUAL "")
		execute_process(COMMAND ${CXXTEST_PYTHON_TESTGEN_EXECUTABLE} --version
		                OUTPUT_VARIABLE _CXXTEST_OUT ERROR_VARIABLE _CXXTEST_OUT
		                RESULT_VARIABLE _CXXTEST_RESULT)
		if (_CXXTEST_RESULT EQUAL 0)
			set(CXXTEST_TESTGEN_INTERPRETER "")
		endif()
	endif()

	if( CXXTEST_FOUND )
		enable_testing()
		add_subdirectory(tests)
	endif()
endif()
