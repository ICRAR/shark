//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2017
// Copyright by UWA (in the framework of the ICRAR)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//

/**
 * @file
 *
 * Type traits for hdf5 data handling
 */

#ifndef SHARK_HDF5_TRAITS
#define SHARK_HDF5_TRAITS

#include <hdf5.h>
#include <H5Cpp.h>

namespace shark {

namespace hdf5 {

//
// Traits for data type conversion
// The references are set in traits.cpp
//
template <typename T>
struct datatype_traits {
};

// std::string doesn't need a native type because we actually implement its
// datawriting separately
template<>
struct datatype_traits<std::string> {
	static const H5::PredType &write_type;
};

// bool doesn't need a native type because it's actually not natively supported
// by HDF5. We instead manually convert them to int later.
template<>
struct datatype_traits<bool> {
	static const H5::PredType &write_type;
};

template<>
struct datatype_traits<float> {
	static const H5::PredType &native_type;
	static const H5::PredType &write_type;
};

template<>
struct datatype_traits<double> {
	static const H5::PredType &native_type;
	static const H5::PredType &write_type;
};

template<>
struct datatype_traits<int> {
	static const H5::PredType &native_type;
	static const H5::PredType &write_type;
};

template<>
struct datatype_traits<unsigned int> {
	static const H5::PredType &native_type;
	static const H5::PredType &write_type;
};

template<>
struct datatype_traits<std::int64_t> {
	static const H5::PredType &native_type;
	static const H5::PredType &write_type;
};

//
// Traits for HDF5 entities
//
template <H5G_obj_t E>
struct entity_traits {
	using rettype = void ;
};

template<>
struct entity_traits<H5G_GROUP> {
	using rettype = H5::Group;
};

template<>
struct entity_traits<H5G_DATASET> {
	using rettype = H5::DataSet;
};

}  // namespace hdf5

}  // namespace shark

#endif // SHARK_HDF5_TRAITS