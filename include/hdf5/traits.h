//
// Type traits for hdf5 data handling
//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2017
// Copyright by UWA (in the framework of the ICRAR)
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307  USA
//

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

template<>
struct datatype_traits<std::string> {
	static const H5::PredType &write_type;
};

template<>
struct datatype_traits<float> {
	static const H5::PredType &write_type;
};

template<>
struct datatype_traits<double> {
	static const H5::PredType &write_type;
};

template<>
struct datatype_traits<int> {
	static const H5::PredType &write_type;
};

template<>
struct datatype_traits<unsigned int> {
	static const H5::PredType &write_type;
};


template<>
struct datatype_traits<long int> {
	static const H5::PredType &write_type;
};

template<>
struct datatype_traits<bool> {
	static const H5::PredType &write_type;
};


//
// Traits for HDF5 entities
//
template <H5G_obj_t E>
struct entity_traits {
	typedef void rettype;
};

template<>
struct entity_traits<H5G_GROUP> {
	typedef H5::Group rettype;
};

template<>
struct entity_traits<H5G_DATASET> {
	typedef H5::DataSet rettype;
};

}  // namespace hdf5

}  // namespace shark

#endif // SHARK_HDF5_TRAITS