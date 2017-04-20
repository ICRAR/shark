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

namespace shark {

namespace hdf5 {

template <typename T>
struct type_traits {
	const static hid_t preferred_type_id = -100;
	const static hid_t type_ids[] = {};
};

template <>
struct type_traits<float> {
#ifdef __BIG_ENDIAN
	const static hid_t preferred_type_id = H5T_IEEE_F32BE;
#else
	const static hid_t preferred_type_id = H5T_IEEE_F32LE;
#endif // __BIG_ENDIAN
	const static hid_t type_ids[] = {H5T_IEEE_F32BE, H5T_IEEE_F32LE};
};

template <>
struct type_traits<double> {
#ifdef __BIG_ENDIAN
	const static hid_t preferred_type_id = H5T_IEEE_F64BE;
#else
	const static hid_t preferred_type_id = H5T_IEEE_F64LE;
#endif // __BIG_ENDIAN
	const static hid_t type_ids[] = {H5T_IEEE_F64BE, H5T_IEEE_F64LE};
};

}  // namespace hdf5

}  // namespace shark

#endif // SHARK_HDF5_TRAITS