//
// HHDF5-realted traits
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

#include <string>

#include "hdf5/traits.h"

namespace shark {

namespace hdf5 {

const H5::PredType &datatype_traits<std::string>::write_type = H5::PredType::C_S1;
const H5::PredType &datatype_traits<float>::write_type = H5::PredType::NATIVE_FLOAT;
const H5::PredType &datatype_traits<double>::write_type = H5::PredType::NATIVE_DOUBLE;
const H5::PredType &datatype_traits<int>::write_type = H5::PredType::NATIVE_INT16;
const H5::PredType &datatype_traits<unsigned int>::write_type = H5::PredType::NATIVE_UINT16;
const H5::PredType &datatype_traits<long int>::write_type = H5::PredType::NATIVE_INT32;
const H5::PredType &datatype_traits<bool>::write_type = H5::PredType::NATIVE_INT8;

}  // namespace hdf5

}  // namespace shark