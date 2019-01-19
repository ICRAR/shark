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
 * HDF5-related traits
 */

#include <string>

#include "hdf5/traits.h"

namespace shark {

namespace hdf5 {

const H5::PredType &datatype_traits<std::string>::write_type = H5::PredType::C_S1;
const H5::PredType &datatype_traits<bool>::write_type = H5::PredType::NATIVE_INT8;

const H5::PredType &datatype_traits<float>::native_type = H5::PredType::NATIVE_FLOAT;
const H5::PredType &datatype_traits<float>::write_type = H5::PredType::NATIVE_FLOAT;
const H5::PredType &datatype_traits<double>::native_type = H5::PredType::NATIVE_DOUBLE;
const H5::PredType &datatype_traits<double>::write_type = H5::PredType::NATIVE_DOUBLE;
const H5::PredType &datatype_traits<int>::native_type = H5::PredType::NATIVE_INT;
const H5::PredType &datatype_traits<int>::write_type = H5::PredType::NATIVE_INT32;
const H5::PredType &datatype_traits<unsigned int>::native_type = H5::PredType::NATIVE_UINT;
const H5::PredType &datatype_traits<unsigned int>::write_type = H5::PredType::NATIVE_UINT32;
const H5::PredType &datatype_traits<std::int64_t>::native_type = H5::PredType::NATIVE_INT64;
const H5::PredType &datatype_traits<std::int64_t>::write_type = H5::PredType::NATIVE_INT64;

}  // namespace hdf5

}  // namespace shark