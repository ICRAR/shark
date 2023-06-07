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
 * C++ Wrappers for HDF5 data types
 */

#include "logging.h"
#include "hdf5/data_type.h"
#include "hdf5/utils.h"

namespace shark {
namespace hdf5 {

DataType::DataType(hid_t handle) : Entity(H5I_DATATYPE, handle) {
}

DataType::~DataType() {
	if (H5Tclose(getId()) < 0) {
		LOG(error) << "H5Tclose() failed";
	}
}

// Always copy the predefined type so that the deconstructor will succeed!
// (H5Tclose() will fail for the predefined types as they are considered immutable)
PredefinedDataType::PredefinedDataType(hid_t handle) : DataType(H5Tcopy(handle)) {}

const PredefinedDataType& PredefinedDataType::C_S1() {
	static const PredefinedDataType t(H5T_C_S1);
	return t;
}

const PredefinedDataType& PredefinedDataType::NATIVE_INT8() {
	static const PredefinedDataType t(H5T_NATIVE_INT8);
	return t;
}

const PredefinedDataType& PredefinedDataType::NATIVE_FLOAT() {
	static const PredefinedDataType t(H5T_NATIVE_FLOAT);
	return t;
}

const PredefinedDataType& PredefinedDataType::NATIVE_DOUBLE() {
	static const PredefinedDataType t(H5T_NATIVE_DOUBLE);
	return t;
}

const PredefinedDataType& PredefinedDataType::NATIVE_INT() {
	static const PredefinedDataType t(H5T_NATIVE_INT);
	return t;
}

const PredefinedDataType& PredefinedDataType::NATIVE_INT32() {
	static const PredefinedDataType t(H5T_NATIVE_INT32);
	return t;
}

const PredefinedDataType& PredefinedDataType::NATIVE_UINT() {
	static const PredefinedDataType t(H5T_NATIVE_UINT);
	return t;
}

const PredefinedDataType& PredefinedDataType::NATIVE_UINT32() {
	static const PredefinedDataType t(H5T_NATIVE_UINT32);
	return t;
}

const PredefinedDataType& PredefinedDataType::NATIVE_INT64() {
	static const PredefinedDataType t(H5T_NATIVE_INT64);
	return t;
}

StringDataType::StringDataType(size_t size) : DataType(H5Tcopy(PredefinedDataType::C_S1().getId())) {
	if (H5Tset_size(getId(), size) < 0) {
		throw hdf5_api_error("H5Tset_size");
	}
}

} // namespace hdf5
} // namespace shark