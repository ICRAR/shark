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

#ifndef SHARK_HDF5_DATA_TYPE_H
#define SHARK_HDF5_DATA_TYPE_H

#include <hdf5.h>
#include "hdf5/entity.h"

namespace shark {
namespace hdf5 {


class DataType : public Entity {
public:
	explicit DataType(hid_t handle);
	~DataType() override;
};

class PredefinedDataType : public DataType {
public:
	static const PredefinedDataType& C_S1();
	static const PredefinedDataType& NATIVE_INT8();
	static const PredefinedDataType& NATIVE_FLOAT();
	static const PredefinedDataType& NATIVE_DOUBLE();
	static const PredefinedDataType& NATIVE_INT();
	static const PredefinedDataType& NATIVE_INT32();
	static const PredefinedDataType& NATIVE_UINT();
	static const PredefinedDataType& NATIVE_UINT32();
	static const PredefinedDataType& NATIVE_INT64();

private:
	explicit PredefinedDataType(hid_t handle);
};

class StringDataType : public DataType {
public:
	explicit StringDataType(size_t size = H5T_VARIABLE);
};

} // namespace hdf5
} // namespace shark

#endif //SHARK_HDF5_DATA_TYPE_H
