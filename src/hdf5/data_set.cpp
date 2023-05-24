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
 * C++ wrappers for dealing with HDF5 datasets
 */

#include "logging.h"
#include "hdf5/data_set.h"
#include "hdf5/data_space.h"
#include "hdf5/group.h"
#include "hdf5/data_type.h"
#include "hdf5/utils.h"

namespace shark {
namespace hdf5 {

DataSet::DataSet(hid_t handle) : Location(H5I_DATASET, handle) {}

DataSet::DataSet(const AbstractGroup& file, const std::string& name) :
		Location(H5I_DATASET, H5Dopen2(file.getId(), name.c_str(), H5P_DEFAULT)) {
}

DataSet::~DataSet() {
	if (H5Dclose(getId()) < 0) {
		LOG(error) << "H5Dclose() failed";
	}
}

DataSet DataSet::create(shark::hdf5::AbstractGroup& parent, const std::string& name, const DataType& dataType,
                        const shark::hdf5::DataSpace& dataSpace) {
	return DataSet(H5Dcreate2(parent.getId(), name.c_str(), dataType.getId(), dataSpace.getId(), H5P_DEFAULT,
	                          H5P_DEFAULT, H5P_DEFAULT));
}

DataSpace DataSet::getSpace() const {
	return DataSpace(*this);
}

DataType DataSet::getDataType() const {
	return DataType(H5Dget_type(getId()));
}

void DataSet::read(void* buf, const DataType& dataType, const DataSpace& memSpace, const DataSpace& fileSpace) const {
	if (H5Dread(getId(), dataType.getId(), memSpace.getId(), fileSpace.getId(), H5P_DEFAULT, buf) < 0) {
		throw hdf5_api_error("H5Dread", "Unable to read from dataset " + getName());
	}
}

void
DataSet::write(const void* buf, const DataType& memDataType, const DataSpace& memSpace, const DataSpace& fileSpace) {
	if (H5Dwrite(getId(), memDataType.getId(), memSpace.getId(), fileSpace.getId(), H5P_DEFAULT,
	             buf) < 0) {
		throw hdf5_api_error("H5Dwrite", "Unable to write to dataset " + getName());
	}
}

} // namespace hdf5
} // namespace shark