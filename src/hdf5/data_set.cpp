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

#include "hdf5/data_set.h"
#include "hdf5/data_space.h"
#include "hdf5/group.h"
#include "hdf5/data_type.h"

namespace shark {
namespace hdf5 {

DataSet::DataSet(std::string name, hid_t handle) : Location(H5I_DATASET, handle), objName(std::move(name)) {}

DataSet::DataSet(const AbstractGroup& file, const std::string& name) :
		Location(H5I_DATASET, H5Dopen(file.getHandle(), name.c_str(), H5P_DEFAULT)),
		objName(name) {
}

DataSet::~DataSet() {
	H5Dclose(getHandle());
}

DataSet DataSet::create(shark::hdf5::AbstractGroup& parent, const std::string& name, const DataType& dataType,
                        const shark::hdf5::DataSpace& dataSpace) {
	return DataSet(name,
	               H5Dcreate(parent.getHandle(), name.c_str(), dataType.getHandle(), dataSpace.getHandle(), H5P_DEFAULT,
	                         H5P_DEFAULT, H5P_DEFAULT));
}

DataSpace DataSet::getSpace() const {
	return DataSpace(*this);
}

const std::string& DataSet::getObjName() const {
	return objName;
}

hid_t DataSet::getDataType() const {
	return H5Dget_type(getHandle());
}

herr_t DataSet::read(void *buf, hid_t dataType, const DataSpace& memSpace, const DataSpace& fileSpace) const {
	return H5Dread(getHandle(), dataType, memSpace.getHandle(), fileSpace.getHandle(), H5P_DEFAULT, buf);
}

herr_t
DataSet::write(const void *buf, const DataType& memDataType, const DataSpace& memSpace, const DataSpace& fileSpace) {
	return H5Dwrite(getHandle(), memDataType.getHandle(), memSpace.getHandle(), fileSpace.getHandle(), H5P_DEFAULT,
	                buf);
}

} // namespace hdf5
} // namespace shark