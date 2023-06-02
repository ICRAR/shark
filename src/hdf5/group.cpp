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
 * C++ wrappers for dealing with HDF5 groups (abstract or concrete)
 */

#include "logging.h"
#include "hdf5/group.h"
#include "hdf5/data_set.h"
#include "hdf5/utils.h"

namespace shark {
namespace hdf5 {

/*** AbstractGroup ***/

AbstractGroup::AbstractGroup(H5I_type_t expectedType, hid_t handle) : Location(expectedType, handle) {}

hsize_t AbstractGroup::getNumObjs() const {
	hsize_t num;
	// TODO Deprecated in 1.8.0, replacement is H5Gget_info()
	if (H5Gget_num_objs(getId(), &num) < 0) {
		throw hdf5_api_error("H5Gget_num_objs");
	}
	return num;
}

std::string AbstractGroup::getObjnameByIdx(hsize_t idx) const {
	auto id = getId();
	return stringFromHdf5Api("H5Gget_objname_by_idx", [id, idx](char* buf, size_t size) {
		// TODO Deprecated in 1.8.0, replacement is H5Lget_name_by_idx()
		return H5Gget_objname_by_idx(id, idx, buf, size);
	});
}

H5G_obj_t AbstractGroup::getObjTypeByIdx(hsize_t idx) const {
	// TODO Deprecated in 1.8.0, replacement is H5Oget_info()
	auto type = H5Gget_objtype_by_idx(getId(), idx);
	if (type < 0) {
		throw hdf5_api_error("H5Gget_objtype_by_idx");
	}
	return type;
}

Group AbstractGroup::openGroup(const std::string& name) const {
	return Group(*this, name);
}

DataSet AbstractGroup::openDataSet(const std::string& name) const {
	return DataSet(*this, name);
}

Group AbstractGroup::createGroup(const std::string& name) {
	return Group::create(*this, name);
}

DataSet AbstractGroup::createDataSet(const std::string& name, const DataType& dataType, const DataSpace& dataSpace) {
	return DataSet::create(*this, name, dataType, dataSpace);
}

/*** Group ***/

Group::Group(hid_t handle) : AbstractGroup(H5I_GROUP, handle) {}

Group::Group(const AbstractGroup& parent, const std::string& name) :
		Group(H5Gopen(parent.getId(), name.c_str(), H5P_DEFAULT)) {
}

Group::~Group() {
	if (H5Gclose(getId()) < 0) {
		LOG(error) << "H5Gclose() failed";
	}
}

Group Group::create(AbstractGroup& parent, const std::string& name) {
	return Group(H5Gcreate(parent.getId(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
}

} // namespace hdf5
} // namespace shark
