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
 * C++ wrappers for abstract HDF5 "locations" (files, groups & datasets)
 */

#include <stdexcept>
#include "hdf5/attribute.h"
#include "hdf5/location.h"
#include "hdf5/utils.h"

namespace shark {
namespace hdf5 {

Location::Location(H5I_type_t expectedType, hid_t handle) : Entity(expectedType, handle) {}

std::string Location::getFileName() const {
	auto id = getId();
	return stringFromHdf5Api("H5Fget_name", [id](char* buf, ssize_t size) {
		return H5Fget_name(id, buf, size);
	});
}

void Location::setComment(const std::string& comment) {
	if (H5Oset_comment(getId(), comment.c_str()) < 0) {
		throw hdf5_api_error("H5Oset_comment");
	}
}

bool Location::attributeExists(const std::string& name) const {
	auto exists = H5Aexists(getId(), name.c_str());
	if (exists < 0) {
		throw hdf5_api_error("H5Aexists");
	}
	return exists > 0;
}

Attribute Location::openAttribute(const std::string& name) const {
	return Attribute(*this, name);
}

Attribute Location::createAttribute(const std::string& name, const DataType& dataType, const DataSpace& dataSpace) {
	return Attribute::create(*this, name, dataType, dataSpace);
}

} // namespace hdf5
} // namespace shark
