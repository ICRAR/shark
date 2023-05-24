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
 * C++ Wrappers for HDF5 entities that are managed with identifiers (hid_t)
 */

#include <stdexcept>
#include <sstream>
#include "hdf5/entity.h"

namespace shark {
namespace hdf5 {

Resource::Resource(H5I_type_t expectedType, hid_t handle) : handle(handle) {
	if (!isValid()) {
		throw std::runtime_error("Invalid handle");
	}

	auto type = H5Iget_type(handle);
	if (type == H5I_BADID) {
		throw std::runtime_error("Unable to determine resource type");
	} else if (type != expectedType) {
		std::ostringstream os;
		os << "Expected resource of type " << expectedType << ", got " << type;
		throw std::runtime_error(os.str());
	}
}

Resource::Resource(const Resource& other) : handle(other.handle) {
	H5Iinc_ref(handle);
}

Resource& Resource::operator=(const Resource& rhs) {
	H5Idec_ref(handle);
	handle = rhs.handle;
	H5Iinc_ref(handle);
	return *this;
}

// Needed to satisfy the linker, but we want to make overriding required for subclasses
// As a result we don't put this in the class definition
// See here for more details: https://stackoverflow.com/a/11437551
Resource::~Resource() = default;

bool Resource::isValid() const {
	auto isValid = H5Iis_valid(handle);
	if (isValid < 0) {
		throw std::runtime_error("Error in H5Iis_valid");
	}

	return isValid > 0;
}

void Resource::setComment(const std::string& comment) {
	H5Oset_comment(getHandle(), comment.c_str());
}

hid_t Resource::getHandle() const {
	return handle;
}

} // shark
} // hdf5