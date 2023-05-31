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
#include <cassert>
#include "hdf5/entity.h"
#include "hdf5/utils.h"

namespace shark {
namespace hdf5 {

Entity::Entity(H5I_type_t expectedType, hid_t id) : id(id) {
	if (!isValid()) {
		throw std::runtime_error("Invalid id");
	}

	auto type = H5Iget_type(id);
	if (type == H5I_BADID) {
		throw std::runtime_error("Unable to determine resource type");
	} else if (type != expectedType) {
		std::ostringstream os;
		os << "Expected resource of type " << expectedType << ", got " << type;
		throw std::runtime_error(os.str());
	}
}

Entity::Entity(const Entity& other) : id(other.id) {
	incrementReference();
}

Entity& Entity::operator=(const Entity& rhs) {
	decrementReference();
	id = rhs.id;
	incrementReference();
	return *this;
}

// Needed to satisfy the linker, but we want to make overriding required for subclasses
// As a result we don't put this in the class definition
// See here for more details: https://stackoverflow.com/a/11437551
Entity::~Entity() = default;

bool Entity::isValid() const {
	auto isValid = H5Iis_valid(id);
	if (isValid < 0) {
		throw std::runtime_error("Error in H5Iis_valid");
	}

	return isValid > 0;
}

void Entity::incrementReference() {
	if (H5Iinc_ref(id) < 0) {
		throw hdf5_api_error("H5Iinc_ref");
	}
}

void Entity::decrementReference() {
	if (H5Idec_ref(id) < 0) {
		throw hdf5_api_error("H5Idec_ref");
	}
}

hid_t Entity::getId() const {
	return id;
}

std::string Entity::getName() const {
	auto _id = getId();
	return stringFromHdf5Api("H5Iget_name", [_id](char* buf, size_t size) {
		return H5Iget_name(_id, buf, size);
	});
}

} // shark
} // hdf5