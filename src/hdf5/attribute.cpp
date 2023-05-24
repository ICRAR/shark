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
 * C++ wrapper for HDF5 attributes
 */

#include "logging.h"
#include "hdf5/attribute.h"

namespace shark {
namespace hdf5 {

Attribute::Attribute(hid_t handle) : Entity(H5I_ATTR, handle) {}

Attribute::Attribute(const Location& group, const std::string& name) :
		Attribute(H5Aopen(group.getId(), name.c_str(), H5P_DEFAULT)) {
}

Attribute::~Attribute() {
	if (H5Aclose(getId()) < 0) {
		LOG(error) << "H5Aclose() failed";
	}
}

Attribute Attribute::create(Location& location, const std::string& name, const DataType& dataType,
                            const DataSpace& dataSpace) {
	return Attribute(
			H5Acreate(location.getId(), name.c_str(), dataType.getId(), dataSpace.getId(), H5P_DEFAULT,
			          H5P_DEFAULT));
}

DataType Attribute::getType() const {
	// H5Aget_type() needs to be closed, so we wrap this in DataType rather
	// than using it directly
	return DataType(H5Aget_type(getId()));
}

template<>
std::string Attribute::read<std::string>() const {
	// Unfortunately different enough that we can't use stringFromHdf5Api
	H5A_info_t info;
	if (H5Aget_info(getId(), &info) < 0) {
		throw hdf5_api_error("H5Aget_info", "Unable to get attribute " + getName() + "info");
	}

	std::string val;
	val.resize(info.data_size); // Size includes null-terminator

	if (H5Aread(getId(), getType().getId(), &val[0]) < 0) {
		throw hdf5_api_error("H5Aread", "Unable to read attribute " + getName());
	}

	val.resize(info.data_size - 1); // Trim included null-terminator
	return val;
}

template<>
void Attribute::write<std::string>(const DataType& dataType, const std::string& val) {
	if (H5Awrite(getId(), dataType.getId(), val.data()) < 0) {
		throw hdf5_api_error("H5Awrite", "Unable to write to attribute " + getName());
	}
}

} // namespace hdf5
} // namespace shark

