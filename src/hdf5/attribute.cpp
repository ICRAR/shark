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

#include "hdf5/attribute.h"

namespace shark {
namespace hdf5 {

Attribute::Attribute(hid_t handle) : Resource(H5I_ATTR, handle) {}

Attribute::Attribute(const Location& group, const std::string& name) :
		Attribute(H5Aopen(group.getHandle(), name.c_str(), H5P_DEFAULT)) {
}

Attribute::~Attribute() {
	H5Aclose(getHandle());
}

Attribute Attribute::create(Location& location, const std::string& name, const DataType& dataType,
                            const DataSpace& dataSpace) {
	return Attribute(
			H5Acreate(location.getHandle(), name.c_str(), dataType.getHandle(), dataSpace.getHandle(), H5P_DEFAULT,
			          H5P_DEFAULT));
}

template<>
std::string Attribute::read<std::string>() const {
	H5A_info_t info;
	H5Aget_info(getHandle(), &info);
	std::string val;
	val.resize(info.data_size); // Size includes null-terminator
	auto type = H5Aget_type(getHandle());
	H5Aread(getHandle(), type, &val[0]);
	val.resize(info.data_size - 1); // Trim included null-terminator
	return val;
}

template<>
void Attribute::write<std::string>(const DataType& dataType, const std::string& val) {
	H5Awrite(getHandle(), dataType.getHandle(), val.data());
}

} // namespace hdf5
} // namespace shark

