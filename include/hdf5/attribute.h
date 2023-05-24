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

#ifndef SHARK_HDF5_ATTRIBUTE_H
#define SHARK_HDF5_ATTRIBUTE_H

#include <hdf5.h>
#include "hdf5/entity.h"
#include "hdf5/location.h"
#include "hdf5/data_space.h"
#include "hdf5/data_type.h"

namespace shark {
namespace hdf5 {

class Attribute : public Resource {
public:
	Attribute(const Location& location, const std::string& name);
	~Attribute() override;

	static Attribute
	create(Location& location, const std::string& name, const DataType& dataType, const DataSpace& dataSpace);

	template<typename T>
	T read() const {
		// Need to close or leak
		hid_t type = H5Aget_type(getHandle());
		T val;
		H5Aread(getHandle(), type, &val);
		return val;
	}

	template<typename T>
	void write(const DataType& dataType, const T& val) {
		H5Awrite(getHandle(), dataType.getHandle(), &val);
	}

private:
	explicit Attribute(hid_t handle);
};

template<>
std::string Attribute::read<std::string>() const;

template<>
void Attribute::write<std::string>(const DataType& dataType, const std::string& val);

} // namespace hdf5
} // namespace shark

#endif //SHARK_HDF5_ATTRIBUTE_H
