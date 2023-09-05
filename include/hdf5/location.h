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

#ifndef SHARK_HDF5_LOCATION_H
#define SHARK_HDF5_LOCATION_H

#include <hdf5.h>
#include "hdf5/entity.h"
#include "hdf5/data_type.h"
#include "hdf5/data_space.h"

namespace shark {
namespace hdf5 {

// Break circular reference
class Attribute;

class Location : public Entity {
public:
	Location(H5I_type_t expectedType, hid_t handle);

	std::string getFileName() const;
	void setComment(const std::string& comment);

	bool attributeExists(const std::string& name) const;
	Attribute openAttribute(const std::string& name) const;
	Attribute createAttribute(const std::string& name, const DataType& dataType, const DataSpace& dataSpace);
};

} // namespace hdf5
} // namespace shark

#endif //SHARK_HDF5_LOCATION_H
