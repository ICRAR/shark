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
 * C++ wrappers for dealing HDF5 groups (abstract or concrete)
 */

#ifndef SHARK_HDF5_GROUP_H
#define SHARK_HDF5_GROUP_H

#include "hdf5/location.h"

namespace shark {
namespace hdf5 {

class Group;

class DataSet;

class AbstractGroup : public Location {
public:
	AbstractGroup(H5I_type_t expectedType, hid_t handle);

	hsize_t getNumObjs() const;
	std::string getObjnameByIdx(hsize_t idx) const;
	H5G_obj_t getObjTypeByIdx(hsize_t idx) const;
	Group openGroup(const std::string& name) const;
	DataSet openDataSet(const std::string& name) const;
	Group createGroup(const std::string& name);
	DataSet createDataSet(const std::string& name, const DataType& dataType, const DataSpace& dataSpace);
};

class Group : public AbstractGroup {
public:
	Group(const AbstractGroup& parent, const std::string& name);
	~Group() override;

	static Group create(AbstractGroup& parent, const std::string& name);

private:
	explicit Group(hid_t handle);
};

} // namespace hdf5
} // namespace shark

#endif //SHARK_HDF5_GROUP_H
