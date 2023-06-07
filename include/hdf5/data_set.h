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

#ifndef SHARK_HDF5_DATA_SET_H
#define SHARK_HDF5_DATA_SET_H

#include <hdf5.h>
#include "hdf5/location.h"

namespace shark {
namespace hdf5 {

class AbstractGroup;

class DataSet : public Location {
public:
	DataSet(const AbstractGroup& file, const std::string& name);
	~DataSet() override;

	static DataSet
	create(AbstractGroup& parent, const std::string& name, const DataType& dataType, const DataSpace& dataSpace);

	DataType getDataType() const;
	DataSpace getSpace() const;

	void read(void* buf, const DataType& dataType, const DataSpace& memSpace, const DataSpace& fileSpace) const;
	void write(const void* buf, const DataType& memDataType, const DataSpace& memSpace, const DataSpace& fileSpace);

private:
	explicit DataSet(hid_t handle);
};

} // namespace hdf5
} // namespace shark

#endif //SHARK_HDF5_DATA_SET_H
