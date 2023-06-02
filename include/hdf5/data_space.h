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
 * C++ wrappers for dealing with HDF5 dataspaces
 */

#ifndef SHARK_HDF5_DATA_SPACE_H
#define SHARK_HDF5_DATA_SPACE_H

#include <vector>
#include <hdf5.h>
#include "hdf5/entity.h"

namespace shark {
namespace hdf5 {

class DataSet;

enum class DataSpaceType {
	Null = H5S_NULL,
	Scalar = H5S_SCALAR,
	Simple = H5S_SIMPLE,
};

enum class HyperslabSelection {
	// More available
	// See https://docs.hdfgroup.org/hdf5/v1_14/group___h5_s.html#ga6adfdf1b95dc108a65bf66e97d38536d
	Set = H5S_SELECT_SET,
};

class DataSpace : public Entity {
public:
	explicit DataSpace(const DataSet& dataSet);
	~DataSpace() override;

	static DataSpace create(const DataSpaceType& type);
	static DataSpace create(std::vector<hsize_t> dimensions);

	int getSimpleExtentNdims() const;
	std::vector<hsize_t> getSimpleExtentDims() const;
	std::vector<hsize_t> getSimpleExtentMaxDims() const;
	void
	selectHyperslab(const HyperslabSelection& op, const std::vector<hsize_t>& start, const std::vector<hsize_t>& stride,
	                const std::vector<hsize_t>& count, const std::vector<hsize_t>& block);

private:
	explicit DataSpace(hid_t handle);
};

} // namespace hdf5
} // namespace shark

#endif //SHARK_HDF5_DATA_SPACE_H
