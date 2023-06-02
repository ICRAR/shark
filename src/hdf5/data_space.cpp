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

#include "logging.h"
#include "hdf5/data_space.h"
#include "hdf5/data_set.h"
#include "hdf5/utils.h"

namespace shark {
namespace hdf5 {

DataSpace::DataSpace(hid_t handle) : Entity(H5I_DATASPACE, handle) {}

DataSpace::DataSpace(const DataSet& dataSet) : Entity(H5I_DATASPACE, H5Dget_space(dataSet.getId())) {
}

DataSpace::~DataSpace() {
	if (H5Sclose(getId()) < 0) {
		LOG(error) << "H5Sclose() failed";
	};
}

DataSpace DataSpace::create(const DataSpaceType& type) {
	return DataSpace(H5Screate(static_cast<H5S_class_t>(type)));
}

DataSpace DataSpace::create(std::vector<hsize_t> dimensions) {
	auto rank = static_cast<int>(dimensions.size());
	return DataSpace(H5Screate_simple(rank, dimensions.data(), nullptr));
}

int DataSpace::getSimpleExtentNdims() const {
	return H5Sget_simple_extent_ndims(getId());
}

std::vector<hsize_t> DataSpace::getSimpleExtentDims() const {
	auto ndims = getSimpleExtentNdims();
	std::vector<hsize_t> dim_sizes(ndims);
	H5Sget_simple_extent_dims(getId(), dim_sizes.data(), nullptr);
	return dim_sizes;
}

std::vector<hsize_t> DataSpace::getSimpleExtentMaxDims() const {
	auto ndims = getSimpleExtentNdims();
	std::vector<hsize_t> dim_sizes(ndims);
	std::vector<hsize_t> max_dim_sizes(ndims);
	H5Sget_simple_extent_dims(getId(), dim_sizes.data(), max_dim_sizes.data());
	return max_dim_sizes;
}

void DataSpace::selectHyperslab(const HyperslabSelection& op, const std::vector<hsize_t>& start,
                                const std::vector<hsize_t>& stride, const std::vector<hsize_t>& count,
                                const std::vector<hsize_t>& block) {
	// Assert all same length & same as rank
	if (H5Sselect_hyperslab(getId(), static_cast<H5S_seloper_t>(op), start.data(), stride.data(), count.data(),
	                        block.data()) < 0) {
		throw hdf5_api_error("H5Sselect_hyperslab");
	}
}

} // namespace hdf5
} // namespace shark