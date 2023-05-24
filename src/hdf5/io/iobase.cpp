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
 * Implementation of the IOBase class for HDF5 handling
 */

#include <sstream>
#include <string>
#include <stdexcept>

#include "hdf5/utils.h"
#include "hdf5/io/iobase.h"
#include "logging.h"
#include "utils.h"

namespace shark {
namespace hdf5 {

IOBase::IOBase(const std::string& filename, const FileOpenMethod& openMethod) : hdf5_file(File(filename, openMethod)) {
}

IOBase::~IOBase() {
	close();
}

void IOBase::close() {
	if (!hdf5_file.has_value()) {
		return;
	}

	// RAII will close the file for us!
	hdf5_file = boost::none;
}

void IOBase::open_file(const std::string& filename, const FileOpenMethod& openMethod) {
	hdf5_file = File(filename, openMethod);
}

std::string IOBase::get_filename() const {
	return hdf5_file->getFileName();
}

DataSpace IOBase::get_nd_dataspace(const DataSet& dataset, unsigned int expected_ndims) const {
	DataSpace space = dataset.getSpace();
	int ndims = space.getSimpleExtentNdims();
	if (ndims != int(expected_ndims)) {
		std::ostringstream os;
		os << ndims << " dimensions found in dataset";
		os << " " << dataset.getName();
		os << ", " << expected_ndims << " expected";
		throw std::runtime_error(os.str());
	}
	return space;
}

DataSpace IOBase::get_scalar_dataspace(const DataSet& dataset) const {
	return get_nd_dataspace(dataset, 0);
}

DataSpace IOBase::get_1d_dataspace(const DataSet& dataset) const {
	return get_nd_dataspace(dataset, 1);
}

DataSpace IOBase::get_2d_dataspace(const DataSet& dataset) const {
	return get_nd_dataspace(dataset, 2);
}

hsize_t IOBase::get_1d_dimsize(const DataSpace& space) const {
	auto dimensions = space.getSimpleExtentDims();
	assert(dimensions.size() == 1);
	return dimensions[0];
}

DataSet IOBase::get_dataset(const std::string& name) const {

	LOG(debug) << "Getting dataset " << name << " on file " << get_filename();

	// The name might contain slashes, so we can navigate through
	// a hierarchy of groups/datasets
	auto parts = tokenize(name, "/");

	return get_dataset(parts);
}

DataSet IOBase::get_dataset(const std::vector<std::string>& path) const {

	// only the attribute name, read directly and come back
	if (path.size() == 1) {
		return hdf5_file->openDataSet(path[0]);
	}

	// else there's a path to follow, go for it!
	Group group = hdf5_file->openGroup(path.front());
	std::vector<std::string> group_paths(path.begin() + 1, path.end() - 1);
	for (auto const& path: group_paths) {
		LOG(debug) << "Getting dataset " << path << " on file " << get_filename();
		group = group.openGroup(path);
	}

	return group.openDataSet(path.back());
}

}  // namespace hdf5
}  // namespace shark