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
 * C++ wrapper for dealing with HDF5 files
 */

#include <stdexcept>
#include "logging.h"
#include "hdf5/file.h"

namespace shark {
namespace hdf5 {

File::File(const std::string& filename, const FileOpenMethod& openMethod) :
		AbstractGroup(H5I_FILE, openOrCreate(filename, openMethod)) {
}

File::~File() {
	if (H5Fclose(getId()) < 0) {
		LOG(error) << "H5Fclose() failed";
	}
}

hid_t File::openOrCreate(const std::string& filename, const FileOpenMethod& openMethod) {
	switch (openMethod) {
		case FileOpenMethod::Read:
			return H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		case FileOpenMethod::ReadWrite:
			return H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		case FileOpenMethod::Truncate:
			return H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		case FileOpenMethod::ExclusiveWrite:
			return H5Fcreate(filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
		default:
			throw std::runtime_error("Unknown open type");
	}
}

} // namespace hdf5
} // namespace shark