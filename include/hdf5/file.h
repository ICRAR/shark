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

#ifndef SHARK_HDF5_FILE_H
#define SHARK_HDF5_FILE_H

#include "hdf5/group.h"

namespace shark {
namespace hdf5 {

enum class FileOpenMethod {
	// Open file for reading only
	Read,
	// Open file for reading and writing
	ReadWrite,
	// Create empty file, deleting it if it already exists
	Truncate,
	// Create empty file, failing if it already exists
	ExclusiveWrite,
};

class File : public AbstractGroup {
public:
	File(const std::string& filename, const FileOpenMethod& openMethod);
	~File() override;

private:
	static hid_t openOrCreate(const std::string& filename, const FileOpenMethod& openMethod);
};

} // namespace hdf5
} // namespace shark

#endif //SHARK_HDF5_FILE_H
