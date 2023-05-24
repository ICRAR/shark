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
 * C++ Wrappers for HDF5 entities that are managed with identifiers (hid_t)
 */

#ifndef SHARK_HDF5_ENTITY_H
#define SHARK_HDF5_ENTITY_H

#include <string>
#include <hdf5.h>

namespace shark {
namespace hdf5 {

class Resource {
public:
	Resource(H5I_type_t expectedType, hid_t handle);
	Resource(const Resource& other);
	Resource& operator=(const Resource& rhs);
	Resource(Resource&&) = default;
	Resource& operator=(Resource&&) = default;
	virtual ~Resource();

	void setComment(const std::string& comment);

	hid_t getHandle() const;
	// TODO getName using H5Iget_name()

private:
	hid_t handle;

	bool isValid() const;
};

} // namespace hdf5
} // namespace shark

#endif //SHARK_HDF5_ENTITY_H
