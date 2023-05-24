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
 * Various utilities for interacting with the HDF5 API
 */

#ifndef SHARK_HDF5_UTILS_H
#define SHARK_HDF5_UTILS_H

#include <string>
#include <stdexcept>
#include <hdf5.h>
#include "exceptions.h"

namespace shark {
namespace hdf5 {

class hdf5_api_error : public exception {
public:
	explicit hdf5_api_error(const std::string& function_name) : exception(
			function_name + "() failure") {};

	hdf5_api_error(const std::string& function_name, const std::string& message) : exception(
			function_name + "() failure: " + message) {};
};

/**
 * Retrieve a string from a HDF5 api call that adheres to the
 * "pass nullptr to get size, then pass buffer with size to write string"
 * pattern.
 *
 * The passed in callback will first be called with (nullptr, 0) and is
 * expected to return the required size of a buffer to hold the true value
 * (excluding null-terminator). Then it will be called with a pointer to
 * the buffer to write to, plus the previously returned size+1 (to account for
 * the null-terminator.
 * @tparam F Type of callable
 * @param f Callable of signature (char*, size_t) -> ssize_t
 * @return The correctly allocated string holding the entire value from the HDF5 api
 */
template<typename F>
std::string stringFromHdf5Api(const std::string& hdf5_fn_name, F&& f) {
	ssize_t size = f(nullptr, 0);
	if (size < 0) {
		throw hdf5_api_error(hdf5_fn_name, "Unable to determine size of string");
	}

	// Ensure there's enough space for the null-terminator that the C API will write
	// (even though C++11 will allocate enough for a null-terminator, it's UB to
	// access it)
	std::string s;
	s.resize(size + 1);

	size = f(&s[0], s.size());
	if (size < 0) {
		throw hdf5_api_error(hdf5_fn_name, "Failed to write string to buffer");
	}

	// Now resize back down to the actual size
	s.resize(size);

	return s;
}

} // namespace hdf5
} // namespace shark

#endif //SHARK_HDF5_UTILS_H
