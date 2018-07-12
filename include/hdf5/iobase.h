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
 * Declaration of the IOBase class for HDF5 handling
 */

#ifndef INCLUDE_HDF5_IOBASE_H_
#define INCLUDE_HDF5_IOBASE_H_

#include <vector>

#include <H5Cpp.h>

// Define handy macros to detect whether we are above 1.8.11 and/or 1.10.0
// These versions introduce some important backward-incompatible changes in the
// C++ API that we need to be aware of if we want to support these versions
#undef HDF5_NEWER_THAN_1_8_11
#undef HDF5_NEWER_THAN_1_10_0
#if HDF5_VERSION_MAJOR == 1 && \
     (HDF5_VERSION_MINOR > 10 || \
      (HDF5_VERSION_MINOR == 10 && HDF5_VERSION_PATCH >= 1))
#define HDF5_NEWER_THAN_1_10_0
#endif
#if HDF5_VERSION_MAJOR == 1 && \
     (HDF5_VERSION_MINOR > 8 || \
      (HDF5_VERSION_MINOR == 8 && HDF5_VERSION_PATCH >= 12))
#define HDF5_NEWER_THAN_1_8_11
#endif

namespace shark {

namespace hdf5 {

/**
 * Base class for objects handling HDF5 I/O
 */
class IOBase {

public:

	/**
	 * Creates a new IOBase instance without opening a file
	 */
	IOBase();

	/**
	 * Opens the given file in the given mode
	 *
	 * @param filename The HDF5 filename
	 * @param flags The mode in which the file will be opened
	 */
	IOBase(const std::string &filename, unsigned int flags);

	/**
	 * Closes the file and destroys this class
	 */
	~IOBase();

	/**
	 * Closes the file
	 */
	void close();

	/**
	 * Returns the filename being handled by this class
	 * @return The filename being handled by this class
	 */
	const std::string get_filename() const;

	/**
	 * Opens the given file in the given mode
	 *
	 * @param filename the HDF5 filename
	 * @param flags The mode in which the file will be opened
	 */
	void open_file(const std::string &filename, unsigned int flags);

protected:

	H5::DataSet get_dataset(const std::string &name) const;
	H5::DataSet get_dataset(const std::vector<std::string> &path) const;
	H5::DataSpace get_scalar_dataspace(const H5::DataSet &dataset) const;
	H5::DataSpace get_1d_dataspace(const H5::DataSet &dataset) const;
	H5::DataSpace get_2d_dataspace(const H5::DataSet &dataset) const;
	hsize_t get_1d_dimsize(const H5::DataSpace &space) const;

	H5::H5File hdf5_file;

private:

	bool opened;
	H5::DataSpace get_nd_dataspace(const H5::DataSet &dataset, unsigned int expected_ndims) const;
};

}  // namespace hdf5

}  // namespace shark



#endif /* INCLUDE_HDF5_IOBASE_H_ */