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

#include <string>
#include <vector>
#include<boost/optional.hpp>

#include "hdf5/data_set.h"
#include "hdf5/data_space.h"
#include "hdf5/file.h"

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
	IOBase() = default;

	/**
	 * Opens the given file in the given mode
	 *
	 * @param filename The HDF5 filename
	 * @param flags The mode in which the file will be opened
	 */
	IOBase(const std::string& filename, const FileOpenMethod& openMethod);

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
	std::string get_filename() const;

	/**
	 * Opens the given file in the given mode
	 *
	 * @param filename the HDF5 filename
	 * @param flags The mode in which the file will be opened
	 */
	void open_file(const std::string& filename, const FileOpenMethod& openMethod);

protected:
	DataSet get_dataset(const std::string& name) const;
	DataSet get_dataset(const std::vector<std::string>& path) const;
	DataSpace get_scalar_dataspace(const DataSet& dataset) const;
	DataSpace get_1d_dataspace(const DataSet& dataset) const;
	DataSpace get_2d_dataspace(const DataSet& dataset) const;
	hsize_t get_1d_dimsize(const DataSpace& space) const;

	boost::optional<File> hdf5_file;

private:
	DataSpace get_nd_dataspace(const DataSet& dataset, unsigned int expected_ndims) const;
};

}  // namespace hdf5

}  // namespace shark



#endif /* INCLUDE_HDF5_IOBASE_H_ */