//
// Header file for the h5reader class
//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2017
// Copyright by UWA (in the framework of the ICRAR)
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307  USA
//

#ifndef SHARK_HDF5_WRITER_H_
#define SHARK_HDF5_WRITER_H_

#include <exception>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>

#include <H5Cpp.h>

#include "utils.h"
#include "hdf5/iobase.h"

namespace shark {

namespace hdf5 {

class Writer : public IOBase {

public:

	/**
	 * Constructs a new Writer object.
	 *
	 * @param filename The name of the HDF5 file to write
	 * @param overwrite Whether existing files should be overwritten or not
	 */
	Writer(const std::string &filename, bool overwrite = true) :
		IOBase(filename, H5F_ACC_RDWR | (overwrite ? H5F_ACC_CREAT|H5F_ACC_TRUNC : H5F_ACC_EXCL)) {}

	template<typename T>
	void write_attribute(const std::string &name, const T &value) const {

		// The name might contains slashes, so we can navigate through
		// a hierarchy of groups/datasets
		std::vector<std::string> parts = tokenize(name, "/");

		// only the attribute name, read directly and come back
		if( parts.size() == 1 ) {
			return _write_attribute<T>(hdf5_file, name, value);
		}

		// else there's a path to follow, go for it!
		const H5::CommonFG &location = hdf5_file;
		std::vector<std::string> path_parts(parts.begin(), parts.end()-1);
		for(auto const &path: path_parts) {
			// not implemented yet
		}

		throw std::runtime_error("write_attribute still not implemented for attributes in groups/datasets");
	}

	template<typename T>
	void write_dataset(const std::string &name, const T &value) const {

		// The name might contains slashes, so we can navigate through
		// a hierarchy of groups/datasets
		std::vector<std::string> parts = tokenize(name, "/");

		// only the attribute name, read directly and come back
		if( parts.size() == 1 ) {
			return _write_dataset<T>(hdf5_file.openDataSet(name), value);
		}

		// else there's a path to follow, go for it!
		const H5::CommonFG &location = hdf5_file;
		std::vector<std::string> path_parts(parts.begin(), parts.end()-1);
		for(auto const &path: path_parts) {
			// not implemented yet
		}

		throw std::runtime_error("write_dataset still not implemented for attributes in groups/datasets");
	}

	template<typename T>
	void write_dataset_v(const std::string &name, const std::vector<T> &values) const {

		// The name might contains slashes, so we can navigate through
		// a hierarchy of groups/datasets
		std::vector<std::string> parts = tokenize(name, "/");

		// only the attribute name, read directly and come back
		if( parts.size() == 1 ) {
			return _write_dataset_v<T>(hdf5_file.openDataSet(name), values);
		}

		// else there's a path to follow, go for it!
		const H5::CommonFG &location = hdf5_file;
		std::vector<std::string> path_parts(parts.begin(), parts.end()-1);
		for(auto const &path: path_parts) {
			// not implemented yet
		}

		throw std::runtime_error("write_dataset_v still not implemented for attributes in groups/datasets");
	}

private:

	template <typename T>
	void _write_attribute(const H5::H5Location &location, const std::string &name, const T &value) {
		H5::Attribute attr = location.openAttribute(name);
		H5::DataType type = attr.getDataType();
		attr.write(type, &value);
		attr.close();
	}

	template<typename T>
	typename std::enable_if<std::is_arithmetic<T>::value, void>::type
	_write_dataset(const H5::DataSet &dataset, const T &data) {

		H5::DataSpace space = get_1d_dataspace(dataset);
		hsize_t dim_size = get_1d_dimsize(space);
		if ( dim_size != 1 ) {
			std::ostringstream os;
			os << "More than 1 element found in dataset " << dataset.getObjName();
			throw std::runtime_error(os.str());
		}

		dataset.write(&data, dataset.getDataType(), space, space);
	}

	template<typename T>
	typename std::enable_if<std::is_arithmetic<T>::value, void>::type
	_write_dataset_v(const H5::DataSet &dataset, const std::vector<T> &data) const {
		H5::DataSpace space = get_1d_dataspace(dataset);
		hsize_t dim_size = get_1d_dimsize(space);
		dataset.write(data.data(), dataset.getDataType(), space, space);
	}

};

}  // namespace hdf5

}  // namespace shark

#endif // SHARK_HDF5_WRITER_H_