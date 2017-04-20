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

#ifndef SHARK_HDF5_READER
#define SHARK_HDF5_READER

#include <exception>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>

#include <H5Cpp.h>

#include "utils.h"

namespace shark {

namespace hdf5 {

class Reader {

public:

	/**
	 * Constructs a new Reader object.
	 *
	 * @param filename The name of the HDF5 file to read
	 */
	Reader(const std::string &filename);

	/**
	 * Destructor. It closes the underlying file.
	 */
	~Reader();

	/**
	 * Returns the name of the file being read
	 * @return The name of the file being read
	 */
	const std::string get_filename() const;

	template<typename T>
	const T read_attribute(const std::string &name) const {

		// The name might contains slashes, so we can navigate through
		// a hierarchy of groups/datasets
		std::vector<std::string> parts = tokenize(name, "/");

		// only the attribute name, read directly and come back
		if( parts.size() == 1 ) {
			return _read_attribute<T>(hdf5_file.get(), name);
		}

		// else there's a path to follow, go for it!
		H5::CommonFG *location = hdf5_file.get();
		std::vector<std::string> path_parts(parts.begin(), parts.end()-1);
		for(auto const &path: path_parts) {
			// not implemented yet
		}

		throw std::runtime_error("read_attribute still not implemented for attributes in groups/datasets");
	}

	template<typename T>
	T read_dataset(const std::string &name) const {

		// The name might contains slashes, so we can navigate through
		// a hierarchy of groups/datasets
		std::vector<std::string> parts = tokenize(name, "/");

		// only the attribute name, read directly and come back
		if( parts.size() == 1 ) {
			return _read_dataset<T>(hdf5_file->openDataSet(name));
		}

		// else there's a path to follow, go for it!
		H5::CommonFG *location = hdf5_file.get();
		std::vector<std::string> path_parts(parts.begin(), parts.end()-1);
		for(auto const &path: path_parts) {
			// not implemented yet
		}

		throw std::runtime_error("read_attribute still not implemented for attributes in groups/datasets");
	}

	template<typename T>
	std::vector<T> read_dataset_v(const std::string &name) const {

		// The name might contains slashes, so we can navigate through
		// a hierarchy of groups/datasets
		std::vector<std::string> parts = tokenize(name, "/");

		// only the attribute name, read directly and come back
		if( parts.size() == 1 ) {
			return _read_dataset_v<T>(hdf5_file->openDataSet(name));
		}

		// else there's a path to follow, go for it!
		H5::CommonFG *location = static_cast<H5::CommonFG *>(hdf5_file.get());
		std::vector<std::string> path_parts(parts.begin(), parts.end()-1);
		for(auto const &path: path_parts) {
			// not implemented yet
		}

		throw std::runtime_error("read_attribute still not implemented for attributes in groups/datasets");
	}

private:
	std::unique_ptr<H5::H5File> hdf5_file;

	template <typename T>
	const T _read_attribute(H5::H5Location *location, const std::string &name) const {
		H5::Attribute attr = location->openAttribute(name);
		H5::DataType type = attr.getDataType();
		T val;
		attr.read(type, &val);
		attr.close();
		return val;
	}

	template<typename T>
	typename std::enable_if<std::is_arithmetic<T>::value, T>::type
	_read_dataset(const H5::DataSet &dataset) const {

		H5::DataSpace space = _get_1d_dataspace(dataset);
		hsize_t dim_size = _get_1d_dimsize(space);
		if ( dim_size != 1 ) {
			std::ostringstream os;
			os << "More than 1 element found in dataset " << dataset.getObjName();
			throw std::runtime_error(os.str());
		}

		T data_out;
		dataset.read(&data_out, dataset.getDataType(), space, space);
		return data_out;
	}

	template<typename T>
	typename std::enable_if<std::is_arithmetic<T>::value, std::vector<T>>::type
	_read_dataset_v(const H5::DataSet &dataset) const {

		H5::DataSpace space = _get_1d_dataspace(dataset);
		hsize_t dim_size = _get_1d_dimsize(space);

		std::vector<T> data(dim_size);
		dataset.read(data.data(), dataset.getDataType(), space, space);
		return data;
	}

	H5::DataSpace _get_1d_dataspace(const H5::DataSet &dataset) const {
		H5::DataSpace space = dataset.getSpace();
		int ndims = space.getSimpleExtentNdims();
		if ( ndims != 1 ) {
			std::ostringstream os;
			os << "More than one dimension found in dataset " << dataset.getObjName();
			throw std::runtime_error(os.str());
		}
		return space;
	}

	hsize_t _get_1d_dimsize(const H5::DataSpace &space) const {
		hsize_t dim_size;
		space.getSimpleExtentDims(&dim_size, NULL);
		return dim_size;
	}
};

}  // namespace hdf5

}  // namespace shark

#endif // SHARK_HDF5_READER