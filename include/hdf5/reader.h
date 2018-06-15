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
#include "hdf5/iobase.h"

namespace shark {

namespace hdf5 {

class Reader : public IOBase {

public:

	/**
	 * Constructs a new Reader object.
	 *
	 * @param filename The name of the HDF5 file to read
	 */
	Reader(const std::string &filename) :
		IOBase(filename, H5F_ACC_RDONLY) {}

	template<typename T>
	const T read_attribute(const std::string &name) const {
		std::string attr_name;
		H5::Attribute attr = get_attribute(name);
		H5::DataType type = attr.getDataType();
		T val;
		attr.read(type, &val);
		attr.close();
		return val;
	}

	template<typename T>
	T read_dataset(const std::string &name) const {
		return _read_dataset<T>(get_dataset(name));
	}

	template<typename T>
	std::vector<T> read_dataset_v(const std::string &name) const {
		return _read_dataset_v<T>(get_dataset(name));
	}

	template<typename T>
	std::vector<T> read_dataset_v_2(const std::string &name) const {
		return _read_dataset_v_2<T>(get_dataset(name));
	}

private:

	H5::Attribute get_attribute(const std::string &name) const;

	template<typename T>
	typename std::enable_if<std::is_arithmetic<T>::value, T>::type
	_read_dataset(const H5::DataSet &dataset) const {
		H5::DataSpace space = get_scalar_dataspace(dataset);
		T data_out;
		dataset.read(&data_out, dataset.getDataType(), space, space);
		return data_out;
	}

	template<typename T>
	typename std::enable_if<std::is_arithmetic<T>::value, std::vector<T>>::type
	_read_dataset_v(const H5::DataSet &dataset) const {

		H5::DataSpace space = get_1d_dataspace(dataset);
		hsize_t dim_size = get_1d_dimsize(space);

		std::vector<T> data(dim_size);
		dataset.read(data.data(), dataset.getDataType(), space, space);
		return data;
	}

	template<typename T>
	typename std::enable_if<std::is_arithmetic<T>::value, std::vector<T>>::type
	_read_dataset_v_2(const H5::DataSet &dataset) const {

		H5::DataSpace space = get_2d_dataspace(dataset);
		hsize_t dim_sizes[2];
		space.getSimpleExtentDims(dim_sizes, NULL);

		std::vector<T> data(dim_sizes[0] * dim_sizes[1]);
		dataset.read(data.data(), dataset.getDataType(), space, space);
		return data;
	}

};

}  // namespace hdf5

}  // namespace shark

#endif // SHARK_HDF5_READER