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
 * Header file for the h5reader class
 */

#ifndef SHARK_HDF5_READER
#define SHARK_HDF5_READER

#include <exception>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>

#include "utils.h"
#include "iobase.h"
#include "hdf5/attribute.h"

namespace shark {

namespace hdf5 {

class Reader : public IOBase {

public:

	/**
	 * Constructs a new Reader object.
	 *
	 * @param filename The name of the HDF5 file to read
	 */
	explicit Reader(const std::string& filename) :
			IOBase(filename, FileOpenMethod::Read) {}

	template<typename T>
	T read_attribute(const std::string& name) const {
		Attribute attr = get_attribute(name);
		return attr.read<T>();
	}

	template<typename T>
	T read_dataset(const std::string& name) const {
		return _read_dataset<T>(get_dataset(name));
	}

	template<typename T>
	std::vector<T> read_dataset_v(const std::string& name) const {
		return _read_dataset_v<T>(get_dataset(name));
	}

	template<typename T>
	std::vector<T> read_dataset_v_2(const std::string& name) const {
		return _read_dataset_v_2<T>(get_dataset(name));
	}

private:
	Attribute get_attribute(const std::string& name) const;

	template<typename T>
	typename std::enable_if<std::is_arithmetic<T>::value, T>::type
	_read_dataset(const DataSet& dataset) const {
		DataSpace space = get_scalar_dataspace(dataset);
		T data_out;
		dataset.read(&data_out, dataset.getDataType(), space, space);
		return data_out;
	}

	template<typename T>
	typename std::enable_if<std::is_arithmetic<T>::value, std::vector<T>>::type
	_read_dataset_v(const DataSet& dataset) const {

		DataSpace space = get_1d_dataspace(dataset);
		hsize_t dim_size = get_1d_dimsize(space);

		std::vector<T> data(dim_size);
		dataset.read(data.data(), dataset.getDataType(), space, space);
		return data;
	}

	template<typename T>
	typename std::enable_if<std::is_arithmetic<T>::value, std::vector<T>>::type
	_read_dataset_v_2(const DataSet& dataset) const {

		DataSpace space = get_2d_dataspace(dataset);
		auto dim_sizes = space.getSimpleExtentDims();
		assert(dim_sizes.size() == 2);

		std::vector<T> data(dim_sizes[0] * dim_sizes[1]);
		dataset.read(data.data(), dataset.getDataType(), space, space);
		return data;
	}

};

}  // namespace hdf5

}  // namespace shark

#endif // SHARK_HDF5_READER