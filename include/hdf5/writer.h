//
// Header file for the hdf5::Writer class
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

#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <H5Cpp.h>

#include "hdf5/iobase.h"
#include "hdf5/traits.h"
#include "exceptions.h"
#include "utils.h"

namespace shark {

namespace hdf5 {

/**
 * Exception thrown when a group or a dataset path is given by the user,
 * but the path refers to an existing, different entity.
 *
 */
class object_exists : public exception {
public:
	object_exists(const std::string &what) : exception(what) {}
};


// How we write attributes, depending on the data type
// Strings need some special treatment, so we specialise for them
template<typename T>
static inline
void _create_and_write_attribute(H5::H5Location &loc, const std::string &name, const T &value) {
	H5::DataType dataType(datatype_traits<T>::write_type);
	auto attr = loc.createAttribute(name, dataType, H5::DataSpace(H5S_SCALAR));
	attr.write(H5::DataType(datatype_traits<T>::native_type), &value);
}

template<>
inline
void _create_and_write_attribute<std::string>(H5::H5Location &loc, const std::string &name, const std::string &value) {
	H5::StrType dataType(H5::PredType::C_S1, value.size());
	hsize_t size = value.size();
	auto attr = loc.createAttribute(name, dataType, H5::DataSpace(H5S_SCALAR));
	attr.write(dataType, value);
}


/**
 * An object that can write data in the form of attributes and datasets into an
 * HDF5 file.
 */
class Writer : public IOBase {

public:

	/**
	 * Constructs a new Writer object.
	 *
	 * @param filename The name of the HDF5 file to write
	 * @param overwrite Whether existing files should be overwritten or not
	 */
	Writer(const std::string &filename, bool overwrite = true) :
		IOBase(filename, H5F_ACC_RDWR | (overwrite ? H5F_ACC_CREAT : H5F_ACC_EXCL)) {}

	template<typename T>
	void write_attribute(const std::string &name, const T &value) {

		std::vector<std::string> parts = tokenize(name, "/");

		// only the attribute name, read directly and come back
		if( parts.size() == 1 ) {
			_create_and_write_attribute<T>(hdf5_file, name, value);
			return;
		}

		// Get the corresponding group/dataset and write the attribute there
		std::vector<std::string> path(parts.begin(), parts.end()-1);
		try {
			auto group = ensure_group(path);
			_create_and_write_attribute<T>(group, parts.back(), value);
		} catch (const object_exists &e) {
			throw;
//			auto dataset = open_dataset(path);
//			_write_attribute<T>(dataset, parts.back(), value);
		}

	}

	template<typename T>
	void write_dataset(const std::string &name, const T &value) {
		std::vector<std::string> parts = tokenize(name, "/");
		H5::DataSpace dataSpace(H5S_SCALAR);
		H5::StrType dataType(H5::PredType::NATIVE_CHAR);
		auto dataset = ensure_dataset(parts, dataType, dataSpace);
		_write_dataset(dataset, value);
	}

	template<typename T>
	void write_dataset_v(const std::string &name, const std::vector<T> &values) {
		std::vector<std::string> parts = tokenize(name, "/");
		const hsize_t size = values.size();
		H5::DataSpace dataSpace(1, &size);
		H5::StrType dataType(H5::PredType::NATIVE_CHAR);
		auto dataset = ensure_dataset(parts, dataType, dataSpace);
		_write_dataset_v(dataset, values);
	}

private:

	template<typename T>
	void _write_dataset(const H5::DataSet &dataset, const T &data) {

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
	void _write_dataset_v(const H5::DataSet &dataset, const std::vector<T> &data) {
		H5::DataSpace space = get_1d_dataspace(dataset);
		hsize_t dim_size = get_1d_dimsize(space);
		dataset.write(data.data(), dataset.getDataType(), space, space);
	}

	H5::Group ensure_group(const std::vector<std::string> &path) const;
	H5::DataSet ensure_dataset(const std::vector<std::string> &path, const H5::DataType &dataType, const H5::DataSpace &dataSpace) const;

};

}  // namespace hdf5

}  // namespace shark

#endif // SHARK_HDF5_WRITER_H_