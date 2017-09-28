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

#include <iterator>
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

// How we construct data types depends on the type
template <typename T>
static inline
H5::DataType _datatype(const T &val)
{
	return H5::DataType(datatype_traits<T>::write_type);
}

template <>
inline
H5::DataType _datatype<std::string>(const std::string &val)
{
	return H5::StrType(H5::PredType::C_S1, val.size());
}

// Overwriting for vectors of values
template <typename T>
static inline
H5::DataType _datatype(const std::vector<T> &val)
{
	return H5::DataType(datatype_traits<T>::write_type);
}

template <>
inline
H5::DataType _datatype<std::string>(const std::vector<std::string> &val)
{
	return H5::StrType(H5::PredType::C_S1, H5T_VARIABLE);
}

// How we write attributes, depending on the data type
// Strings need some special treatment, so we specialise for them
template <typename T>
static inline
void _write_attribute(const H5::Attribute &attr, const H5::DataType &dataType, const T &val)
{
	attr.write(dataType, &val);
}

template <>
inline
void _write_attribute<std::string>(const H5::Attribute &attr, const H5::DataType &dataType, const std::string &val)
{
	attr.write(dataType, val);
}

// How we write datasets, depending on the data type
template <typename T>
static inline
void _write_dataset(const H5::DataSet &dataset, const H5::DataType &dataType, const H5::DataSpace &dataSpace, const T &val)
{
	H5::DataType mem_dataType(datatype_traits<T>::native_type);
	dataset.write(&val, mem_dataType, dataSpace, dataSpace);
}

// Specialization for bools, which we need to convert into int first
template <>
inline
void _write_dataset<bool>(const H5::DataSet &dataset, const H5::DataType &dataType, const H5::DataSpace &dataSpace, const bool &val)
{
	int int_val = static_cast<int>(val);
	H5::DataType mem_dataType(H5::PredType::NATIVE_INT);
	dataset.write(&int_val, mem_dataType, dataSpace, dataSpace);
}

// Specialization for strings, for which there is a specific .write method
template <>
inline
void _write_dataset<std::string>(const H5::DataSet &dataset, const H5::DataType &dataType, const H5::DataSpace &dataSpace, const std::string &val)
{
	dataset.write(val, dataType, dataSpace, dataSpace);
}

// Overwriting of _write_datasets for vectors
template <typename T>
static inline
void _write_dataset(const H5::DataSet &dataset, const H5::DataType &dataType, const H5::DataSpace &dataSpace, const std::vector<T> &vals)
{
	H5::DataType mem_dataType(datatype_traits<T>::native_type);
	dataset.write(vals.data(), mem_dataType, dataSpace, dataSpace);
}

// and specializing for vectors of strings
template <>
inline
void _write_dataset<std::string>(const H5::DataSet &dataset, const H5::DataType &dataType, const H5::DataSpace &dataSpace, const std::vector<std::string> &vals)
{
	std::vector<const char *> c_strings;
	std::transform(vals.begin(), vals.end(), std::back_inserter(c_strings), [](const std::string &s) {
		return s.c_str();
	});
	dataset.write(c_strings.data(), dataType, dataSpace, dataSpace);
}

template<typename T>
static inline
void _create_and_write_attribute(H5::H5Location &loc, const std::string &name, const T &value) {
	H5::DataType dataType = _datatype<T>(value);
	auto attr = loc.createAttribute(name, dataType, H5::DataSpace(H5S_SCALAR));
	_write_attribute<T>(attr, dataType, value);
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
	Writer(const std::string &filename, bool overwrite = true);

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
			// name doesn't point to a group but to an existing dataset
			// if this one fails we give up
			auto dataset = get_dataset(path);
			_create_and_write_attribute(dataset, parts.back(), value);
		}

	}

	template<typename T>
	void write_dataset(const std::string &name, const T &value) {
		H5::DataSpace dataSpace(H5S_SCALAR);
		H5::DataType dataType = _datatype<T>(value);
		auto dataset = ensure_dataset(tokenize(name, "/"), dataType, dataSpace);
		_write_dataset(dataset, dataType, dataSpace, value);
	}

	template<typename T>
	void write_dataset(const std::string &name, const std::vector<T> &values) {
		const hsize_t size = values.size();
		H5::DataSpace dataSpace(1, &size);
		H5::DataType dataType = _datatype<T>(values);
		auto dataset = ensure_dataset(tokenize(name, "/"), dataType, dataSpace);
		_write_dataset(dataset, dataType, dataSpace, values);
	}

private:

	H5::Group ensure_group(const std::vector<std::string> &path) const;
	H5::DataSet ensure_dataset(const std::vector<std::string> &path, const H5::DataType &dataType, const H5::DataSpace &dataSpace) const;

};

}  // namespace hdf5

}  // namespace shark

#endif // SHARK_HDF5_WRITER_H_