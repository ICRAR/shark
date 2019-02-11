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
 * Header file for the hdf5::Writer class
 */

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
#include "naming_convention.h"
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

// Overwriting for vectors of values
template <typename T>
static inline
H5::DataType _datatype(const std::vector<std::vector<T>> &val)
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

// Overwriting of _write_datasets for vectors of vectors
template <typename T>
static inline
void _write_dataset(const H5::DataSet &dataset, const H5::DataType &dataType, const H5::DataSpace &fDataSpace, const std::vector<std::vector<T>> &vals)
{

	H5::DataType mem_dataType(datatype_traits<T>::native_type);

	// Find out maximum dimensions of the dataset
	hsize_t dataset_max_dims[2], dataset_dims[2];
	fDataSpace.getSimpleExtentDims(dataset_dims, dataset_max_dims);

	// We iterate over each inner vector and write it as a new row
	// in the dataset. For this, we set the parameters for selecting where the
	// data will be written into the file. The fstart will thus keep changing
	// to write each row vector consecutively
	hsize_t fcount[2] = {1, dataset_max_dims[1]};
	hsize_t fstart[2] = {0, 0};
	hsize_t fstride[2] = {1, 1};
	hsize_t fblock[2] = {1, 1};

	// The data from the vector should only be of length dataset_max_dims[1]
	H5::DataSpace mDataSpace(1, dataset_max_dims + 1);
	for(auto &row: vals) {
		fDataSpace.selectHyperslab(H5S_SELECT_SET, fcount, fstart, fstride, fblock);
		dataset.write(row.data(), mem_dataType, mDataSpace, fDataSpace);
		fstart[0]++;
	}
}

template<typename AttributeHolder, typename T>
static inline
void _create_and_write_attribute(AttributeHolder &dataset, const std::string &name, const T &value) {
	H5::DataType dataType = _datatype<T>(value);
	auto attr = dataset.createAttribute(name, dataType, H5::DataSpace(H5S_SCALAR));
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
	Writer(const std::string &filename, bool overwrite = true,
		naming_convention group_naming_convention = naming_convention::SNAKE_CASE,
		naming_convention dataset_naming_convention = naming_convention::SNAKE_CASE,
		naming_convention attr_naming_convention = naming_convention::SNAKE_CASE);

	void set_comment(H5::DataSet &dataset, const std::string &comment)
	{
		if (comment.empty()) {
			return;
		}
#ifdef HDF5_NEWER_THAN_1_8_11
		dataset.setComment(comment);
#endif

		// Follow naming convention, "comment" works with snake_case and lowerCamelCase
		auto comment_attr_name = "comment";
		if (attr_naming_convention == naming_convention::CAMEL_CASE) {
			comment_attr_name = "Comment";
		}
		_create_and_write_attribute(dataset, comment_attr_name, comment);
	}

	template<typename T>
	void write_attribute(const std::string &name, const T &value) {

		std::vector<std::string> parts = tokenize(name, "/");

		// Cannot attach attributes directly to files
		if (parts.size() == 1) {
			std::ostringstream os;
			os << "cannot attach attribute " << name << " directly to HDF5 file";
			throw invalid_argument(os.str());
		}

		std::vector<std::string> path(parts.begin(), parts.end()-1);
		auto &attr_name = parts.back();
		check_attr_name(attr_name);

		// Get the corresponding group/dataset and write the attribute there
		try {
			auto group = ensure_group(path);
			_create_and_write_attribute(group, attr_name, value);
		} catch (const object_exists &) {
			// name doesn't point to a group but to an existing dataset
			// if this one fails we give up
			auto dataset = get_dataset(path);
			_create_and_write_attribute(dataset, attr_name, value);
		}

	}

	template<typename T>
	void write_dataset(const std::string &name, const T &value, const std::string &comment = NO_COMMENT) {
		H5::DataSpace dataSpace(H5S_SCALAR);
		H5::DataType dataType = _datatype<T>(value);
		auto dataset = ensure_dataset(tokenize(name, "/"), dataType, dataSpace);
		set_comment(dataset, comment);
		_write_dataset(dataset, dataType, dataSpace, value);
	}

	template<typename T>
	void write_dataset(const std::string &name, const std::vector<T> &values, const std::string &comment = NO_COMMENT) {
		const hsize_t size = values.size();
		H5::DataSpace dataSpace(1, &size);
		H5::DataType dataType = _datatype<T>(values);
		auto dataset = ensure_dataset(tokenize(name, "/"), dataType, dataSpace);
		set_comment(dataset, comment);
		_write_dataset(dataset, dataType, dataSpace, values);
	}

	template<typename T>
	void write_dataset(const std::string &name, const std::vector<std::vector<T>> &values, const std::string &comment = NO_COMMENT) {
		if (values.empty()) {
			return;
		}
		const hsize_t sizes[] = {values.size(), values[0].size()};
		H5::DataSpace dataSpace(2, sizes);
		H5::DataType dataType = _datatype<T>(values);
		auto dataset = ensure_dataset(tokenize(name, "/"), dataType, dataSpace);
		set_comment(dataset, comment);
		_write_dataset(dataset, dataType, dataSpace, values);
	}

private:

	H5::Group ensure_group(const std::vector<std::string> &path) const;
	H5::DataSet ensure_dataset(const std::vector<std::string> &path, const H5::DataType &dataType, const H5::DataSpace &dataSpace) const;

	naming_convention group_naming_convention;
	naming_convention dataset_naming_convention;
	naming_convention attr_naming_convention;

	void check_group_name(const std::string &group_name) const;
	void check_dataset_name(const std::string &dataset_name) const;
	void check_attr_name(const std::string &attr_name) const;

	static const std::string NO_COMMENT;
};

}  // namespace hdf5

}  // namespace shark

#endif // SHARK_HDF5_WRITER_H_