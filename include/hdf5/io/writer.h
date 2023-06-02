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

#include "hdf5/data_type.h"
#include "hdf5/data_set.h"
#include "hdf5/attribute.h"
#include "iobase.h"
#include "traits.h"
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
	explicit object_exists(const std::string& what) : exception(what) {}
};

// How we construct data types depends on the type
template<typename T>
static inline
DataType _datatype(const T& val) {
	return DataType(datatype_traits<T>::write_type);
}

template<>
inline
DataType _datatype<std::string>(const std::string& val) {
	// Allow space for null terminator
	return StringDataType(val.size() + 1);
}

// Overwriting for vectors of values
template<typename T>
static inline
DataType _datatype(const std::vector<T>& val) {
	return DataType(datatype_traits<T>::write_type);
}

// Overwriting for vectors of values
template<typename T>
static inline
DataType _datatype(const std::vector<std::vector<T>>& val) {
	return DataType(datatype_traits<T>::write_type);
}

template<>
inline
DataType _datatype<std::string>(const std::vector<std::string>& val) {
	return StringDataType();
}

// How we write datasets, depending on the data type
template<typename T>
static inline
void _write_dataset(DataSet& dataset, const DataType& dataType, const DataSpace& dataSpace, const T& val) {
	DataType mem_dataType(datatype_traits<T>::native_type);
	dataset.write(&val, mem_dataType, dataSpace, dataSpace);
}

// Specialization for bools, which we need to convert into int first
template<>
inline
void _write_dataset<bool>(DataSet& dataset, const DataType& dataType, const DataSpace& dataSpace, const bool& val) {
	int int_val = static_cast<int>(val);
	DataType mem_dataType(PredefinedDataType::NATIVE_INT());
	dataset.write(&int_val, mem_dataType, dataSpace, dataSpace);
}

// Specialization for strings, for which there is a specific .write method
template<>
inline
void _write_dataset<std::string>(DataSet& dataset, const DataType& dataType, const DataSpace& dataSpace,
                                 const std::string& val) {
	dataset.write(val.data(), dataType, dataSpace, dataSpace);
}

// Overwriting of _write_datasets for vectors
template<typename T>
static inline
void
_write_dataset(DataSet& dataset, const DataType& dataType, const DataSpace& dataSpace, const std::vector<T>& vals) {
	DataType mem_dataType(datatype_traits<T>::native_type);
	dataset.write(vals.data(), mem_dataType, dataSpace, dataSpace);
}

// and specializing for vectors of strings
template<>
inline
void _write_dataset<std::string>(DataSet& dataset, const DataType& dataType, const DataSpace& dataSpace,
                                 const std::vector<std::string>& vals) {
	std::vector<const char*> c_strings;
	std::transform(vals.begin(), vals.end(), std::back_inserter(c_strings), [](const std::string& s) {
		return s.c_str();
	});
	dataset.write(c_strings.data(), dataType, dataSpace, dataSpace);
}

// Overwriting of _write_datasets for vectors of vectors
template<typename T>
static inline
void _write_dataset(DataSet& dataset, const DataType& dataType, DataSpace& fDataSpace,
                    const std::vector<std::vector<T>>& vals) {
	DataType mem_dataType(datatype_traits<T>::native_type);

	// Find out maximum dimensions of the dataset
	auto dataset_max_dims = fDataSpace.getSimpleExtentMaxDims();

	// We iterate over each inner vector and write it as a new row
	// in the dataset. For this, we set the parameters for selecting where the
	// data will be written into the file. The fstart will thus keep changing
	// to write each row vector consecutively
	std::vector<hsize_t> fcount = {1, dataset_max_dims[1]};
	std::vector<hsize_t> fstart = {0, 0};
	std::vector<hsize_t> fstride = {1, 1};
	std::vector<hsize_t> fblock = {1, 1};

	// The data from the vector should only be of length dataset_max_dims[1]
	auto mDataSpace = DataSpace::create({dataset_max_dims[1]});
	for (auto& row: vals) {
		fDataSpace.selectHyperslab(HyperslabSelection::Set, fstart, fstride, fcount, fblock);
		dataset.write(row.data(), mem_dataType, mDataSpace, fDataSpace);
		fstart[0]++;
	}
}

template<typename AttributeHolder, typename T>
static inline
void _create_and_write_attribute(AttributeHolder& dataset, const std::string& name, const T& value) {
	DataType dataType = _datatype<T>(value);
	auto attr = dataset.createAttribute(name, dataType, DataSpace::create(DataSpaceType::Scalar));
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
	explicit Writer(const std::string& filename, bool overwrite = true,
	                naming_convention group_naming_convention = naming_convention::SNAKE_CASE,
	                naming_convention dataset_naming_convention = naming_convention::SNAKE_CASE,
	                naming_convention attr_naming_convention = naming_convention::SNAKE_CASE);

	void set_comment(DataSet& dataset, const std::string& comment) {
		if (comment.empty()) {
			return;
		}

		dataset.setComment(comment);

		// Follow naming convention, "comment" works with snake_case and lowerCamelCase
		auto comment_attr_name = "comment";
		if (attr_naming_convention == naming_convention::CAMEL_CASE) {
			comment_attr_name = "Comment";
		}
		_create_and_write_attribute(dataset, comment_attr_name, comment);
	}

	template<typename T>
	void write_attribute(const std::string& name, const T& value) {

		std::vector<std::string> parts = tokenize(name, "/");

		// Cannot attach attributes directly to files
		if (parts.size() == 1) {
			std::ostringstream os;
			os << "cannot attach attribute " << name << " directly to HDF5 file";
			throw invalid_argument(os.str());
		}

		std::vector<std::string> path(parts.begin(), parts.end() - 1);
		auto& attr_name = parts.back();
		check_attr_name(attr_name);

		// Get the corresponding group/dataset and write the attribute there
		try {
			auto group = get_or_create_group(path);
			_create_and_write_attribute(group, attr_name, value);
		} catch (const object_exists&) {
			// name doesn't point to a group but to an existing dataset
			// if this one fails we give up
			auto dataset = get_dataset(path);
			_create_and_write_attribute(dataset, attr_name, value);
		}

	}

	template<typename T>
	void write_dataset(const std::string& name, const T& value, const std::string& comment = NO_COMMENT) {
		auto dataSpace = DataSpace::create(DataSpaceType::Scalar);
		DataType dataType = _datatype<T>(value);
		auto dataset = get_or_create_dataset(tokenize(name, "/"), dataType, dataSpace);
		set_comment(dataset, comment);
		_write_dataset(dataset, dataType, dataSpace, value);
	}

	template<typename T>
	void write_dataset(const std::string& name, const std::vector<T>& values, const std::string& comment = NO_COMMENT) {
		auto dataSpace = DataSpace::create({values.size()});
		DataType dataType = _datatype<T>(values);
		auto dataset = get_or_create_dataset(tokenize(name, "/"), dataType, dataSpace);
		set_comment(dataset, comment);
		_write_dataset(dataset, dataType, dataSpace, values);
	}

	template<typename T>
	void write_dataset(const std::string& name, const std::vector<std::vector<T>>& values,
	                   const std::string& comment = NO_COMMENT) {
		if (values.empty()) {
			return;
		}
		auto dataSpace = DataSpace::create({values.size(), values[0].size()});
		DataType dataType = _datatype<T>(values);
		auto dataset = get_or_create_dataset(tokenize(name, "/"), dataType, dataSpace);
		set_comment(dataset, comment);
		_write_dataset(dataset, dataType, dataSpace, values);
	}

private:

	Group get_or_create_group(const std::vector<std::string>& path);

	DataSet
	get_or_create_dataset(const std::vector<std::string>& path, const DataType& dataType, const DataSpace& dataSpace);

	naming_convention group_naming_convention;
	naming_convention dataset_naming_convention;
	naming_convention attr_naming_convention;

	void check_group_name(const std::string& group_name) const;

	void check_dataset_name(const std::string& dataset_name) const;

	void check_attr_name(const std::string& attr_name) const;

	static const std::string NO_COMMENT;
};

}  // namespace hdf5

}  // namespace shark

#endif // SHARK_HDF5_WRITER_H_