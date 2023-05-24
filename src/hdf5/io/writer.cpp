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
 * Implementation of Writer class methods
 */

#include <sstream>
#include <utility>

#include "hdf5/io/writer.h"
#include "exceptions.h"


namespace shark {
namespace hdf5 {

const std::string Writer::NO_COMMENT;

Writer::Writer(const std::string& filename, bool overwrite,
               naming_convention group_naming_convention, naming_convention dataset_naming_convention,
               naming_convention attr_naming_convention) :
		IOBase(filename, overwrite ? FileOpenMethod::Truncate : FileOpenMethod::ExclusiveWrite),
		group_naming_convention(group_naming_convention),
		dataset_naming_convention(dataset_naming_convention),
		attr_naming_convention(attr_naming_convention) {
}

static
void _check_entity_name(const std::string& name, const char* entity_type, naming_convention convention) {
	if (!follows_convention(name, convention)) {
		std::ostringstream os;
		os << entity_type << " name " << name << " does not follow the " << convention << " naming convention";
		throw invalid_argument(os.str());
	}
}

void Writer::check_group_name(const std::string& group_name) const {
	_check_entity_name(group_name, "Group", group_naming_convention);
}

void Writer::check_dataset_name(const std::string& dataset_name) const {
	_check_entity_name(dataset_name, "Dataset", dataset_naming_convention);
}

void Writer::check_attr_name(const std::string& attr_name) const {
	_check_entity_name(attr_name, "Attribute", attr_naming_convention);
}

template<H5G_obj_t E>
inline static
typename entity_traits<E>::rettype
get_entity(const AbstractGroup& file_or_group, const std::string& name);

template<>
inline
Group get_entity<H5G_GROUP>(const AbstractGroup& file_or_group, const std::string& name) {
	return file_or_group.openGroup(name);
}

template<>
inline
DataSet get_entity<H5G_DATASET>(const AbstractGroup& file_or_group, const std::string& name) {
	return file_or_group.openDataSet(name);
}

template<H5G_obj_t E, typename ... Ts>
inline static
typename entity_traits<E>::rettype
create_entity(AbstractGroup& file_or_group, const std::string& name, Ts&& ...create_args);

template<>
inline
Group create_entity<H5G_GROUP>(AbstractGroup& file_or_group, const std::string& name) {
	return file_or_group.createGroup(name);
}

template<>
inline
DataSet create_entity<H5G_DATASET>(AbstractGroup& file_or_group, const std::string& name, const DataType& dataType,
                                   const DataSpace& dataSpace) {
	return file_or_group.createDataSet(name, dataType, dataSpace);
}

template<H5G_obj_t E, typename ... Ts>
typename entity_traits<E>::rettype
get_or_create_entity(AbstractGroup& file_or_group, const std::string& name, Ts&& ...create_args) {
	// Loop through subobjects and find entity
	for (hsize_t i = 0; i < file_or_group.getNumObjs(); i++) {

		auto objtype = file_or_group.getObjTypeByIdx(i);
		auto objname = file_or_group.getObjnameByIdx(i);

		// Name already used, sorry!
		if (name == objname && objtype != E) {
			std::ostringstream os;
			os << "Name " << name << " is already used by an object";
			throw object_exists(os.str());
		}

			// entity exists!
		else if (name == objname && objtype == E) {
			return get_entity<E>(file_or_group, name);
		}
	}

	// Nothing found, create new entity
	return create_entity<E>(file_or_group, name, std::forward<Ts>(create_args)...);
}

Group Writer::get_or_create_group(const std::vector<std::string>& path) {
	if (path.size() == 1) {
		check_group_name(path[0]);
		return get_or_create_entity<H5G_GROUP>(hdf5_file.value(), path[0]);
	}

	Group group = get_or_create_entity<H5G_GROUP>(hdf5_file.value(), path.front());
	std::vector<std::string> group_paths(path.begin() + 1, path.end());
	for (auto& part: group_paths) {
		group = get_or_create_entity<H5G_GROUP>(group, part);
		check_group_name(part);
	}
	return group;
}

DataSet Writer::get_or_create_dataset(const std::vector<std::string>& path, const DataType& dataType,
                                      const DataSpace& dataSpace) {
	if (path.size() == 1) {
		check_dataset_name(path[0]);
		return get_or_create_entity<H5G_DATASET>(hdf5_file.value(), path[0], dataType, dataSpace);
	}

	std::vector<std::string> group_paths(path.begin(), path.end() - 1);
	auto& dataset_name = path.back();
	check_dataset_name(dataset_name);
	Group group = get_or_create_group(group_paths);
	return get_or_create_entity<H5G_DATASET>(group, dataset_name, dataType, dataSpace);
}

}  // namespace hdf5

}  // namespace shark