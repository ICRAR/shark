//
// Implementation of Writer class methods
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

#include <utility>

#include "hdf5/writer.h"
#include "exceptions.h"
#include "logging.h"


namespace shark {

namespace hdf5 {


template <H5G_obj_t E>
inline static
typename entity_traits<E>::rettype
get_entity(const H5::CommonFG &file_or_group, const std::string &name);

template <> inline
H5::Group get_entity<H5G_GROUP>(const H5::CommonFG &file_or_group, const std::string &name)
{
	return file_or_group.openGroup(name);
}

template <> inline
H5::DataSet get_entity<H5G_DATASET>(const H5::CommonFG &file_or_group, const std::string &name)
{
	return file_or_group.openDataSet(name);
}

template <H5G_obj_t E, typename ... Ts>
inline static
typename entity_traits<E>::rettype
create_entity(const H5::CommonFG &file_or_group, const std::string &name, Ts&&...create_args);

template <> inline
H5::Group create_entity<H5G_GROUP>(const H5::CommonFG &file_or_group, const std::string &name)
{
	return file_or_group.createGroup(name);
}

template <> inline
H5::DataSet create_entity<H5G_DATASET>(const H5::CommonFG &file_or_group, const std::string &name, const H5::DataType &dataType, const H5::DataSpace &dataSpace)
{
	return file_or_group.createDataSet(name, dataType, dataSpace);
}

template <H5G_obj_t E, typename ... Ts>
typename entity_traits<E>::rettype
ensure_entity(const H5::CommonFG &file_or_group, const std::string &name, Ts&&...create_args)
{
	// Group exists?
	bool exists = false;
	for(hsize_t i = 0; i < file_or_group.getNumObjs(); i++) {

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

H5::Group Writer::ensure_group(const std::vector<std::string> &path) const
{
	if (path.size() == 1) {
		return ensure_entity<H5G_GROUP>(hdf5_file, path[0]);
	}

	H5::Group group =  ensure_entity<H5G_GROUP>(hdf5_file, path.front());
	std::vector<std::string> group_paths(path.begin() + 1, path.end());
	for(auto part: path) {
		group = ensure_entity<H5G_GROUP>(group, part);
	}
	return group;
}

H5::DataSet Writer::ensure_dataset(const std::vector<std::string> &path, const H5::DataType &dataType, const H5::DataSpace &dataSpace) const
{
	if (path.size() == 1) {
		return ensure_entity<H5G_DATASET>(hdf5_file, path[0], dataType, dataSpace);
	}

	std::vector<std::string> group_paths(path.begin(), path.end() - 1);
	H5::Group group = ensure_group(group_paths);
	return ensure_entity<H5G_DATASET>(group, path.back(), dataType, dataSpace);
}


}  // namespace hdf5

}  // namespace shark