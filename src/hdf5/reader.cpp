//
// Implementation of Reader class methods
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

#include <stdexcept>
#include <string>
#include <vector>

#include "exceptions.h"
#include "logging.h"
#include "hdf5/reader.h"

using namespace std;

namespace shark {

namespace hdf5 {

class attribute_not_found : public std::exception {};
class name_not_found : public std::exception {};


H5::Attribute Reader::get_attribute(const string &name) const {
	LOG(debug) << "Getting attribute " << name << " from file " << get_filename();
	std::vector<std::string> parts = tokenize(name, "/");

	try {
		return _get_attribute(hdf5_file, parts);
	} catch (const attribute_not_found &) {
		std::ostringstream os;
		os << "Attribute " << name << " doesn't exist in " << get_filename();
		throw invalid_data(os.str());
	} catch (const name_not_found &) {
		std::ostringstream os;
		os << "Name " << name << " names no object in file " << get_filename();
		throw invalid_data(os.str());
	}
}

H5::Attribute Reader::_get_attribute(const H5::CommonFG &file_or_group, const std::vector<std::string> &parts) const
{

	// This is the attribute name
	if (parts.size() == 1) {
		// both file and groups derive from H5Location too
		return _get_attribute(dynamic_cast<const H5::H5Location &>(file_or_group), parts[0]);
	}

	auto n_groups = file_or_group.getNumObjs();

	const auto path = parts.front();
	for(hsize_t i = 0; i < n_groups; i++) {

		auto objname = file_or_group.getObjnameByIdx(i);
		if (objname != path) {
			continue;
		}

		auto objtype = file_or_group.getObjTypeByIdx(i);
		if (objtype == H5G_GROUP) {
			std::vector<std::string> subparts(parts.begin() + 1, parts.end());
			return _get_attribute(file_or_group.openGroup(objname), subparts);
		}
		else if (objtype == H5G_DATASET) {
			std::vector<std::string> subparts(parts.begin() + 1, parts.end());
			return _get_attribute(file_or_group.openDataSet(objname), parts.back());
		}
	}

	throw name_not_found();
}

H5::Attribute Reader::_get_attribute(const H5::H5Location &l, const std::string attr_name) const {

	LOG(debug) << "Getting attribute " << attr_name << " from file " << get_filename();

	if (!l.attrExists(attr_name)) {
		throw attribute_not_found();
	}

	return l.openAttribute(attr_name);
}

}  // namespace hdf5

}  // namespace shark