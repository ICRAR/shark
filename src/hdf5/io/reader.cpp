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

/**
 * @file
 *
 * Implementation of Reader class methods
 */

#include <stdexcept>
#include <string>
#include <vector>

#include "exceptions.h"
#include "logging.h"
#include "hdf5/io/reader.h"

namespace shark {
namespace hdf5 {

class attribute_not_found : public std::exception {
};

class name_not_found : public std::exception {
};

static Attribute _get_attribute(const Location& l, const std::string& attr_name) {
	if (!l.attributeExists(attr_name)) {
		throw attribute_not_found();
	}
	return l.openAttribute(attr_name); // l.openAttribute(attr_name);
}

static Attribute _get_attribute(const AbstractGroup& file_or_group, const std::vector<std::string>& parts) {
	// This is the attribute name
	if (parts.size() == 1) {
		// This is a group (we don't support attributes in files)
		// This throws cast exception, do we want a custom error for this?
		const auto& attribute_holder = dynamic_cast<const Location&>(file_or_group);
		return _get_attribute(attribute_holder, parts[0]);
	}

	auto n_groups = file_or_group.getNumObjs();

	auto& path = parts.front();
	for (hsize_t i = 0; i < n_groups; i++) {

		auto objname = file_or_group.getObjnameByIdx(i);
		if (objname != path) {
			continue;
		}

		auto objtype = file_or_group.getObjTypeByIdx(i);
		if (objtype == H5G_GROUP) {
			std::vector<std::string> subparts(parts.begin() + 1, parts.end());
			return _get_attribute(file_or_group.openGroup(objname), subparts);
		} else if (objtype == H5G_DATASET) {
			std::vector<std::string> subparts(parts.begin() + 1, parts.end());
			return _get_attribute(file_or_group.openDataSet(objname), parts.back());
		}
	}

	throw name_not_found();
}

Attribute Reader::get_attribute(const std::string& name) const {
	LOG(debug) << "Getting attribute " << name << " from file " << get_filename();
	std::vector<std::string> parts = tokenize(name, "/");
	if (parts.size() == 1) {
		throw invalid_argument("attribute name does not name a group or dataset");
	}

	try {
		return _get_attribute(hdf5_file.value(), parts);
	} catch (const attribute_not_found&) {
		std::ostringstream os;
		os << "Attribute " << name << " doesn't exist in " << get_filename();
		throw invalid_data(os.str());
	} catch (const name_not_found&) {
		std::ostringstream os;
		os << "Name " << name << " names no object in file " << get_filename();
		throw invalid_data(os.str());
	}
}

}  // namespace hdf5

}  // namespace shark