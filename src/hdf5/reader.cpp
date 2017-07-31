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

#include <string>
#include <vector>

#include "hdf5/reader.h"

using namespace std;

namespace shark {

namespace hdf5 {

H5::DataSet Reader::get_dataset(const string &name) const {

	// The name might contains slashes, so we can navigate through
	// a hierarchy of groups/datasets
	const vector<string> parts = tokenize(name, "/");

	return get_dataset(parts);
}

H5::DataSet Reader::get_dataset(const std::vector<std::string> &path) const {

	// only the attribute name, read directly and come back
	if( path.size() == 1 ) {
		return hdf5_file.openDataSet(path[0]);
	}

	// else there's a path to follow, go for it!
	H5::Group group = hdf5_file.openGroup(path.front());
	vector<string> group_paths(path.begin() + 1, path.end() - 1);
	for(auto const &path: group_paths) {
		group = group.openGroup(path);
	}

	return group.openDataSet(path.back());
}

H5::Attribute Reader::get_attribute(const string &name) const {

	// The name might contains slashes, so we can navigate through
	// a hierarchy of groups/datasets
	std::vector<std::string> parts = tokenize(name, "/");

	// only the attribute name, read directly and come back
	if( parts.size() == 1 ) {
		return hdf5_file.openAttribute(name);
	}

	// otherwise the path looks like
	//  [group1/.../groupN/]dataset/attribute
	const auto attr_name = parts.back();
	parts.pop_back();

	auto dataset = get_dataset(parts);
	return dataset.openAttribute(attr_name);
}

}  // namespace hdf5

}  // namespace shark