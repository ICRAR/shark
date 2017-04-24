//
// Importer options class
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

#include "importer/options.h"

using namespace std;

namespace shark {

namespace importer {

Options::Options(const string &filename) :
	shark::Options(filename),
	descendants_format(ASCII),
	descendants_file("descendants.txt"),
	tree_format(TREES_VELOCIRAPTOR),
	tree_dir("."),
	first_snapshot(0),
	last_snapshot(0)
{
	load("input.tree_dir", tree_dir);
	load("input.tree_format", tree_format);
	load("input.descendants_format", descendants_format);
	load("input.descendants", descendants_file);
	load("input.first_snapshot", first_snapshot);
	load("input.last_snapshot", last_snapshot);
}

} // namespace importer


namespace detail {

template <>
importer::Options::tree_format_t Helper<importer::Options::tree_format_t>::get(const std::string &name, const std::string &value) {
	if ( value == "velociraptor" ) {
		return importer::Options::TREES_VELOCIRAPTOR;
	}
	else if ( value == "nifty" ) {
		return importer::Options::TREES_NIFTY;
	}
	std::ostringstream os;
	os << "tree_format option value invalid: " << value;
	throw invalid_option(os.str());
}

} // namespace detail

} // namespace shark