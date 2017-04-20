//
// Options class implementation
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

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <string>

#include "options.h"
#include "utils.h"

using namespace std;

namespace shark {

Options::Options() :
	descendants_file("descendants.txt"),
	descendants_format(DESCENDANTS_HDF5),
	tree_dir("."),
	tree_format(TREES_VELOCIRAPTOR),
	first_snapshot(0),
	last_snapshot(1)
{
	// no-op
}

const Options Options::from_file(const string &name)
{
	Options opts;

	ifstream f = open_file(name);
	string line;
	while ( getline(f, line) ) {

		trim(line);

		// Skip blanks and comments
		if ( line.size() == 0 ) {
			continue;
		}
		if ( line[0] == '#' ) {
			continue;
		}

		string name, value, equals;
		istringstream iss(line);
		iss >> name >> equals >> value;

		if ( name == "tree_dir" ) {
			opts.tree_dir = value;
		}
		else if ( name == "tree_format" ) {
			lower(value);
			if ( value == "velociraptor" ) {
				opts.tree_format = Options::TREES_VELOCIRAPTOR;
			}
			else if ( value == "nifty" ) {
				opts.tree_format = Options::TREES_NIFTY;
			}
			else {
				ostringstream os;
				os << "tree_format option value invalid: " << value;
				throw invalid_option(os.str());
			}
		}
		else if ( name == "descendants_format" ) {
			lower(value);
			if ( value == "hdf5" ) {
				opts.descendants_format = Options::DESCENDANTS_HDF5;
			}
			else if ( value == "ascii" ) {
				opts.descendants_format = Options::DESCENDANTS_ASCII;
			}
			else {
				ostringstream os;
				os << "descendants_format option value invalid: " << value;
				throw invalid_option(os.str());
			}
		}
		else if ( name == "descendants") {
			opts.descendants_file = value;
		}
		else if ( name == "first_snapshot" ) {
			opts.first_snapshot = read_int(name, value);
		}
		else if ( name == "last_snapshot" ) {
			opts.last_snapshot = read_int(name, value);
		}
		else {
			cout << "Ignoring unknown option: " << name << endl;
		}
	}

	return opts;

}

int Options::read_int(const std::string &name, const std::string &value)
{
	try {
		return stoi(value);
	} catch (const invalid_argument &e) {
		ostringstream os;
		os << "Invalid value for option " << name << ": " << value << ". An integer was expected";
		throw invalid_option(os.str());
	}
}

}  // namespace shark