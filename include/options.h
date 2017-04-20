//
// Header file for the options class
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

#ifndef SHARK_OPTIONS
#define SHARK_OPTIONS

#include <string>
#include <stdexcept>

class invalid_option : public std::runtime_error {
public:
	invalid_option(const std::string &what) : std::runtime_error(what) {}
};

namespace shark {

class Options {

public:

	enum tree_format_t {
		TREES_VELOCIRAPTOR,
		TREES_NIFTY
	};

	enum descendants_format_t {
		DESCENDANTS_HDF5,
		DESCENDANTS_ASCII
	};

	/**
	 * Constructor that provides default values for the options.
	 */
	Options();

	/**
	 * Reads the given file and parses out the options contained in it.
	 *
	 * @param name The name of the options file
	 * @return An Options object populated with the options contained in the file
	 */
	static const Options from_file(const std::string &name);

	descendants_format_t descendants_format;
	std::string descendants_file;
	tree_format_t tree_format;
	std::string tree_dir;
	int first_snapshot;
	int last_snapshot;

private:
	static int read_int(const std::string &name, const std::string &value);
};

}  // namespace shark

#endif // SHARK_OPTIONS