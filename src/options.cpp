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
#include <set>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <string>

#include "options.h"
#include "utils.h"

using namespace std;

namespace shark {

Options::Options(const string &name)
{

	ifstream f = open_file(name);
	string line;
	string option_group;

	while ( getline(f, line) ) {

		trim(line);

		// Skip blanks and comments
		if ( line.size() == 0 ) {
			continue;
		}
		if ( line[0] == '#' ) {
			continue;
		}

		if ( line[0] == '[' ) {
			if ( line[line.size() - 1] != ']' ) {
				ostringstream os;
				os << "Invalid group definition: " << line;
				throw invalid_option(os.str());
			}
			option_group = line.substr(1, line.size() - 2);
			continue;
		}

		string name, value, equals;
		istringstream iss(line);
		iss >> name >> equals >> value;

		if ( option_group.size() == 0 ) {
			cerr << "WARNING: No option group defined for option " << name << endl;
		}

		name = option_group + '.' + name;
		options[name] = value;
	}

}


namespace detail {

template <>
std::string Helper<std::string>::get(const std::string &name, const std::string &value) {
	return value;
}

template <>
Options::file_format_t Helper<Options::file_format_t>::get(const std::string &name, const std::string &value) {
	std::string lowered(value);
	lower(lowered);
	if ( lowered == "hdf5" ) {
		return Options::HDF5;
	}
	else if ( lowered == "ascii" ) {
		return Options::ASCII;
	}

	std::ostringstream os;
	os << name << " option value invalid: " << value << ". Either hdf5 or ascii was expected";
	throw invalid_option(os.str());
}

template<>
int Helper<int>::get(const std::string &name, const std::string &value) {
	try {
		return std::stoi(value);
	} catch (const std::invalid_argument &e) {
		std::ostringstream os;
		os << "Invalid value for option " << name << ": " << value << ". An integer was expected";
		throw invalid_option(os.str());
	}
}

template<>
float Helper<float>::get(const std::string &name, const std::string &value) {
	try {
		return std::stof(value);
	} catch (const std::invalid_argument &e) {
		std::ostringstream os;
		os << "Invalid value for option " << name << ": " << value << ". A float was expected";
		throw invalid_option(os.str());
	}
}

template<>
double Helper<double>::get(const std::string &name, const std::string &value) {
	try {
		return std::stod(value);
	} catch (const std::invalid_argument &e) {
		std::ostringstream os;
		os << "Invalid value for option " << name << ": " << value << ". A double was expected";
		throw invalid_option(os.str());
	}
}

template<>
bool Helper<bool>::get(const std::string &name, const std::string &value) {
	try {
		bool bool_val = true;
		std::istringstream is(value);
		is >> std::boolalpha;
		is >> bool_val;
		printf("%d\n", bool_val);
		return bool_val;
	} catch (const std::exception &e) {
		std::ostringstream os;
		os << "Invalid value for option " << name << ": " << value << ". A boolean (true/false) was expected";
		throw invalid_option(os.str());
	}
}

template<>
std::vector<double> Helper<std::vector<double>>::get(const std::string &name, const std::string &value) {
	try {
		std::vector<std::string> values_as_str = tokenize(value, " ");
		std::vector<double> values;
		for(auto value_as_str: values_as_str) {
			values.push_back(std::stod(value_as_str));
		}
		return values;
	} catch (const std::invalid_argument &e) {
		std::ostringstream os;
		os << "Invalid value for option " << name << ": " << value << ". A double was expected";
		throw invalid_option(os.str());
	}
}

template<>
std::set<int> Helper<std::set<int>>::get(const std::string &name, const std::string &value) {
	try {
		std::vector<std::string> values_as_str = tokenize(value, " ");
		std::set<int> values;
		for(auto value_as_str: values_as_str) {
			values.insert(std::stod(value_as_str));
		}
		return values;
	} catch (const std::invalid_argument &e) {
		std::ostringstream os;
		os << "Invalid value for option " << name << ": " << value << ". A double was expected";
		throw invalid_option(os.str());
	}
}

template<>
std::vector<int> Helper<std::vector<int>>::get(const std::string &name, const std::string &value) {
	try {
		std::vector<std::string> values_as_str = tokenize(value, " ");
		std::vector<int> values;
		for(auto value_as_str: values_as_str) {
			values.push_back(std::stod(value_as_str));
		}
		return values;
	} catch (const std::invalid_argument &e) {
		std::ostringstream os;
		os << "Invalid value for option " << name << ": " << value << ". A double was expected";
		throw invalid_option(os.str());
	}
}

}  // namespace detail

}  // namespace shark
