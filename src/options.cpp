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
 * Options class implementation
 */

#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <string>

#include "logging.h"
#include "naming_convention.h"
#include "options.h"
#include "utils.h"

namespace shark {

Options::Options(const std::string &fname)
{
	add_file(fname);
}

void Options::add_file(const std::string &fname) {

	LOG(info) << "Loading options from " << fname;
	std::ifstream f = open_file(fname);
	std::string line;
	std::string option_group;

	while ( std::getline(f, line) ) {

		trim(line);

		// Skip blanks and comments
		if (line.empty()) {
			continue;
		}
		if ( line[0] == '#' ) {
			continue;
		}

		if ( line[0] == '[' ) {
			if ( line[line.size() - 1] != ']' ) {
				std::ostringstream os;
				os << "Invalid group definition: " << line;
				throw invalid_option(os.str());
			}
			option_group = line.substr(1, line.size() - 2);
			continue;
		}

		std::string name;
		std::string value;
		parse_option(line, name, value);

		if (option_group.empty()) {
			LOG(warning) << "No option group defined for option " << name;
		}

		std::string full_name = option_group;
		full_name += '.';
		full_name += name;
		store_option(full_name, value);
	}

}

void Options::add(const std::string &optspec)
{
	std::string name;
	std::string value;
	parse_option(optspec, name, value);
	store_option(name, value);
}

void Options::check_valid_name(const std::string &name)
{
	auto tokens = tokenize(name, ".");
	for (auto &token: tokens) {
		if (!follows_convention(token, naming_convention::SNAKE_CASE)) {
			std::ostringstream os;
			os << "A part of option " << name << " does not follow the ";
			os << "snake_case naming convention: " << token;
			throw invalid_option(os.str());
		}
	}

}

void Options::store_option(const std::string &name, const std::string &value)
{
	check_valid_name(name);
	if (options.find(name) == options.end()) {
		LOG(info) << "Loading new option: " << name << " = " << value;
	}
	else {
		LOG(info) << "Overwriting option: " << name << " = " << value;
	}
	options[name] = value;
}

void Options::parse_option(const std::string &optspec, std::string &name, std::string &value)
{
	auto tokens = tokenize(optspec, "=");
	if (tokens.size() < 2) {
		std::ostringstream os;
		os << "Option " << optspec << " has no value (should be name = value)";
		throw invalid_option(os.str());
	}

	name = tokens[0];
	value = tokens[1];
	trim(name);
	trim(value);

	// Empty after trimming? no way
	if (name.empty()) {
		std::ostringstream os;
		os << "Option without name: " << value;
		throw invalid_option(os.str());
	}
	if (value.empty()) {
		std::ostringstream os;
		os << "Option without value: " << name;
		throw invalid_option(os.str());
	}

	// Remove possible comment in value
	std::string::size_type hash;
	if ((hash = value.find('#')) != std::string::npos) {
		value = value.substr(0, hash);
		trim(value);
	}

}

template <typename T>
struct type_name {
	constexpr static const char *name = "unknown";
};

template <>
struct type_name<unsigned int> {
	constexpr static const char *name = "unsigned int";
};

template <>
struct type_name<int> {
	constexpr static const char *name = "int";
};

template <>
struct type_name<float> {
	constexpr static const char *name = "float";
};

template <>
struct type_name<double> {
	constexpr static const char *name = "double";
};

template <typename T>
static inline
T _from_string(const std::string &val);

template <>
float _from_string<float>(const std::string &val)
{
	return std::stof(val);
}

template <>
double _from_string<double>(const std::string &val)
{
	return std::stod(val);
}

template <>
int _from_string<int>(const std::string &val)
{
	return std::stoi(val);
}

template <>
unsigned int _from_string<unsigned int>(const std::string &val)
{
	return std::stoul(val);
}

template <typename T>
T _builtin_from_string(const std::string &name, const std::string &val)
{
	try {
		return _from_string<T>(val);
	} catch (const std::invalid_argument &) {
		std::ostringstream os;
		os << "Invalid value for option " << name << ": " << val << ". "
		   << type_name<T>::name << " value was expected";
		throw invalid_option(os.str());
	}
}

template <typename Cont>
typename std::enable_if<std::is_integral<typename Cont::value_type>::value, Cont>::type
_read_ranges(const std::string &name, const std::string &value, const std::string &sep = " ")
{

	using T = typename Cont::value_type;

	Cont values;
	for(auto value_or_range: tokenize(value, sep)) {

		trim(value_or_range);
		if (value_or_range.empty()) {
			continue;
		}

		// a dash found neither at the beginning, nor at the end
		// means that we have a range specification
		auto pos = value_or_range.find_last_of('-');
		if (pos != 0 && pos != value_or_range.size() && pos != std::string::npos) {
			auto first_s = value_or_range.substr(0, pos);
			auto last_s = value_or_range.substr(pos + 1);

			auto first = _builtin_from_string<T>(name, first_s);
			auto last = _builtin_from_string<T>(name, last_s);

			// Can't find a more intelligent way of doing this, sorry...
			if (first < last) {
				for(auto i = first; i <= last; i++) {
					std::inserter(values, values.end()) = i;
				}
			}
			else {
				for(auto i = first; i >= last; i--) {
					std::inserter(values, values.end()) = i;
				}
			}
			continue;
		}

		// A normal value
		std::inserter(values, values.end()) = _builtin_from_string<T>(name, value_or_range);
	}

	return values;
}


template <>
std::string Options::get<std::string>(const std::string &name, const std::string &value) const {
	return value;
}

template <>
Options::file_format_t Options::get<Options::file_format_t>(const std::string &name, const std::string &value) const {
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
int Options::get<int>(const std::string &name, const std::string &value) const {
	return _builtin_from_string<int>(name, value);
}


template<>
unsigned int Options::get<unsigned int>(const std::string &name, const std::string &value) const {
	return _builtin_from_string<unsigned int>(name, value);
}

template<>
float Options::get<float>(const std::string &name, const std::string &value) const {
	return _builtin_from_string<float>(name, value);
}

template<>
double Options::get<double>(const std::string &name, const std::string &value) const {
	return _builtin_from_string<double>(name, value);
}

template<>
bool Options::get<bool>(const std::string &name, const std::string &value) const {
	try {
		bool bool_val = true;
		std::istringstream is(value);
		is >> std::boolalpha;
		is >> bool_val;
		return bool_val;
	} catch (const std::exception &) {
		std::ostringstream os;
		os << "Invalid value for option " << name << ": " << value << ". A boolean (true/false) was expected";
		throw invalid_option(os.str());
	}
}

template<>
std::vector<double> Options::get<std::vector<double>>(const std::string &name, const std::string &value) const {
	std::vector<std::string> values_as_str = tokenize(value, " ");
	std::vector<double> values;
	std::transform(values_as_str.begin(), values_as_str.end(), std::back_inserter(values), [&name](const std::string &s) {
		return _builtin_from_string<double>(name, s);
	});
	return values;
}

template<>
std::set<int> Options::get<std::set<int>>(const std::string &name, const std::string &value) const {
	return _read_ranges<std::set<int>>(name, value);
}

template<>
std::vector<int> Options::get<std::vector<int>>(const std::string &name, const std::string &value) const {
	return _read_ranges<std::vector<int>>(name, value);
}

template<>
std::vector<unsigned int> Options::get<std::vector<unsigned int>>(const std::string &name, const std::string &value) const {
	return _read_ranges<std::vector<unsigned int>>(name, value);
}

}  // namespace shark
