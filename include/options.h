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

#ifndef SHARK_OPTIONS_H_
#define SHARK_OPTIONS_H_

#include <map>
#include <string>
#include <sstream>

#include "exceptions.h"
#include "logging.h"
#include "utils.h"

namespace shark {

namespace detail {

template <typename T>
class Helper {
public:
	static
	T get(const std::string &name, const std::string &value);
};

}  // namespace detail

/**
 * Base class to load user-provided options.
 */
class Options {

public:

	typedef std::map<std::string, std::string> options_t;

	enum file_format_t {
		HDF5,
		ASCII
	};

	/**
	 * A ctor that reads options from a file
	 *
	 * @param filename The name of the options file
	 */
	Options(const std::string &filename);

	/**
	 * Read the value of option `name`, if present, and set it in `value_holder`
	 * of type `T`.
	 *
	 * @param name The full option name (e.g., "group.name")
	 * @param value_holder The variable where the value will be stored
	 * @param optional Whether the option is mandatory or not. Defaults to `false`
	 * @tparam T The type of the object returned by this methods
	 */
	template <typename T>
	void load(const std::string &name, T &value_holder, bool mandatory = false) const {
		if ( mandatory or options.find(name) != options.end() ) {
			value_holder = get<T>(name);
		}
	}

protected:

	/**
	 * Read the value of option `name` and return it as an object of type `T`
	 * @param name The full option name (e.g., "group.name")
	 * @tparam T The type of the object returned by this method.
	 * @return
	 */
	template <typename T>
	T get(const std::string &name) const {
		options_t::const_iterator it = options.find(name);
		if ( it == options.end() ) {
			throw missing_option(name);
		}

		LOG(debug) << "Loading option " << name << " = " << it->second;
		return detail::Helper<T>::get(name, it->second);
	}

	options_t options;

};

}  // namespace shark

#endif // SHARK_OPTIONS_H_
