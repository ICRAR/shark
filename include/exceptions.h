//
// Common exception definitions
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

#ifndef SHARK_EXCEPTIONS_H
#define SHARK_EXCEPTIONS_H

#include <exception>
#include <string>

#include "components.h"

namespace shark {

/**
 * The mother of all SHArk exceptions
 */
class exception : public std::exception {

public:

	exception(const std::string &what) :
		std::exception(),
		_what(what)
	{}

	virtual const char* what() const noexcept {
		return _what.c_str();
	}

private:
	std::string _what;
};

/**
 * An exception indicating that an invalid option value has been given by the
 * user.
 */
class invalid_option : public exception {
public:
	invalid_option(const std::string &what) : exception(what) {}
};

/**
 * An exception indicating that a required option value is missing
 */
class missing_option : public exception {
public:
	missing_option(const std::string &what) : exception(what) {}
};

/**
 * An exception indicating that invalid data has been encountered
 */
class invalid_data : public exception {
public:
	invalid_data(const std::string &what) : exception(what) {};
};

/**
 * An exception indicating that a structural component was expected
 * but not found in the simulation data.
 */
class component_not_found : public invalid_data {
public:
	component_not_found(const std::string &what) : invalid_data(what) {};
};

/**
 * An exception indicating that a Halo was expected but not found
 */
class halo_not_found : public component_not_found {
public:
	halo_not_found(const std::string &what, Halo::id_t halo_id) :
		component_not_found(what),
		halo_id(halo_id) {}

	/**
	 * The ID of the Halo that could not be found.
	 */
	Halo::id_t halo_id;
};

/**
 * An exception indicating that a Halo was expected but not found
 */
class subhalo_not_found : public component_not_found {
public:
	subhalo_not_found(const std::string &what, Subhalo::id_t subhalo_id) :
		component_not_found(what),
		subhalo_id(subhalo_id) {}

	/**
	 * The ID of the Subhalo that could not be found.
	 */
	Subhalo::id_t subhalo_id;
};

/**
 * An exception indicating a mathematical error
 */
class math_error: public exception {
public:
	math_error(const std::string &what) : exception(what) {};
};

}  // namespace shark

#endif // SHARK_EXCEPTIONS_H