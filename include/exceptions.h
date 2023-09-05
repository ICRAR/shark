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
 * Common exception definitions
 */

#ifndef SHARK_EXCEPTIONS_H
#define SHARK_EXCEPTIONS_H

#include <exception>
#include <string>
#include <utility>

#include "components.h"

namespace shark {

/**
 * The mother of all SHArk exceptions
 */
class exception : public std::exception {

public:

	explicit exception(std::string what) :
		_what(std::move(what))
	{}

	const char* what() const noexcept override {
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
	explicit invalid_option(std::string what) : exception(std::move(what))
	{}
};

/**
 * An exception indicating that an invalid argument value has been given to a
 * function or a class. Like std::invalid_argument, but inside shark::
 */
class invalid_argument : public exception {
public:
	explicit invalid_argument(std::string what) : exception(std::move(what))
	{}
};

/**
 * An exception indicating that a required option value is missing
 */
class missing_option : public exception {
public:
	explicit missing_option(std::string what) : exception(std::move(what))
	{}
};

/**
 * An exception indicating that invalid data has been encountered
 */
class invalid_data : public exception {
public:
	explicit invalid_data(std::string what) : exception(std::move(what))
	{};
};

/**
 * An exception indicating that a structural component was expected
 * but not found in the simulation data.
 */
class component_not_found : public invalid_data {
public:
	explicit component_not_found(std::string what) : invalid_data(std::move(what))
	{};
};

/**
 * An exception indicating that a Halo was expected but not found
 */
class halo_not_found : public component_not_found {
public:
	halo_not_found(std::string what, halo_id_t halo_id) :
		component_not_found(std::move(what)),
		halo_id(halo_id) {}

	/**
	 * The ID of the Halo that could not be found.
	 */
	halo_id_t halo_id;
};

/**
 * An exception indicating that a Halo was expected but not found
 */
class subhalo_not_found : public component_not_found {
public:
	subhalo_not_found(std::string what, subhalo_id_t subhalo_id) :
		component_not_found(std::move(what)),
		subhalo_id(subhalo_id) {}

	/**
	 * The ID of the Subhalo that could not be found.
	 */
	subhalo_id_t subhalo_id;
};

/**
 * An exception indicating a mathematical error
 */
class math_error: public exception {
public:
	explicit math_error(std::string what) : exception(std::move(what))
	{};
};

/**
 * An exception indicating that a GSL error occurred
 */
class gsl_error : public math_error {
public:
	gsl_error(const char *reason, const char *file, int line, int gsl_errno, const char *errmsg) :
		math_error(std::string("GSL error at ") + file + ":" +
		           std::to_string(line) + ": " + errmsg +
		           " (" + std::to_string(gsl_errno) + ")"),
		reason(reason),
		file(file),
		line(line),
		gsl_errno(gsl_errno)
	{
	}

	gsl_error(int gsl_errno, const char *errmsg) :
		math_error(std::string("GSL error: ") + errmsg +
		           " (" + std::to_string(gsl_errno) + ")"),
		gsl_errno(gsl_errno)
	{
	}

	std::string get_reason() {
		return reason;
	}

	std::string get_file() {
		return file;
	}

	int get_line() {
		return line;
	}

	int get_gsl_errno() {
		return gsl_errno;
	}

private:
	std::string reason;
	std::string file;
	int line = 0;
	int gsl_errno;
};

}  // namespace shark

#endif // SHARK_EXCEPTIONS_H