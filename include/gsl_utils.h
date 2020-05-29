//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2020
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
 */

#ifndef SHARK_GSL_UTILS_H
#define SHARK_GSL_UTILS_H


#include <cassert>
#include <utility>

#include <gsl/gsl_errno.h>

#include "exceptions.h"

namespace shark {

/**
 * Takes a GSL error number and turns it into a gsl_error exception
 * @param status the GSL status, must be something other than GSL_SUCCESS
 */
inline gsl_error to_gsl_error(int status)
{
	assert(status != GSL_SUCCESS);
	throw gsl_error(status, gsl_strerror(status));
}


/**
 * Invokes a GSL function and throws an exception if the return code is anything
 * other than GSL_SUCESS
 * @param func The GSL function to invoke
 * @param args The arguments to pass to the GSL function
 */
template<typename Func, typename ... Args>
void gsl_invoke(Func && func, Args ... args)
{
	int status = func(std::forward<Args>(args)...);
	if (status != GSL_SUCCESS) {
		throw to_gsl_error(status);
	}
}

} // namespace shark

#endif // SHARK_GSL_UTILS_H