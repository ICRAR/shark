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
 * Solver class headers
 */

#ifndef SHARK_SOLVER_H_
#define SHARK_SOLVER_H_

#include <memory>

#include <gsl/gsl_roots.h>

#include "utils.h"

namespace shark {

///
/// A class that find the roots of arbitrary functions
///
class Root_Solver {

public:

	using func_x = double (*)(double x, void *);

	explicit Root_Solver(size_t max_intervals) {}

	/// Copy constructor
	Root_Solver(const Root_Solver &other);

	double root_solver_function(func_x t, void *params, double from, double to, double epsabs, double epsrel);

//private:
//	std::unique_ptr<gsl_root_fsolver, gsl_root_fsolver_free> root_solver =
//			std::unique_ptr<gsl_root_fsolver, gsl_root_fsolver_free> (gsl_root_fsolver_alloc (gsl_root_fsolver_brent));

};

}  // namespace shark

#endif // SHARK_SOLVER_H_
