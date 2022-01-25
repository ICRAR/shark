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
 * Integrator class implementation
 */

#include <string>
#include <gsl/gsl_integration.h>

#include "gsl_utils.h"
#include "logging.h"
#include "root_solver.h"

namespace shark {


/*struct root_solver_data_t {
	Root_Solver::func_x func;
	void *params;
	int status;
};*/


double Root_Solver::root_solver_function(func_x f, void *params, double from, double to, double epsabs, double epsrel)
{

	int status;
	int iter = 0, max_iter = 1000;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double result = 0;

	gsl_function F;
	F.function = f;
	F.params = params;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set (s, &F, from, to);

	do
	{
		iter++;
		status = gsl_root_fsolver_iterate (s);
		result = gsl_root_fsolver_root (s);
		from = gsl_root_fsolver_x_lower (s);
		to = gsl_root_fsolver_x_upper (s);
		status = gsl_root_test_interval (from, to, epsabs, epsrel);
	}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free (s);

	return result;
}


}  // namespace shark
