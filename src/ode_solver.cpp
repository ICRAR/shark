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
 * ODE Solver class implementation
 */

#include <sstream>

#include "exceptions.h"
#include "logging.h"
#include "ode_solver.h"

#include <gsl/gsl_errno.h>

using namespace std;

namespace shark {

ODESolver::ODESolver(ode_evaluator evaluator, size_t dimension, double precision, void *params) :
	ode_system(),
	driver()
{
	// "42" is a dummy hstart, we need something != 0
	ode_system = std::unique_ptr<gsl_odeiv2_system>(new gsl_odeiv2_system{evaluator, nullptr, dimension, params});
	driver.reset(gsl_odeiv2_driver_alloc_y_new(ode_system.get(), gsl_odeiv2_step_rkck, 42, 0, precision));
}

void ODESolver::evolve(std::vector<double> &y, double delta_t)
{
	double t0 = 0;
	double t1 = t0 + delta_t;
	gsl_odeiv2_driver_reset_hstart(driver.get(), delta_t);
	int status = gsl_odeiv2_driver_apply(driver.get(), &t0, t1, y.data());

	// TODO: add compiler-dependent likelihood macro
	if (status == GSL_SUCCESS) {
		return;
	}

	//TEST: forcing integration to finish regardless of accuracy issue in three cases.

	ostringstream os;
	os << "Error while solving ODE system: ";
	if (status == GSL_FAILURE) {
		os << "step size decreases below machine precision ";
		LOG(warning) << "ODE: step size decreases below machine precision. Will force integration to finish regardless of desired accuracy not reached.";
		return;
	}
	if (status == GSL_ENOPROG) {
		os << "step size dropped below minimum value";
		LOG(warning) << "ODE:step size dropped below minimum value. Will force integration to finish regardless of desired accuracy not reached.";
		return;
	}
	else if (status == GSL_EBADFUNC) {
		os << "user function signaled an error";
		gsl_odeiv2_driver_reset(driver.get());
	}
	else if (status == GSL_EMAXITER) {
		os << "maximum number of steps reached";
		LOG(warning) << "ODE:maximum number of steps reached. Will force integration to finish regardless of desired accuracy not reached.";
		return;
	}
	else {
		os << "unexpected GSL error: " << gsl_strerror(status);
	}

	throw math_error(os.str());
}

unsigned long int ODESolver::num_evaluations()
{
	return driver->n;
}

}  // namespace shark
