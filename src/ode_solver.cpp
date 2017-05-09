//
// ODE Solver class implementation
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

#include <vector>
#include <sstream>

#include "exceptions.h"
#include "ode_solver.h"

#include <gsl/gsl_errno.h>

using namespace std;

namespace shark {

ODESolver::ODESolver(std::vector<double> y0, double t0, double delta_t, ode_evaluator evaluator, double precision) :
	y(y0),
	t0(t0),
	t(t0),
	delta_t(delta_t),
	step(0)
{
	system = {evaluator, NULL, y0.size(), NULL};
	driver = gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rk8pd, precision, precision, 0.0);
}

ODESolver::~ODESolver() {
	gsl_odeiv2_driver_free(driver);
}

std::vector<double> ODESolver::evolve() {

	step++;
	double t_i = t0 + step*delta_t;
	int status = gsl_odeiv2_driver_apply(driver, &t, t_i, y.data());

	// TODO: add compiler-dependent likelihood macro
	if (status == GSL_SUCCESS) {
		return y;
	}

	ostringstream os;
	os << "Error while solving ODE system: ";
	if (status == GSL_FAILURE) {
		os << "step size decreases below machine precision";
	}
	if (status == GSL_ENOPROG) {
		os << "step size dropped below minimum value";
	}
	else if (status == GSL_EBADFUNC) {
		os << "user function signaled an error";
		gsl_odeiv2_driver_reset(driver);
	}
	else if (status == GSL_EMAXITER) {
		os << "maximum number of steps reached";
	}
	else {
		os << "unknown error " << status;
	}
	throw math_error(os.str());
}

}  // namespace shark