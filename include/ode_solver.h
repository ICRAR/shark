//
// ODE solver class definition
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

#ifndef SHARK_ODE_SOLVER_H_
#define SHARK_ODE_SOLVER_H_

#include <vector>

#include <gsl/gsl_odeiv2.h>

namespace shark {

/**
 * A solver of ODE systems
 *
 * The ODE system solved by this class is defined in terms of an evaluation
 * function, a list of initial values `y0`, a time zero `t0` and a `delta_t`
 * parameter. After defining the solver, each call to the `next` method will
 * evolve the system, evaluating it at `t = t + delta_t`, with `t` starting at
 * `t0`, and returning the new values.
 */
class ODESolver {

public:

	/**
	 * The definition that ODE evaluators must follow
	 * @param
	 * @param y
	 * @param f
	 * @param A pointer to any user-provided data
	 * @return
	 */
	typedef int (*ode_evaluator)(double, const double y[], double f[], void *);

	/**
	 * Creates a new ODESolver
	 *
	 * @param y0 The initial values for the ODE system. It should contain as
	 * many values as those produced by `evaluator`
	 * @param t0 The time associated with the initial values `y0`.
	 * @param delta_t The time difference used to evolve the system on each
	 * evaluation
	 * @param evaluator The function evaluating the system at time `t`
	 * @param precision The precision to use for the adaptive step sizes.
	 */
	ODESolver(std::vector<double> y0, double t0, double delta_t,
	          ode_evaluator evaluator, double precision = 1e-6);

	/**
	 * Destructs this solver and frees up all resources associated with it
	 */
	~ODESolver();

	/**
	 * Steps the time axis in `delta_t`, and evaluates the underlying ODE system
	 * in `t`.
	 *
	 * @return The `y` values for the evaluation of the system at `t`.
	 */
	std::vector<double> next();

	/**
	 * Returns the current time `t` at which the system is sitting
	 * @return
	 */
	inline double current_t() {
		return t;
	}

private:
	std::vector<double> y;
	double t;
	double t0;
	double delta_t;
	unsigned int step;
	gsl_odeiv2_system system;
	gsl_odeiv2_driver *driver;
};

}  // namespace shark

#endif // SHARK_ODE_SOLVER_H_