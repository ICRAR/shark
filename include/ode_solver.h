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
 * ODE solver class definition
 */

#ifndef SHARK_ODE_SOLVER_H_
#define SHARK_ODE_SOLVER_H_

#include <memory>
#include <vector>

#include <gsl/gsl_odeiv2.h>

#include "utils.h"

namespace shark {

/// A deleter of gsl_odeiv2_driver objects
using gsl_odeiv2_driver_deleter = deleter<gsl_odeiv2_driver, gsl_odeiv2_driver_free>;

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
	using ode_evaluator = int (*)(double, const double y[], double f[], void *);

	/// Simple structure holding a user-provided evaluation function and params
	struct evaluator_and_params {
		ode_evaluator evaluator;
		void *user_params;
	};

	/**
	 * Creates a new ODESolver
	 *
	 * @param evaluator The function that evaluates the ODE system
	 * @param dimension The dimensionality of the ODE system to solve
	 * @param precision The precision to use for the adaptive step sizes.
	 * @param params The parameters to pass down to ``evaluator``
	 */
	ODESolver(ode_evaluator evaluator, size_t dimension, double precision, void *params);

	/**
	 * Evolves the ODE system from 0 to ``delta_t``
	 *
	 * @param y The values of the system at ``t = 0``. After returning the vector
	 *  contains the values at ``delta_t``.
	 * @param delta_t The amount of time the system is evolved for
	 */
	void evolve(std::vector<double> &y, double delta_t);

	/**
	 * Returns the number of times that the internal ODE system has been
	 * evaluated so far.
	 *
	 * @return The number of times the ODE system has been evaluated.
	 */
	std::size_t num_evaluations();

private:
	std::unique_ptr<gsl_odeiv2_system> ode_system;
	std::unique_ptr<gsl_odeiv2_driver, gsl_odeiv2_driver_deleter> driver;
	evaluator_and_params wrapped_params;
};

}  // namespace shark

#endif // SHARK_ODE_SOLVER_H_