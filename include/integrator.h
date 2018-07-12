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
 * Integrator class headers
 */

#include <memory>

#include <gsl/gsl_integration.h>

#ifndef SHARK_INTEGRATOR_H_
#define SHARK_INTEGRATOR_H_

namespace shark {

///
/// A class that integrates functions through different ranges
///
class Integrator {

public:

	typedef double (*func_t)(double x, void *);

	///
	/// Creates a new Integrator that will integrate using at most
	/// `max_intervals` intervals internally.
	///
	/// @param max_intervals
	///
	Integrator(size_t max_intervals);

	// Copy/movy constructors, destructor
	Integrator(const Integrator &other);
	Integrator(Integrator &&other);
	~Integrator();

	///
	/// Integrates function `f` with parameters `params` between `from` and `to`
	/// using the indicated error tolerances.
	///
	double integrate(func_t f, void *params, double from, double to, double epsabs, double epsrel);

	///
	/// Returns the number of internal intervals used during all integrations
	/// so far, or since the last call to reset_num_intervals.
	/// The number of intervals is an indication of how many times the functions
	/// being integrated have been called.
	///
	unsigned long int get_num_intervals();

	///
	/// Reset the number of intervals count.
	///
	void reset_num_intervals();

private:
	std::unique_ptr<gsl_integration_workspace> workspace;
	size_t max_intervals;
	unsigned long int num_intervals;

	void init_gsl_objects();
};

}  // namespace shark

#endif // SHARK_INTEGRATOR_H_