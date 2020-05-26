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

#include "gsl_utils.h"
#include "integrator.h"


namespace shark {

Integrator::Integrator(size_t max_intervals) :
	max_intervals(max_intervals),
	num_intervals(0)
{
	init_gsl_objects();
}

Integrator::Integrator(const Integrator &other) :
	max_intervals(other.max_intervals),
	num_intervals(other.num_intervals)
{
	init_gsl_objects();
}

void Integrator::init_gsl_objects()
{
	workspace.reset(gsl_integration_workspace_alloc(max_intervals));
}

double Integrator::integrate(func_t f, void *params, double from, double to, double epsabs, double epsrel)
{
	gsl_function F;
	F.function = f;
	F.params = params;

	double result, abserr;
	// Adopt a 15 point Gauss-Kronrod rule
	int key = GSL_INTEG_GAUSS15;
	gsl_invoke(gsl_integration_qag, &F, from, to, epsabs, epsrel, max_intervals, key, workspace.get(), &result, &abserr);

	num_intervals += workspace->size;
	return result;
}

std::size_t Integrator::get_num_intervals()
{
	return num_intervals;
}

void Integrator::reset_num_intervals()
{
	num_intervals = 0;
}

}  // namespace shark
