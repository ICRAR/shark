//
// Integrator class implementation
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

#include "integrator.h"


namespace shark {

Integrator::Integrator(size_t max_intervals) :
	workspace(),
	max_intervals(max_intervals),
	num_intervals(0)
{
	init_gsl_objects();
}

Integrator::Integrator(Integrator &&other) :
	workspace(),
	max_intervals(other.max_intervals),
	num_intervals(other.num_intervals)
{
	std::swap(workspace, other.workspace);
}

Integrator::Integrator(const Integrator &other) :
	workspace(),
	max_intervals(other.max_intervals),
	num_intervals(other.num_intervals)
{
	init_gsl_objects();
}

Integrator::~Integrator()
{
	if (workspace) {
		gsl_integration_workspace_free(workspace.release());
	}
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
	// Adopt a 15 point Gauss-Kronrod rule (key = 1).
	int key = 1;
	//gsl_integration_qags(&F, from, to, epsabs, epsrel, max_intervals, workspace.get(), &result, &abserr);
	gsl_integration_qag(&F, from, to, epsabs, epsrel, max_intervals, key, workspace.get(), &result, &abserr);
	num_intervals += workspace->size;
	return result;
}

unsigned long int Integrator::get_num_intervals()
{
	return num_intervals;
}

void Integrator::reset_num_intervals()
{
	num_intervals = 0;
}

}  // namespace shark
