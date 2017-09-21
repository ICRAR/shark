//
// Integrator class headers
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

	Integrator(size_t max_samples);

	// Copy/movy constructors, destructor
	Integrator(const Integrator &other);
	Integrator(Integrator &&other);
	~Integrator();

	double integrate(func_t f, void *params, double from, double to, double epsabs, double epsrel);

private:
	std::unique_ptr<gsl_integration_workspace> workspace;
	size_t max_samples;

	void init_gsl_objects();
};

}  // namespace shark

#endif // SHARK_INTEGRATOR_H_