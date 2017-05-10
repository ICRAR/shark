//
// System class definition
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

#ifndef SHARK_SYSTEM_H_
#define SHARK_SYSTEM_H_

#include <memory>
#include <sstream>
#include <stdexcept>

#include <gsl/gsl_odeiv2.h>

#include "ode_solver.h"

namespace shark {

template <int NC>
class System {

public:
	System(double t0, double delta_t, double ode_solver_precision, ODESolver::ode_evaluator evaluator) :
		ode_system(std::shared_ptr<gsl_odeiv2_system>(new gsl_odeiv2_system{evaluator, NULL, NC, this})),
		t0(t0),
		delta_t(delta_t),
		ode_solver_precision(ode_solver_precision)
	{
		// no-opsrc/utils.cpp
	}

	virtual ~System()
	{
		// no-op
	}

	ODESolver get_solver(const std::vector<double> &y0) const {
		if (y0.size() != NC) {
			std::ostringstream os;
			os << "# initial values != ODE components: " << y0.size() << " != " << NC;
			throw std::invalid_argument(os.str());
		}
		return ODESolver(y0, t0, delta_t, ode_solver_precision, ode_system);
	}

private:
	std::shared_ptr<gsl_odeiv2_system> ode_system;
	double t0;
	double delta_t;
	double ode_solver_precision;
};

class BasicSystem : public System<6> {
public:
	BasicSystem(double t0, double delta_t, double ode_solver_precision);
};

class BasicSystemWithSatellites : public System<7> {
public:
	BasicSystemWithSatellites(double t0, double delta_t, double ode_solver_precision);
};

}  // namespace shark

#endif // SHARK_SYSTEM_H_