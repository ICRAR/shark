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

#include "components.h"
#include "ode_solver.h"
#include "parameters.h"
#include "stellar_feedback.h"
#include "star_formation.h"
#include "gas_cooling.h"

namespace shark {

template <int NC>
class PhysicalModel {

public:
	PhysicalModel(double ode_solver_precision, ODESolver::ode_evaluator evaluator,
			GasCooling gas_cooling) :
		ode_system(std::shared_ptr<gsl_odeiv2_system>(new gsl_odeiv2_system{evaluator, NULL, NC, this})),
		ode_solver_precision(ode_solver_precision),
		gas_cooling(gas_cooling)
	{
		// no-opsrc/utils.cpp
	}

	virtual ~PhysicalModel()
	{
		// no-op
	}

	ODESolver get_solver(double t0, double delta_t, const std::vector<double> &y0) {
		if (y0.size() != NC) {
			std::ostringstream os;
			os << "# initial values != ODE components: " << y0.size() << " != " << NC;
			throw std::invalid_argument(os.str());
		}
		return ODESolver(y0, t0, delta_t, ode_solver_precision, ode_system);
	}

	void evolve_galaxy(Subhalo &subhalo, Galaxy &galaxy)
	{
		double mcool = gas_cooling.cooling_rate(subhalo.halo_gas.mass, subhalo.Mvir, subhalo.Vvir, subhalo.halo_gas.mass_metals);
		std::vector<double> y0 = from_galaxy(subhalo, galaxy);
		std::vector<double> y1 = get_solver(t0, delta_t, y0).evolve();
		to_galaxy(y1, subhalo, galaxy);
	}

	virtual std::vector<double> from_galaxy(Subhalo &subhalo, Galaxy &galaxy) = 0;
	virtual void to_galaxy(const std::vector<double> &y, Subhalo &subhalo, Galaxy &galaxy) = 0;

private:
	std::shared_ptr<gsl_odeiv2_system> ode_system;
	double ode_solver_precision;
	GasCooling gas_cooling;
};

class BasicPhysicalModel : public PhysicalModel<6> {
public:
	BasicPhysicalModel(double ode_solver_precision,
			GasCooling gas_cooling,
			StellarFeedback stellar_feedback,
			StarFormation star_formation,
			RecyclingParameters recycling_parameters);

	virtual std::vector<double> from_galaxy(Subhalo &subhalo, Galaxy &galaxy) = 0;
	virtual void to_galaxy(const std::vector<double> &y, Subhalo &subhalo, Galaxy &galaxy) = 0;

	StellarFeedback stellar_feedback;
	StarFormation star_formation;
	RecyclingParameters recycling_parameters;
};

}  // namespace shark

#endif // SHARK_SYSTEM_H_
