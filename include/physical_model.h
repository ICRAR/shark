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
 * Physical model definition
 */

#ifndef SHARK_SYSTEM_H_
#define SHARK_SYSTEM_H_

#include <memory>
#include <sstream>
#include <stdexcept>

#include <gsl/gsl_odeiv2.h>
#include <recycling.h>
#include "components.h"
#include "gas_cooling.h"
#include "numerical_constants.h"
#include "ode_solver.h"
#include "stellar_feedback.h"
#include "star_formation.h"

namespace shark {

template <int NC>
class PhysicalModel {

public:

	/**
	 * The set of parameters passed down to the ODESolver. It includes the
	 * physical model itself, the galaxy and subhalo being evolved on each call,
	 * and other various values.
	 */
	struct solver_params {
		PhysicalModel<NC> &model;
		bool burst;
		double rgas;
		double rstar;
		double mcoolrate;
		double jcold_halo;
		double delta_t;
		double redshift;
		double vsubh;
		double vgal;
	};

	PhysicalModel(
			double ode_solver_precision,
			ODESolver::ode_evaluator evaluator,
			GasCooling gas_cooling) :
		params {*this, false, 0., 0., 0., 0., 0., 0., 0.},
		starburst_params {*this, true, 0., 0., 0., 0., 0., 0., 0.},
		ode_solver(evaluator, NC, ode_solver_precision, &params),
		starburst_ode_solver(evaluator, NC, ode_solver_precision, &starburst_params),
		ode_values(NC), starburst_ode_values(NC),
		gas_cooling(std::move(gas_cooling)),
		galaxy_ode_evaluations(0),
		galaxy_starburst_ode_evaluations(0)
	{
		// no-op
	}

	virtual ~PhysicalModel() = default;

	void evolve_galaxy(Subhalo &subhalo, Galaxy &galaxy, double z, double delta_t)
	{
		/**
		 * Parameters that are needed as input in the ode_solver:
		 * mcoolrate: gas cooling rate onto galaxy [Msun/Gyr/h]
		 * rgas: half-gas mass radius of the disk [Mpc/h]
		 * vgal: disk velocity at rgas [km/s]
		 * rstar: half-stellar mass radius of the disk [Mpc/h]
		 * vsubh: virial velocity of the host subhalo [km/s]
		 * jcold_halo: specific angular momentum of the cooling gas [Msun/h Mpc/h km/s]
		 * burst: boolean parameter indicating if this is a starburst or not.
		 */

		// Define cooling rate only in the case galaxy is central.
		params.mcoolrate = 0;
		if (galaxy.galaxy_type == Galaxy::CENTRAL) {
			params.mcoolrate = gas_cooling.cooling_rate(subhalo, galaxy, z, delta_t);
		}

		params.rgas = galaxy.disk_gas.rscale; //gas scale radius.
		params.vgal = galaxy.disk_gas.sAM / galaxy.disk_gas.rscale * constants::EAGLEJconv;

		// Catch cases where gas disk doesn't exist yet.
		if (params.rgas <= 0) {
			//In this case assign a scalelength due to the cooling gas.
			params.rgas = subhalo.cold_halo_gas.sAM / galaxy.vmax * constants::EAGLEJconv;
			params.vgal = galaxy.vmax;
		}

		params.rstar      = galaxy.disk_stars.rscale; //stellar scale radius.
		params.vsubh      = subhalo.Vvir;
		params.jcold_halo = subhalo.cold_halo_gas.sAM;
		params.delta_t = delta_t;
		params.redshift = z;

		from_galaxy(ode_values, subhalo, galaxy);
		ode_solver.evolve(ode_values, delta_t);
		galaxy_ode_evaluations += ode_solver.num_evaluations();
		to_galaxy(ode_values, subhalo, galaxy, delta_t);
	}

	void evolve_galaxy_starburst(Subhalo &subhalo, Galaxy &galaxy, double z, double delta_t, bool from_galaxy_merger)
	{

		/**
		 * Parameters that are needed as input in the ode_solver:
		 * mcoolrate: gas cooling rate onto galaxy [Msun/Gyr/h]. In the case of starbursts, this is \equiv 0
		 * rgas: half-gas mass radius of the bulge [Mpc/h]
		 * vgal: bulge velocity at rgas [km/s]
		 * rstar: half-stellar mass radius of the bulge [Mpc/h]
		 * vsubh: virial velocity of the host subhalo [km/s]
		 * jcold_halo: specific angular momentum of the cooling gas [Msun/h Mpc/h km/s]
		 * burst: boolean parameter indicating if this is a starburst or not.
		 */

		starburst_params.rgas = galaxy.bulge_gas.rscale; //gas scale radius.
		starburst_params.rstar = galaxy.bulge_stars.rscale; //stellar scale radius.
		starburst_params.vsubh = subhalo.Vvir;
		starburst_params.vgal = galaxy.bulge_gas.sAM / galaxy.bulge_gas.rscale;
		starburst_params.delta_t = delta_t;
		starburst_params.redshift = z;

		from_galaxy_starburst(starburst_ode_values, subhalo, galaxy);
		starburst_ode_solver.evolve(starburst_ode_values, delta_t);
		galaxy_starburst_ode_evaluations += starburst_ode_solver.num_evaluations();
		to_galaxy_starburst(starburst_ode_values, subhalo, galaxy, delta_t, from_galaxy_merger);
	}

	virtual void from_galaxy(std::vector<double> &y, const Subhalo &subhalo, const Galaxy &galaxy) = 0;
	virtual void to_galaxy(const std::vector<double> &y, Subhalo &subhalo, Galaxy &galaxy, double delta_t) = 0;

	virtual void from_galaxy_starburst(std::vector<double> &y, const Subhalo &subhalo, const Galaxy &galaxy) = 0;
	virtual void to_galaxy_starburst(const std::vector<double> &y, Subhalo &subhalo, Galaxy &galaxy, double delta_t, bool from_galaxy_merger) = 0;

	unsigned long int get_galaxy_ode_evaluations() {
		return galaxy_ode_evaluations;
	}

	unsigned long int get_galaxy_starburst_ode_evaluations() {
		return galaxy_starburst_ode_evaluations;
	}

	void reset_ode_evaluations() {
		galaxy_ode_evaluations = 0;
		galaxy_starburst_ode_evaluations = 0;
	}

private:
	solver_params params;
	solver_params starburst_params;
	ODESolver ode_solver;
	ODESolver starburst_ode_solver;
	std::vector<double> ode_values;
	std::vector<double> starburst_ode_values;
	GasCooling gas_cooling;
	unsigned long int galaxy_ode_evaluations;
	unsigned long int galaxy_starburst_ode_evaluations;
};

class BasicPhysicalModel : public PhysicalModel<17> {
public:
	BasicPhysicalModel(double ode_solver_precision,
			GasCooling gas_cooling,
			StellarFeedback stellar_feedback,
			StarFormation star_formation,
			RecyclingParameters recycling_parameters,
			GasCoolingParameters gas_cooling_parameters);

	void from_galaxy(std::vector<double> &y, const Subhalo &subhalo, const Galaxy &galaxy);
	void to_galaxy(const std::vector<double> &y, Subhalo &subhalo, Galaxy &galaxy, double delta_t);

	void from_galaxy_starburst(std::vector<double> &y, const Subhalo &subhalo, const Galaxy &galaxy);
	void to_galaxy_starburst(const std::vector<double> &y, Subhalo &subhalo, Galaxy &galaxy, double delta_t, bool from_galaxy_merger);

	StellarFeedback stellar_feedback;
	StarFormation star_formation;
	RecyclingParameters recycling_parameters;
	GasCoolingParameters gas_cooling_parameters;

	void reset_ode_evaluations() {
		PhysicalModel::reset_ode_evaluations();
		star_formation.reset_integration_intervals();
	}

	unsigned long int get_star_formation_integration_intervals() {
		return star_formation.get_integration_intervals();
	}

};

}  // namespace shark

#endif // SHARK_SYSTEM_H_
