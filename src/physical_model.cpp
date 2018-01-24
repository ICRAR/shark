//
// System classes implementation
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

#include <cmath>
#include <memory>

#include "physical_model.h"
#include "logging.h"
#include "numerical_constants.h"

namespace shark {

static
int basic_physicalmodel_evaluator(double t, const double y[], double f[], void *data) {

	/*
	 * f[0]: stellar mass of galaxy.
	 * f[1]: cold gas mass of galaxy.
	 * f[2]: hot gas mass;
	 * f[3]: ejected gas mass;
	 * f[4]: metals locked in the stellar mass of galaxies.
	 * f[5]: metals locked in the cold gas mass of galaxies.
	 * f[6]: metals locked in the hot gas mass.
	 * f[7]: metals locked in the ejected gas mass.
	 * f[8]: total stellar mass formed (without recycling included).
	 */

	auto params= reinterpret_cast<BasicPhysicalModel::solver_params *>(data);
	BasicPhysicalModel &model = dynamic_cast<BasicPhysicalModel &>(params->model);

	double R = model.recycling_parameters.recycle; /*recycling fraction of newly formed stars*/

	double yield = model.recycling_parameters.yield; /*yield of newly formed stars*/

	double mcoolrate = params->mcoolrate; /*cooling rate in units of Msun/Gyr*/

	double SFR = model.star_formation.star_formation_rate(y[1], y[0], params->rgas, params->rstar, params->redshift, params->burst);

	double beta = model.stellar_feedback.outflow_rate(SFR, params->v, params->redshift); /*mass loading parameter*/

	double zcold = model.gas_cooling_parameters.pre_enrich_z; /*cold gas metallicity*/

	if(y[1] > 0 && y[5] > 0) {
		zcold = y[5] / y[1];
	}

	double zhot = model.gas_cooling_parameters.pre_enrich_z; /*hot gas metallicity*/
	if(y[2] > 0 && y[6] > 0) {
		zhot = y[6] / y[2];
	}

	double rsub = 1.0- R;

	f[0] = SFR * rsub;
	f[1] = mcoolrate - (rsub + beta) * SFR;
	f[2] = - mcoolrate;
	f[3] = beta * SFR;
	f[4] = rsub * zcold * SFR;
	f[5] = mcoolrate * zhot + SFR * (yield - (rsub + beta) * zcold);
	f[6] = - mcoolrate * zhot;

	f[7] = beta * zcold * SFR;

	// Keeps track of total stellar mass formed.
	f[8] = SFR;


	return 0;
}

BasicPhysicalModel::BasicPhysicalModel(
		double ode_solver_precision,
		GasCooling gas_cooling,
		StellarFeedback stellar_feedback,
		StarFormation star_formation,
		RecyclingParameters recycling_parameters,
		GasCoolingParameters gas_cooling_parameters) :
	PhysicalModel(ode_solver_precision, basic_physicalmodel_evaluator, gas_cooling),
	stellar_feedback(stellar_feedback),
	star_formation(star_formation),
	recycling_parameters(recycling_parameters),
	gas_cooling_parameters(gas_cooling_parameters)
{
	// no-op
}

std::vector<double> BasicPhysicalModel::from_galaxy(const Subhalo &subhalo, const Galaxy &galaxy)
{

	std::vector<double> y(9);

	y[0] = galaxy.disk_stars.mass;
	y[1] = galaxy.disk_gas.mass;
	y[2] = subhalo.cold_halo_gas.mass; //This is the component that has the cooling gas.
	y[3] = subhalo.ejected_galaxy_gas.mass;
	y[4] = galaxy.disk_stars.mass_metals;
	y[5] = galaxy.disk_gas.mass_metals;
	y[6] = subhalo.cold_halo_gas.mass_metals;
	y[7] = subhalo.ejected_galaxy_gas.mass_metals;

	// Variable to keep track of total stellar mass formed.
	y[8] = 0;
	return y;
}

void BasicPhysicalModel::to_galaxy(const std::vector<double> &y, Subhalo &subhalo, Galaxy &galaxy, double delta_t)
{
	using namespace constants;

	/*In the case of starbursts one should be using the bulge instead of the disk
	 * properties.*/


	/* Check unrealistic cases*/
	if(y[0] < galaxy.disk_stars.mass){
		std::ostringstream os;
		os << "Galaxy decreased its stellar mass after disk star formation process.";
		throw invalid_argument(os.str());
	}

	galaxy.disk_stars.mass 					= y[0];
	galaxy.disk_gas.mass   					= y[1];
	subhalo.cold_halo_gas.mass 				= y[2];
	subhalo.ejected_galaxy_gas.mass 		= y[3];
	galaxy.disk_stars.mass_metals 			= y[4];
	galaxy.disk_gas.mass_metals 			= y[5];
	subhalo.cold_halo_gas.mass_metals 		= y[6];
	subhalo.ejected_galaxy_gas.mass_metals 	= y[7];

	// Calculate average SFR.
	galaxy.sfr_disk                         = y[8]/delta_t;

	/**
	 * Check that metallicities are not negative. If they are, mass in metals is set to zero.
	 */
	if(galaxy.disk_stars.mass_metals < tolerance){
		galaxy.disk_stars.mass_metals = 0;
	}
	if(galaxy.disk_gas.mass_metals < tolerance){
		galaxy.disk_gas.mass_metals = 0;
	}
	if(subhalo.cold_halo_gas.mass_metals < tolerance){
		subhalo.cold_halo_gas.mass_metals = 0;
	}
	if(subhalo.ejected_galaxy_gas.mass_metals < tolerance){
		subhalo.ejected_galaxy_gas.mass_metals = 0;
	}

	/**
	 * Check that masses are not negative. If they are, mass and metals are set to zero.
	 *
	 */
	if(galaxy.disk_stars.mass < tolerance){
		galaxy.disk_stars.mass = 0;
		galaxy.disk_stars.mass_metals = 0;
	}
	if(galaxy.disk_gas.mass < tolerance){
		galaxy.disk_gas.mass = 0;
		galaxy.disk_gas.mass_metals = 0;
	}
	if(subhalo.cold_halo_gas.mass < tolerance){
		subhalo.cold_halo_gas.mass = 0;
		subhalo.cold_halo_gas.mass_metals = 0;
	}
	if(subhalo.ejected_galaxy_gas.mass < tolerance){
		subhalo.ejected_galaxy_gas.mass = 0;
		subhalo.ejected_galaxy_gas.mass_metals = 0;
	}

}


std::vector<double> BasicPhysicalModel::from_galaxy_starburst(const Subhalo &subhalo, const Galaxy &galaxy)
{

	std::vector<double> y(9);

	y[0] = galaxy.bulge_stars.mass;
	y[1] = galaxy.bulge_gas.mass;
	y[2] = 0.0; //This is the component that has the cooling gas.
	y[3] = subhalo.ejected_galaxy_gas.mass;
	y[4] = galaxy.bulge_stars.mass_metals;
	y[5] = galaxy.bulge_gas.mass_metals;
	y[6] = 0.0; //This is the component that has the cooling gas mass in metals.
	y[7] = subhalo.ejected_galaxy_gas.mass_metals;

	// Variable to keep track of total stellar mass created.
	y[8] = 0;

	return y;
}

void BasicPhysicalModel::to_galaxy_starburst(const std::vector<double> &y, Subhalo &subhalo, Galaxy &galaxy, double delta_t)
{
	using namespace constants;

	/* Check unrealistic cases*/
	if(y[0] < galaxy.bulge_stars.mass){
		std::ostringstream os;
		os << "Galaxy decreased its stellar mass after burst of star formation.";
		throw invalid_argument(os.str());
	}

	/*In the case of starbursts one should be using the bulge instead of the disk
	 * properties.*/
	galaxy.bulge_stars.mass 				= y[0];
	galaxy.bulge_gas.mass   				= y[1];
	subhalo.ejected_galaxy_gas.mass 		= y[3];
	galaxy.bulge_stars.mass_metals 			= y[4];
	galaxy.bulge_gas.mass_metals 			= y[5];
	subhalo.ejected_galaxy_gas.mass_metals 	= y[7];

	// Calculate average SFR
	galaxy.sfr_bulge                        = y[8]/delta_t;

	/**
	 * Check that metallicities are not negative. If they are, mass in metals is set to zero.
	 */
	if(galaxy.bulge_stars.mass_metals < tolerance){
		galaxy.bulge_stars.mass_metals = 0;
	}
	if(galaxy.bulge_gas.mass_metals < tolerance){
		galaxy.bulge_gas.mass_metals = 0;
	}
	if(subhalo.ejected_galaxy_gas.mass_metals < tolerance){
		subhalo.ejected_galaxy_gas.mass_metals = 0;
	}

	/**
	 * Check that masses are not negative. If they are, mass and metals are set to zero.
	 *
	 */
	if(galaxy.bulge_stars.mass < tolerance){
		galaxy.bulge_stars.mass = 0;
		galaxy.bulge_stars.mass_metals = 0;
	}
	if(galaxy.bulge_gas.mass < tolerance){
		galaxy.bulge_gas.mass = 0;
		galaxy.bulge_gas.mass_metals = 0;
	}
	if(subhalo.ejected_galaxy_gas.mass < tolerance){
		subhalo.ejected_galaxy_gas.mass = 0;
		subhalo.ejected_galaxy_gas.mass_metals = 0;
	}

}

}  // namespace shark
