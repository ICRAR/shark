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
 * Physical model classes implementation
 */

#include <cmath>
#include <memory>
#include <vector>

#include "logging.h"
#include "numerical_constants.h"
#include "physical_model.h"

namespace shark {

static
int basic_physicalmodel_evaluator(double t, const double y[], double f[], void *data) {

	/** Functions describing the time derivatives of:
	 * mass exchange:
	 * f[0]: stellar mass of galaxy.
	 * f[1]: cold gas mass of galaxy.
	 * f[2]: cold gas in the halo (the one cooling).
	 * f[3]: hot gas mass.
	 * f[4]: ejected gas mass.
	 * f[5]: lost gas mass (only relevant during QSO outflows).
	 *
	 * metals exchange:
	 * f[6]: metals locked in the stellar mass of galaxies.
	 * f[7]: metals locked in the cold gas mass of galaxies.
	 * f[8]: metals locked in the cold gas mass of the halo.
	 * f[9]: metals locked in the hot halo gas reservoir.
	 * f[10]: metals locked in the ejected gas mass reservoir.
	 * f[11]: metals locked in the lost gas mass reservoir (only relevant for QSO feedback).
	 *
	 * tracking total mass and angular momentum formed.
	 * f[12]: total stellar mass formed (without recycling included).
	 * f[13]: total stellar angular momentum of the galaxy component.
	 *
	 * angular momentum exchange:
	 * f[14]: stellar angular momentum of the galaxy component.
	 * f[15]: gas angular momentum of the galaxy component.
	 * f[16]: angular momentum of the cold halo gas component.
	 * f[17]: angular momentum of the hot gas component.
	 * f[18]: angular momentum of the ejected gas component.
	 */

	auto params = static_cast<BasicPhysicalModel::solver_params *>(data);
	auto &model = static_cast<BasicPhysicalModel &>(params->model);

	double R = model.recycling_parameters.recycle; /*recycling fraction of newly formed stars*/

	double yield = model.recycling_parameters.yield; /*yield of newly formed stars*/

	double mcoolrate = params->mcoolrate; /*cooling rate in units of Msun/Gyr*/

	// Define minimum gas metallicities.
	double zcold = model.gas_cooling_parameters.pre_enrich_z; /*cold gas minimum metallicity*/
	double zhot = model.gas_cooling_parameters.pre_enrich_z; /*cooling gas minimum metallicity*/

	// Define angular momentum parameters.
	double jgas = 2.0 * params->vgal * params->rgas / constants::RDISK_HALF_SCALE; /*current sAM of the cold gas*/
	double jrate = 0; /*variable that saves the angular momentum transfer rate from gas to stars*/

	// Define current gas metallicity and angular momentum.
	if(y[1] > 0 && y[7] > 0) {
		zcold = y[7] / y[1];
		jgas  = y[15] / y[1];
	}

	// Define metallicity of cooling gas.
	if(params->zcool > zhot){
		zhot = params->zcool;
	}

	// compute modified yield if user sets that option
	double yield_eff = 0.0;
	if(model.recycling_parameters.evolving_yield){
		yield_eff  = yield - zcold*0.25; //from Robotham et al. (2019)
	}
	else{
		yield_eff  = yield;
	}

	// Calculate SFR.
	double SFR   = model.star_formation.star_formation_rate(y[1], y[0], params->rgas, params->rstar, zcold, params->redshift, params->burst, params->vgal, jrate, jgas);

	// Initialize mass loading and angular momentum loading parameters related to star formation.
	double beta1 = 0, beta2 = 0;
	double betaj_1 = 0, betaj_2 = 0;

	// Initialize mass loading related to QSO outflows;
	double beta_qso1 = 0, beta_qso2 = 0;

	// Calculate mass and angular momentum loading from stellar feedback process.
	model.stellar_feedback.outflow_rate(SFR, params->vsubh, params->vgal, params->redshift, beta1, beta2, betaj_1, betaj_2); /*mass loading parameter*/

	// Calculate the mass and metal loading from QSO feedback.
	model.agn_feedback.qso_outflow_rate(y[1], params->smbh, zcold, params->vgal, SFR, y[0]+y[1], params->rstar, beta_qso1, beta_qso2);

	// Retained fraction.
	double rsub = 1.0-R;

	// Mass transfer equations.
	f[0] = SFR * rsub;
	f[1] = mcoolrate - (rsub + beta1 + beta_qso1) * SFR;
	f[2] = - mcoolrate;
	f[3] = (beta1 + beta_qso1 - (beta2 + beta_qso2)) * SFR;
	f[4] = beta2 * SFR;
	f[5] = beta_qso2 * SFR;

	// Metallicity transfer equations.
	f[6] = rsub * zcold * SFR;
	f[7] = mcoolrate * zhot + SFR * (yield_eff - (rsub + beta1 + beta_qso1) * zcold);
	f[8] = - mcoolrate * zhot;
	f[9] = ((beta1 + beta_qso1) - (beta2 + beta_qso2) ) * SFR * zcold;
	f[10] = beta2 * zcold * SFR;
	f[11] = beta_qso2  * zcold * SFR;

	// Keeps track of total stellar mass formed and the metals locked up in it..
	f[12] = SFR;
	f[13] = zcold * SFR;

	// Solve angular momentum equations.
	f[14] = rsub * jrate;
	f[15] = mcoolrate * params->jcold_halo - (rsub + betaj_1) * jrate;
	f[16] = - mcoolrate * params->jcold_halo;
	f[17] = (betaj_1 - betaj_2) * jrate;
	f[18] = betaj_2 * jrate;

	return 0;
}

BasicPhysicalModel::BasicPhysicalModel(
		double ode_solver_precision,
		GasCooling gas_cooling,
		StellarFeedback stellar_feedback,
		StarFormation star_formation,
		AGNFeedback agn_feedback,
		RecyclingParameters recycling_parameters,
		GasCoolingParameters gas_cooling_parameters,
		AGNFeedbackParameters agn_parameters) :
	PhysicalModel(ode_solver_precision, basic_physicalmodel_evaluator, std::move(gas_cooling)),
	stellar_feedback(stellar_feedback),
	star_formation(std::move(star_formation)),
	agn_feedback(std::move(agn_feedback)),
	recycling_parameters(recycling_parameters),
	gas_cooling_parameters(std::move(gas_cooling_parameters)),
	agn_parameters(agn_parameters)
{
	// no-op
}

void BasicPhysicalModel::from_galaxy(std::vector<double> &y, const Subhalo &subhalo, const Galaxy &galaxy)
{

	/** Variables introduced to solve ODE equations.
	 * y[0]: stellar mass of galaxy.
	 * y[1]: cold gas mass of galaxy.
	 * y[2]: cold gas in the halo (the one cooling).
	 * y[3]: hot gas mass;
	 * y[4]: ejected gas mass;
	 * y[5]: lost gas mass.
	 *
	 * y[6]: metals locked in the stellar mass of galaxies.
	 * y[7]: metals locked in the cold gas mass of galaxies.
	 * y[8]: metals locked in the cold gas mass of the halo.
	 * y[9]: metals locked in the hot halo gas reservoir.
	 * y[10]: metals locked in the ejected gas mass reservoir.
	 * y[11]: metals locked in the lost has mass reservoir.
	 *
	 * y[12]: total stellar mass formed (without recycling included).
	 * y[13]: total stellar mass in metals formed (without recycling included).
	 *
	 * Equations dealing with angular momentum:
	 * y[14]: total stellar angular momentum of the bulge.
	 * y[15]: total gas angular momentum of the bulge.
	 * y[16]: total cold gas angular momentum of the cold halo gas component.
	 * y[17]: total angular momentum of the hot gas component.
	 * y[18]: total angular momentum of the ejected gas component.
	 *
	 */

	// Define mass inputs.
	y[0] = galaxy.disk_stars.mass;
	y[1] = galaxy.disk_gas.mass;
	y[2] = subhalo.cold_halo_gas.mass; //This is the component that has the cooling gas.
	y[3] = subhalo.hot_halo_gas.mass;
	y[4] = subhalo.ejected_galaxy_gas.mass;
	y[5] = subhalo.lost_galaxy_gas.mass;

	// Define mass in metals inputs.
	y[6] = galaxy.disk_stars.mass_metals;
	y[7] = galaxy.disk_gas.mass_metals;
	y[8] = subhalo.cold_halo_gas.mass_metals;
	y[9] = subhalo.hot_halo_gas.mass_metals;
	y[10] = subhalo.ejected_galaxy_gas.mass_metals;
	y[11] = subhalo.lost_galaxy_gas.mass_metals;

	// Variable to keep track of total stellar mass and metals formed in this SF episode.
	y[12] = 0;
	y[13] = 0;

	// Equations of angular momentum exchange. Input total angular momentum.
	y[14] = galaxy.disk_stars.sAM * galaxy.disk_stars.mass;
	y[15] = galaxy.disk_gas.sAM * galaxy.disk_gas.mass;
	y[16] = subhalo.cold_halo_gas.sAM * subhalo.cold_halo_gas.mass;
	y[17] = subhalo.hot_halo_gas.sAM * subhalo.hot_halo_gas.mass;
	y[18] = subhalo.ejected_galaxy_gas.sAM * subhalo.ejected_galaxy_gas.mass;
}

void BasicPhysicalModel::to_galaxy(const std::vector<double> &y, Subhalo &subhalo, Galaxy &galaxy, double delta_t)
{
	using namespace constants;

	/* Check unrealistic cases*/
	if(y[0] < galaxy.disk_stars.mass){
		std::ostringstream os;
		os << "Galaxy decreased its stellar mass after disk star formation process.";
		throw invalid_argument(os.str());
	}

	// Assign new masses.
	galaxy.disk_stars.mass 			= y[0];
	galaxy.disk_gas.mass   			= y[1];
	subhalo.cold_halo_gas.mass 		= y[2];
	subhalo.hot_halo_gas.mass               = y[3];
	subhalo.ejected_galaxy_gas.mass 	= y[4];
	subhalo.lost_galaxy_gas.mass		= y[5];

	// Assign new mass in metals.
	galaxy.disk_stars.mass_metals 			= y[6];
	galaxy.disk_gas.mass_metals 			= y[7];
	subhalo.cold_halo_gas.mass_metals 		= y[8];
	subhalo.hot_halo_gas.mass_metals                = y[9];
	subhalo.ejected_galaxy_gas.mass_metals 	        = y[10];
	subhalo.lost_galaxy_gas.mass_metals		= y[11];

	// Calculate average SFR and metallicity of newly formed stars.
	galaxy.sfr_disk                         += y[12]/delta_t;
	galaxy.sfr_z_disk                       += y[13]/delta_t;

	// Equations of angular momentum exchange. Input total angular momentum.
	// Redefine angular momentum ONLY if the new value is > 0.
	if(y[15] > 0){

		// Assign new specific angular momenta.
		galaxy.disk_gas.sAM            = y[15] / galaxy.disk_gas.mass;
		subhalo.cold_halo_gas.sAM      = y[16] / subhalo.cold_halo_gas.mass;
		subhalo.hot_halo_gas.sAM       = y[17] / subhalo.hot_halo_gas.mass;
		subhalo.ejected_galaxy_gas.sAM = y[18] / subhalo.ejected_galaxy_gas.mass;

		// Assign new sizes based on new AM.
		galaxy.disk_gas.rscale   = galaxy.disk_gas.sAM   / galaxy.vmax * constants::EAGLEJconv;

		// check for unrealistic cases.
		if (std::isnan(galaxy.disk_gas.sAM) || std::isnan(galaxy.disk_gas.rscale)) {
			throw invalid_argument("rgas or sAM are NaN, cannot continue at physical model");
		}

	}
	if (y[14] > 0){

		// Assign new specific angular momenta.
		galaxy.disk_stars.sAM          = y[14] / galaxy.disk_stars.mass;

		// Assign new sizes based on new AM.
		galaxy.disk_stars.rscale = galaxy.disk_stars.sAM / galaxy.vmax * constants::EAGLEJconv;

		// check for unrealistic cases.
		if(galaxy.disk_stars.rscale <= constants::tolerance && galaxy.disk_stars.mass > 0){
			std::ostringstream os;
			os << "Galaxy with extremely small size, rdisk_stars < 1e-10, in physical model";
			throw invalid_argument(os.str());
		}
	}

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
	if(subhalo.hot_halo_gas.mass_metals < tolerance){
		subhalo.hot_halo_gas.mass_metals = 0;
	}
	if(subhalo.ejected_galaxy_gas.mass_metals < tolerance){
		subhalo.ejected_galaxy_gas.mass_metals = 0;
	}
	if(subhalo.lost_galaxy_gas.mass_metals < tolerance){
		subhalo.lost_galaxy_gas.mass_metals = 0;
	}

	/**
	 * Check that metallicities are not bigger than the gas mass. If they are, mass in metals to gas mass.
	 */
	if(galaxy.disk_gas.mass_metals > galaxy.disk_gas.mass){
		galaxy.disk_gas.mass_metals = galaxy.disk_gas.mass;
	}
	if(subhalo.cold_halo_gas.mass_metals > subhalo.cold_halo_gas.mass){
		subhalo.cold_halo_gas.mass_metals = subhalo.cold_halo_gas.mass;
	}
	if(subhalo.hot_halo_gas.mass_metals > subhalo.hot_halo_gas.mass){
		subhalo.hot_halo_gas.mass_metals = subhalo.hot_halo_gas.mass;
	}
	if(subhalo.ejected_galaxy_gas.mass_metals > subhalo.ejected_galaxy_gas.mass){
		subhalo.ejected_galaxy_gas.mass_metals = subhalo.ejected_galaxy_gas.mass;
	}
	if(subhalo.lost_galaxy_gas.mass_metals > subhalo.lost_galaxy_gas.mass){
		subhalo.lost_galaxy_gas.mass_metals = subhalo.lost_galaxy_gas.mass;
	}

	/**
	 * Check that masses are not negative. If they are, mass and metals are set to zero.
	 *
	 */
	if(galaxy.disk_stars.mass < tolerance){
		galaxy.disk_stars.restore_baryon();
	}
	if(galaxy.disk_gas.mass < tolerance){
		galaxy.disk_gas.restore_baryon();
	}
	if(subhalo.cold_halo_gas.mass < tolerance){
		subhalo.cold_halo_gas.restore_baryon();
	}
	if(subhalo.hot_halo_gas.mass < tolerance){
		subhalo.hot_halo_gas.restore_baryon();
	}
	if(subhalo.ejected_galaxy_gas.mass < tolerance){
		subhalo.ejected_galaxy_gas.restore_baryon();
	}
	if(subhalo.lost_galaxy_gas.mass < tolerance){
		subhalo.lost_galaxy_gas.restore_baryon();
	}

	/* Check unrealistic cases*/
	if(galaxy.disk_gas.mass < galaxy.disk_gas.mass_metals || subhalo.hot_halo_gas.mass < subhalo.hot_halo_gas.mass_metals || subhalo.ejected_galaxy_gas.mass < subhalo.ejected_galaxy_gas.mass_metals){
		std::ostringstream os;
		os << "Galaxy has more gas mass in metals that total gas mass.";
		throw invalid_argument(os.str());
	}
}


void BasicPhysicalModel::from_galaxy_starburst(std::vector<double> &y, const Subhalo &subhalo, const Galaxy &galaxy)
{
	/** Variables introduced to solve ODE equations.
	 * y[0]: stellar mass of galaxy.
	 * y[1]: cold gas mass of galaxy.
	 * y[2]: cold gas in the halo (the one cooling).
	 * y[3]: hot gas mass;
	 * y[4]: ejected gas mass;
	 * y[5]: lost gas mass.
	 *
	 * y[6]: metals locked in the stellar mass of galaxies.
	 * y[7]: metals locked in the cold gas mass of galaxies.
	 * y[8]: metals locked in the cold gas mass of the halo.
	 * y[9]: metals locked in the hot halo gas reservoir.
	 * y[10]: metals locked in the ejected gas mass reservoir.
	 * y[11]: metals locked in the lost has mass reservoir.
	 *
	 * y[12]: total stellar mass formed (without recycling included).
	 * y[13]: total stellar mass in metals formed (without recycling included).
	 *
	 * Equations dealing with angular momentum:
	 * y[14]: total stellar angular momentum of the bulge.
	 * y[15]: total gas angular momentum of the bulge.
	 * y[16]: total cold gas angular momentum of the cold halo gas component.
	 * y[17]: total angular momentum of the hot gas component.
	 * y[18]: total angular momentum of the ejected gas component.
	 *
	 */

	// Define mass inputs.
	y[0] = galaxy.bulge_stars.mass;
	y[1] = galaxy.bulge_gas.mass;
	y[2] = 0; //there is no gas cooling.
	y[3] = subhalo.hot_halo_gas.mass;
	y[4] = subhalo.ejected_galaxy_gas.mass;
	y[5] = subhalo.lost_galaxy_gas.mass;

	// Define mass in metals inputs.
	y[6] = galaxy.bulge_stars.mass_metals;
	y[7] = galaxy.bulge_gas.mass_metals;
	y[8] = 0; //there is no gas cooling.
	y[9] = subhalo.hot_halo_gas.mass_metals;
	y[10] = subhalo.ejected_galaxy_gas.mass_metals;
	y[11] = subhalo.lost_galaxy_gas.mass_metals;

	// Variable to keep track of total stellar mass and metals formed in this SF episode.
	y[12] = 0;
	y[13] = 0;

	// Equations of angular momentum exchange are ignored in the case of starbursts.
}

void BasicPhysicalModel::to_galaxy_starburst(const std::vector<double> &y, Subhalo &subhalo, Galaxy &galaxy, double delta_t, bool from_galaxy_merger)
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

	// Accumulated burst stellar mass in the corresponding baryon budget depending on triggering mechanism:
	if(from_galaxy_merger){
		galaxy.galaxymergers_burst_stars.mass                 += y[0] -  galaxy.bulge_stars.mass;
		galaxy.galaxymergers_burst_stars.mass_metals          += y[6] -  galaxy.bulge_stars.mass_metals;
		// Calculate average SFR and metallicity of newly formed stars.
		galaxy.sfr_bulge_mergers                              += y[12]/delta_t;
		galaxy.sfr_z_bulge_mergers                            += y[13]/delta_t;
	}
	else{
		galaxy.diskinstabilities_burst_stars.mass             += y[0] -  galaxy.bulge_stars.mass;
		galaxy.diskinstabilities_burst_stars.mass_metals      += y[6] -  galaxy.bulge_stars.mass_metals;
		// Calculate average SFR and metallicity of newly formed stars.
		galaxy.sfr_bulge_diskins                              += y[12]/delta_t;
		galaxy.sfr_z_bulge_diskins                            += y[13]/delta_t;
	}

	// Assign new masses.
	galaxy.bulge_stars.mass 		= y[0];
	galaxy.bulge_gas.mass   		= y[1];
	subhalo.hot_halo_gas.mass               = y[3];
	subhalo.ejected_galaxy_gas.mass 	= y[4];
	subhalo.lost_galaxy_gas.mass 		= y[5];

	// Assign new mass in metals.
	galaxy.bulge_stars.mass_metals 		= y[6];
	galaxy.bulge_gas.mass_metals 		= y[7];
	subhalo.hot_halo_gas.mass_metals        = y[9];
	subhalo.ejected_galaxy_gas.mass_metals 	= y[10];
	subhalo.lost_galaxy_gas.mass_metals	= y[11];

	// Equations of angular momentum exchange are ignored in the case of starbursts.

	/**
	 * Check that metallicities are not negative. If they are, mass in metals is set to zero.
	 */
	if(galaxy.bulge_stars.mass_metals < tolerance){
		galaxy.bulge_stars.mass_metals = 0;
	}
	if(galaxy.bulge_gas.mass_metals < tolerance){
		galaxy.bulge_gas.mass_metals = 0;
	}
	if(subhalo.hot_halo_gas.mass_metals < tolerance){
		subhalo.hot_halo_gas.mass_metals = 0;
	}
	if(subhalo.ejected_galaxy_gas.mass_metals < tolerance){
		subhalo.ejected_galaxy_gas.mass_metals = 0;
	}
	if(subhalo.lost_galaxy_gas.mass_metals < tolerance){
		subhalo.lost_galaxy_gas.mass_metals = 0;
	}

	/**
	 * Check that metallicities are not bigger than the gas mass. If they are, mass in metals to gas mass.
	 */
	if(galaxy.bulge_gas.mass_metals > galaxy.bulge_gas.mass){
		galaxy.bulge_gas.mass_metals = galaxy.bulge_gas.mass;
	}
	if(subhalo.hot_halo_gas.mass_metals > subhalo.hot_halo_gas.mass){
		subhalo.hot_halo_gas.mass_metals = subhalo.hot_halo_gas.mass;
	}
	if(subhalo.ejected_galaxy_gas.mass_metals > subhalo.ejected_galaxy_gas.mass){
		subhalo.ejected_galaxy_gas.mass_metals = subhalo.ejected_galaxy_gas.mass;
	}
	if(subhalo.lost_galaxy_gas.mass_metals > subhalo.lost_galaxy_gas.mass){
		subhalo.lost_galaxy_gas.mass_metals = subhalo.lost_galaxy_gas.mass;
	}

	/**
	 * Check that masses are not negative. If they are, mass and metals are set to zero.
	 *
	 */
	if(galaxy.bulge_stars.mass < tolerance){
		galaxy.bulge_stars.restore_baryon();
	}
	if(galaxy.bulge_gas.mass < tolerance){
		galaxy.bulge_gas.restore_baryon();
	}
	if(subhalo.hot_halo_gas.mass < tolerance){
		subhalo.hot_halo_gas.restore_baryon();
	}
	if(subhalo.ejected_galaxy_gas.mass < tolerance){
		subhalo.ejected_galaxy_gas.restore_baryon();
	}
	if(subhalo.lost_galaxy_gas.mass < tolerance){
		subhalo.lost_galaxy_gas.restore_baryon();
	}

	/* Check unrealistic cases*/
	if(galaxy.bulge_gas.mass < galaxy.bulge_gas.mass_metals || subhalo.hot_halo_gas.mass < subhalo.hot_halo_gas.mass_metals || subhalo.ejected_galaxy_gas.mass < subhalo.ejected_galaxy_gas.mass_metals){
		std::ostringstream os;
		os << "Galaxy has more gas mass in metals that total gas mass.";
		throw invalid_argument(os.str());
	}

}

}  // namespace shark
