//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2018
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


#include <cmath>
#include <gsl/gsl_errno.h>


#include "environment.h"
#include "subhalo.h"
#include "galaxy.h"
#include "cosmology.h"
#include "dark_matter_halos.h"
#include "simulation.h"

namespace shark {

struct galaxy_properties_for_root_solver {
	double z;
	SubhaloPtr primary;
	Subhalo &secondary;
};

EnvironmentParameters::EnvironmentParameters(const Options &options)
{
	options.load("environment.stripping", stripping, true);
	options.load("environment.gradual_stripping_halo", gradual_stripping_halo);
	options.load("environment.gradual_stripping_ism", gradual_stripping_ism);
	options.load("environment.tidal_stripping", tidal_stripping);
	options.load("environment.minimum_halo_mass_fraction", minimum_halo_mass_fraction);
}

Environment::Environment(const EnvironmentParameters &parameters,
		DarkMatterHalosPtr darkmatterhalos,
		CosmologyPtr cosmology,
		CosmologicalParameters cosmo_params,
		SimulationParameters simparams
		):
	parameters(parameters),
	darkmatterhalos(darkmatterhalos),
	cosmology(cosmology),
	cosmo_params(cosmo_params),
	simparams(simparams),
	root_solver(1000)
	{
	//no-opt
}

void Environment::process_satellite_subhalo_environment(Subhalo &satellite_subhalo, SubhaloPtr &central_subhalo, double z){

	// Ejected gas is moved to the budget of ejected gas of the central, as this gas escaped
	// the subhalo of the satellite. This is always the case.
	//central_subhalo.ejected_galaxy_gas += satellite_subhalo.ejected_galaxy_gas;
	central_subhalo->lost_galaxy_gas += satellite_subhalo.lost_galaxy_gas + satellite_subhalo.ejected_galaxy_gas;
	central_subhalo->stellar_halo += satellite_subhalo.stellar_halo;
	central_subhalo->mean_galaxy_making_stellar_halo += satellite_subhalo.mean_galaxy_making_stellar_halo;

	satellite_subhalo.ejected_galaxy_gas.restore_baryon();
	satellite_subhalo.lost_galaxy_gas.restore_baryon();
	satellite_subhalo.stellar_halo.restore_baryon();
	satellite_subhalo.mean_galaxy_making_stellar_halo = 0;

	// Assume halo gas of subhalos that will disappear in the next snapshot is fully stripped.
	if(satellite_subhalo.infall_t == z){
		satellite_subhalo.transfer_halo_gas_to(central_subhalo);
	}
	else{
		// Treatment for all satellite subhalos.
		// Remove hot halo gas only if stripping is applied.
		if(parameters.stripping){

			if(satellite_subhalo.hot_halo_gas.mass > 0 || satellite_subhalo.cold_halo_gas.mass > 0){
				if(parameters.gradual_stripping_halo){
					//first check whether the function is positive at Rvir_infall. In that case, the satellite subhalo experiences no stripping:
					//second, check whether function is negative at Rvir_infall/100. In that case assume all hot gas is stripped.
					double mhot_removed = 0;
					double r_rps = 0;

					auto func_rvir = ram_pressure_stripping_hot_gas(central_subhalo, satellite_subhalo, satellite_subhalo.rvir_infall, z);
					auto func_rvirdiv100 = ram_pressure_stripping_hot_gas(central_subhalo, satellite_subhalo, satellite_subhalo.rvir_infall/100, z);

					if(func_rvir > 0){
						r_rps = satellite_subhalo.hot_halo_gas_r_rps;
					}
					else if (func_rvir < 0 && func_rvirdiv100 > 0){
						r_rps = process_ram_pressure_stripping(central_subhalo, satellite_subhalo, z);
					}
					else if (func_rvir < 0 && func_rvirdiv100 < 0){
						r_rps = 0;
					}

					auto mhot_tot = satellite_subhalo.cold_halo_gas.mass + satellite_subhalo.hot_halo_gas.mass + satellite_subhalo.hot_halo_gas_stripped.mass;

					// If the ram-pressure stripping radius has decreased from previous timesteps, then compute how much new gas is lost.
					if(r_rps < satellite_subhalo.hot_halo_gas_r_rps && r_rps > 0){
						// 1. compute hot gas outside r_rps
						// 2. update satellite subhalo r_rps
						// 3. update hot gas that has been stripped.
						// 4. move that extra hot gas to central subhalo's reservoir.

						mhot_removed = mhot_tot * (1 - std::pow(r_rps/satellite_subhalo.rvir_infall,2)) - satellite_subhalo.hot_halo_gas_stripped.mass;
						if(mhot_removed < 0){
							mhot_removed = 0;
						}

						if(mhot_removed > satellite_subhalo.hot_halo_gas.mass + satellite_subhalo.cold_halo_gas.mass) {
							mhot_removed = satellite_subhalo.hot_halo_gas.mass + satellite_subhalo.cold_halo_gas.mass;
						}


						// remove halo gas in the proportion of the cold-to-hot halo mass ratio
						auto frac_gas_phase = satellite_subhalo.hot_halo_gas.mass / (satellite_subhalo.hot_halo_gas.mass + satellite_subhalo.cold_halo_gas.mass);
						float metals_removed_hot = 0;
						float metals_removed_cold = 0;

						if(satellite_subhalo.cold_halo_gas.mass > 0){
							auto zcold = satellite_subhalo.cold_halo_gas.mass_metals / satellite_subhalo.cold_halo_gas.mass;
							metals_removed_cold = zcold * mhot_removed * (1 - frac_gas_phase);

							satellite_subhalo.cold_halo_gas.mass -= mhot_removed * (1 - frac_gas_phase);
							satellite_subhalo.cold_halo_gas.mass_metals -= metals_removed_cold;

							if(satellite_subhalo.cold_halo_gas.mass < 0){
								satellite_subhalo.cold_halo_gas.restore_baryon();
							}

						}
						auto zhot = satellite_subhalo.hot_halo_gas.mass_metals / satellite_subhalo.hot_halo_gas.mass;
						metals_removed_hot = zhot * mhot_removed * frac_gas_phase;

						// track removed hot halo gas and ram-pressure stripping radius.
						satellite_subhalo.hot_halo_gas_r_rps = r_rps;
						satellite_subhalo.hot_halo_gas_stripped.mass += mhot_removed;
						satellite_subhalo.hot_halo_gas_stripped.mass_metals += metals_removed_cold + metals_removed_hot;

						// remove hot gas from satellite subhalo
						satellite_subhalo.hot_halo_gas.mass -= mhot_removed * frac_gas_phase;
						satellite_subhalo.hot_halo_gas.mass_metals -= metals_removed_hot;

						// check mass is larger than 0. If not, then restore baryon.
						if(satellite_subhalo.hot_halo_gas.mass < 0){
							satellite_subhalo.hot_halo_gas.restore_baryon();
						}

						central_subhalo->hot_halo_gas.mass += mhot_removed;
						central_subhalo->hot_halo_gas.mass_metals += metals_removed_cold + metals_removed_hot;

					}
					else if(r_rps == 0){
						satellite_subhalo.transfer_halo_gas_to(central_subhalo);
					}
				}
				else{
					satellite_subhalo.transfer_halo_gas_to(central_subhalo);
				}
			}

		}
	}

	// Remove part of the stellar content of satellite if relevant.
	// Cases to consider:
	// 1. A satellite subhalo that has a central galaxy type=1
	// 2. A satellite subhalo that does not have a central. In this case we have to check the type II satellites
	// in the halo of the central to see if they have been stripped accordingly.
	if(parameters.tidal_stripping){

		// mass in metals to be tidally stripped is computed inside the function remove_tidal_stripped_stars.
		BaryonBase lost_stellar;

		if(satellite_subhalo.type1_galaxy()){
			// Apply here model of Errani et al. (2015) with rstar/a=0.2.
			float ratio_mass = satellite_subhalo.Mvir / satellite_subhalo.Mvir_infall;
			// Apply a maximum stripping of 99% in halo mass and on the other and a maximum of 1 in the ratio..
			if(ratio_mass < parameters.minimum_halo_mass_fraction){
				ratio_mass = parameters.minimum_halo_mass_fraction;
			}
			if(ratio_mass > 1){
				ratio_mass = 1;
			}

			// compute how much has been lost since galaxy infall
			float ratio_sm = std::pow(2, 3.57) * std::pow(ratio_mass, 2.06) / std::pow( 1 + ratio_mass, 3.57);
			lost_stellar.mass = (1 - ratio_sm) * satellite_subhalo.star_central_infall.mass;

			if(lost_stellar.mass < 0){
				lost_stellar.restore_baryon();
			}

			// now remove what has already been lost in previous snapshots to compute what should be subtracted from
			// the central galaxy in the satellite subhalo now.
			lost_stellar.mass -= satellite_subhalo.type1_galaxy()->stars_tidal_stripped.mass;

			if(lost_stellar.mass > 0){
				lost_stellar = remove_tidal_stripped_stars(central_subhalo, *satellite_subhalo.type1_galaxy(), lost_stellar);

				// add the stripped material to central subhalo
				central_subhalo->stellar_halo += lost_stellar;
			}
		}
		// check all type 2 galaxies in the central subhalo and strip only those that haven't been stripped yet.
		if(central_subhalo->type2_galaxies_count() > 0){
			float ratio_sm = std::pow(2.0, 3.57) * std::pow(parameters.minimum_halo_mass_fraction, 2.06) / std::pow( 1 + parameters.minimum_halo_mass_fraction, 3.57);
			// in this case loop over type II satellites and if they have not been stripped
			// then strip a fraction of the satellite mass.
			for (auto &satellite: central_subhalo->type2_galaxies() ){
				// If no stripping has occurred in this type II then assume 99% loss of DM.
				// is stars_tidal_stripped.mass>0 is because stripping should have occurred already.
				if(satellite.stars_tidal_stripped.mass == 0){
					// compute how much has been lost since galaxy infall
					lost_stellar.mass = (1 - ratio_sm) * satellite.stellar_mass();

					lost_stellar = remove_tidal_stripped_stars(central_subhalo, satellite, lost_stellar);

					// add the stripped material to central subhalo
					central_subhalo->stellar_halo += lost_stellar;
				}
			}
		}

		// consistency check:
		BaryonBase total_mshalo;
		for (auto &gal: central_subhalo->galaxies){
			total_mshalo += gal.stars_tidal_stripped;
		}
	}

}

BaryonBase Environment::remove_tidal_stripped_stars(SubhaloPtr &subhalo, Galaxy &galaxy, BaryonBase lost_stellar){


	if(lost_stellar.mass > 0){

		BaryonBase lost_cold_gas;

		//check that lost mass does not exceed the total stellar mass of the galaxy to be stripped.
		if(lost_stellar.mass > galaxy.stellar_mass()){
			lost_stellar.mass = galaxy.stellar_mass();
			lost_stellar.mass_metals = galaxy.stellar_mass_metals();
			subhalo->mean_galaxy_making_stellar_halo += galaxy.stellar_mass() * lost_stellar.mass;

			galaxy.disk_stars.restore_baryon();
			galaxy.bulge_stars.restore_baryon();
		}
		else{
			subhalo->mean_galaxy_making_stellar_halo += galaxy.stellar_mass() * lost_stellar.mass;

			// strip first the disk of the galaxy and then the bulge:
			if(lost_stellar.mass < galaxy.disk_stars.mass){
				// in this case we strip material from the disk but not the bulge
				float frac_lost = lost_stellar.mass/galaxy.disk_stars.mass;
				lost_stellar.mass_metals = frac_lost * galaxy.disk_stars.mass_metals;

				// compute how much is lost of the cold gas reservoir
				lost_cold_gas.mass = frac_lost * galaxy.disk_gas.mass;
				lost_cold_gas.mass_metals = frac_lost * galaxy.disk_gas.mass_metals;

				// add lost cold gas reservoir to hot halo and remove it from disk.
				subhalo->hot_halo_gas += lost_cold_gas;
				galaxy.disk_stars -= lost_stellar;
				galaxy.disk_gas -= lost_cold_gas;
			}
			else{
				// in this case we strip all the disk and remove a fraction of the bulge.
				BaryonBase lost_bulge;

				// compute mass lost in bulge
				lost_bulge.mass = lost_stellar.mass - galaxy.disk_stars.mass;
				float frac_lost = lost_bulge.mass / galaxy.bulge_stars.mass;
				lost_bulge.mass_metals = frac_lost * galaxy.bulge_stars.mass_metals;
				lost_stellar.mass_metals = galaxy.disk_stars.mass_metals + lost_bulge.mass_metals;

				//compute loss of gas in the bulge;
				lost_cold_gas.mass = frac_lost * galaxy.bulge_gas.mass;
				lost_cold_gas.mass_metals = frac_lost * galaxy.bulge_gas.mass_metals;

				// compute mass lost to the hot gas
				subhalo->hot_halo_gas += (galaxy.disk_gas + lost_cold_gas);
				galaxy.disk_stars.restore_baryon();
				galaxy.disk_gas.restore_baryon();
				//remove mass from bulge
				galaxy.bulge_stars -= lost_bulge;
				galaxy.bulge_gas -= lost_cold_gas;
			}

			// sanity checks
			if(galaxy.disk_stars.mass < 0){
				galaxy.disk_stars.restore_baryon();
			}
			if(galaxy.bulge_stars.mass < 0){
				galaxy.bulge_stars.restore_baryon();
			}
			if(galaxy.disk_gas.mass < 0){
				galaxy.disk_gas.restore_baryon();
			}
			if(galaxy.bulge_gas.mass < 0){
				galaxy.bulge_gas.restore_baryon();
			}
		}

		// accumulate what has been lost to tidal stripping in this galaxy.
		galaxy.stars_tidal_stripped += lost_stellar;

	}

	return lost_stellar;

}

double Environment::process_ram_pressure_stripping(const SubhaloPtr &primary,
		Subhalo &secondary,
		double z){

	galaxy_properties_for_root_solver props = {
		z,
		primary,
		secondary
	};

	struct EnvironmentProcessAndProps {
		Environment *environment;
		galaxy_properties_for_root_solver *props;
	};

	auto f = [](double r, void *ctx) -> double {
		auto *env_and_props = static_cast<EnvironmentProcessAndProps *>(ctx);
		return env_and_props->environment->ram_pressure_stripping_hot_gas(env_and_props->props->primary, env_and_props->props->secondary, r, env_and_props->props->z);
	};

	EnvironmentProcessAndProps env_and_props = {this, &props};
	double result = 0;
	try{
		result = root_solver.root_solver_function(f, &env_and_props, env_and_props.props->secondary.rvir_infall/100.0, env_and_props.props->secondary.rvir_infall, 0, parameters.Accuracy_RPS);
	} catch (gsl_error &e) {
		auto gsl_errno = e.get_gsl_errno();
		std::ostringstream os;
		os << "RPS of hot gas failed with GSL error number " << gsl_errno << ": ";
		os << gsl_strerror(gsl_errno) << ", reason=" << e.get_reason();
		LOG(warning) << os.str();
	}

	// Avoid negative values.
	if(result < 0){
		result = 0.0;
	}

	result = cosmology->physical_to_comoving_size(result, z);

	return result;


}
double Environment::ram_pressure_stripping_hot_gas(const SubhaloPtr &primary,
		Subhalo &secondary,
		double r,
		double z){
	/* here we evaluate the function that needs to go to zero to find the best R_RPS solution for the hot gas
	this comes from the following equation:
	rho_hot_cen(Rsat) * vsat^2.0 = G * Msat(<r)*mhot_sat / (8 * Rvir_sat * r^3)
	Here, Rsat is the distance from the halo centre to the satellite galaxy; vsat is the relative velocity of the subhalo and r what we need to solve for.
	*/

	// Find the objects physical position
	double conversion_factor = cosmo_params.Hubble_h * (1 +  z);
	double xrel = (secondary.position.x - primary->position.x) / conversion_factor;
	double yrel = (secondary.position.y - primary->position.y) / conversion_factor;
	double zrel = (secondary.position.z - primary->position.z) / conversion_factor;
	double rsat = std::sqrt(xrel*xrel + yrel*yrel + zrel*zrel);

	// Find the physical velocity + the hubble flow
	double hubble_flow = rsat * cosmology->hubble_parameter(z);
	double vxrel = secondary.velocity.x - primary->velocity.x + hubble_flow;
	double vyrel = secondary.velocity.y - primary->velocity.y + hubble_flow;
	double vzrel = secondary.velocity.z - primary->velocity.z + hubble_flow;
	double vrel = std::sqrt(vxrel*vxrel + vyrel*vyrel + vzrel*vzrel);

	auto rvir_prim = darkmatterhalos->halo_virial_radius(primary->host_halo, z);
	auto rho_cen = primary->hot_halo_gas.mass / (shark::constants::PI4 * std::pow(rvir_prim,2) * rsat) / 1e18 ; //in Msun/pc^3
	// Here use the sum of the current hot halo gas plus what has been stripped. This assumed that the density of gas is only affected by cooling
	// and not ram pressure stripping.
	auto enc_mass = darkmatterhalos->enclosed_total_mass(secondary, z, r);
	double func = shark::constants::G * enc_mass *
			(secondary.hot_halo_gas.mass + secondary.hot_halo_gas_stripped.mass) / (8 * secondary.rvir_infall * std::pow(r,3)) / 1e18 -
			rho_cen * std::pow(vrel,2);

	return func;
}


} // namespace shark
