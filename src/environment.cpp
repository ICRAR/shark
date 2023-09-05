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
	double ram_pressure;
	double x_low;
	bool halo_strip;
	bool ism_strip;
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
	options.load("environment.alpha_rps_halo", alpha_rps_halo);

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

	// Ejected gas of the satellite subhalo is moved to the lost galaxy gas of the central. This is because Wright et al. (2020) found that 
	// the ejected gas of halos that become satellites make a negligible contribution to the gas accretion of halos.
	central_subhalo->lost_galaxy_gas += satellite_subhalo.lost_galaxy_gas + satellite_subhalo.ejected_galaxy_gas;
	central_subhalo->stellar_halo += satellite_subhalo.stellar_halo;
	central_subhalo->mean_galaxy_making_stellar_halo += satellite_subhalo.mean_galaxy_making_stellar_halo;

	satellite_subhalo.ejected_galaxy_gas.restore_baryon();
	satellite_subhalo.lost_galaxy_gas.restore_baryon();
	satellite_subhalo.stellar_halo.restore_baryon();
	satellite_subhalo.mean_galaxy_making_stellar_halo = 0;

	double ram_press = 0;
	double r_rps = 0;

	auto satellite_galaxy = satellite_subhalo.type1_galaxy();

	// Assume halo gas of subhalos that will disappear in the next snapshot is fully stripped.
	if(parameters.stripping && satellite_subhalo.infall_t == z){
		satellite_subhalo.transfer_halo_gas_to(central_subhalo);
	}

	if(parameters.stripping && satellite_subhalo.infall_t != z){
		// If I'm computing gradial ram pressure stripping of any form, then compute the ram pressure the satellite feels.
		if(parameters.gradual_stripping_halo || parameters.gradual_stripping_ism){
			ram_press = ram_pressure(central_subhalo, satellite_subhalo, z);
		}
	}

	if(parameters.stripping && satellite_subhalo.infall_t != z){
		// Treatment for all satellite subhalos.
		// Remove hot halo gas only if stripping is applied.

		auto mhot_tot = satellite_subhalo.cold_halo_gas.mass + satellite_subhalo.hot_halo_gas.mass + satellite_subhalo.hot_halo_gas_stripped.mass;
		if(mhot_tot > 0){
			if(parameters.gradual_stripping_halo){
				//first check whether the function is positive at Rvir_infall. In that case, the satellite subhalo experiences no stripping:
				//second, check whether function is negative at Rvir_infall/100. In that case assume all hot gas is stripped.
				auto func_rvir = ram_pressure_stripping_hot_gas(central_subhalo, satellite_subhalo, satellite_subhalo.rvir_infall, z, ram_press);
				auto func_rvirdiv100 = ram_pressure_stripping_hot_gas(central_subhalo, satellite_subhalo, satellite_subhalo.rvir_infall/100, z, ram_press);

				if(func_rvir > 0){
					r_rps = satellite_subhalo.hot_halo_gas_r_rps;
				}
				else if (func_rvir < 0 && func_rvirdiv100 > 0){
					r_rps = process_ram_pressure_stripping_gas(central_subhalo, satellite_subhalo, z, ram_press, true, false);
				}
				else if (func_rvir < 0 && func_rvirdiv100 < 0){
					r_rps = 0;
				}

				// If the ram-pressure stripping radius has decreased from previous timesteps, then compute how much new gas is lost.
				if(r_rps < satellite_subhalo.hot_halo_gas_r_rps && r_rps > 0){
					// 1. compute hot gas outside r_rps
					// 2. update satellite subhalo r_rps
					// 3. update hot gas that has been stripped.
					// 4. move that extra hot gas to central subhalo's reservoir.

					auto mhot_removed = mhot_tot * (1 - std::pow(r_rps/satellite_subhalo.rvir_infall,2)) - satellite_subhalo.hot_halo_gas_stripped.mass;
					if(mhot_removed < 0){
						mhot_removed = 0;
					}

					if(mhot_removed > satellite_subhalo.hot_halo_gas.mass + satellite_subhalo.cold_halo_gas.mass) {
						mhot_removed = satellite_subhalo.hot_halo_gas.mass + satellite_subhalo.cold_halo_gas.mass;
					}

					// remove halo gas in the proportion of the cold-to-hot halo mass ratio
					auto frac_gas_phase = satellite_subhalo.hot_halo_gas.mass / (satellite_subhalo.hot_halo_gas.mass + satellite_subhalo.cold_halo_gas.mass);
					auto metals_removed_hot = remove_gas(satellite_subhalo.hot_halo_gas, mhot_removed, frac_gas_phase);
					auto metals_removed_cold = remove_gas(satellite_subhalo.cold_halo_gas, mhot_removed, (1-frac_gas_phase));


					// track removed hot halo gas and ram-pressure stripping radius.
					satellite_subhalo.hot_halo_gas_r_rps = r_rps;
					satellite_subhalo.hot_halo_gas_stripped.mass += mhot_removed;
					satellite_subhalo.hot_halo_gas_stripped.mass_metals += metals_removed_cold + metals_removed_hot;

					central_subhalo->hot_halo_gas.mass += mhot_removed;
					central_subhalo->hot_halo_gas.mass_metals += metals_removed_cold + metals_removed_hot;

				}
				else if(r_rps == 0){
					// track stripping of gas
					satellite_subhalo.hot_halo_gas_stripped += (satellite_subhalo.hot_halo_gas + satellite_subhalo.cold_halo_gas);
					satellite_subhalo.hot_halo_gas_r_rps = satellite_subhalo.rvir_infall/100;

					//now transfer gas mass
					satellite_subhalo.transfer_halo_gas_to(central_subhalo);
				}
			}
			else{
				// instantaneous stripping assumed here. All halo gas is transferred.
				satellite_subhalo.transfer_halo_gas_to(central_subhalo);
			}
		}

		// Compute ram-pressure stripping of ISM gas if this option is turned on by the user.
		if(parameters.gradual_stripping_ism){
			auto mgal_gas = satellite_galaxy->disk_gas.mass + satellite_galaxy->bulge_gas.mass;
			auto mgal_gas_metals = satellite_galaxy->disk_gas.mass_metals + satellite_galaxy->bulge_gas.mass_metals;

			if(mgal_gas > 0){
				//first check whether the function is positive at Rvir_infall. In that case, the satellite subhalo experiences no stripping:
				//second, check whether function is negative at Rvir_infall/100. In that case assume all hot gas is stripped.
				auto func_rvir = ram_pressure_stripping_galaxy_gas(satellite_galaxy, satellite_subhalo.rvir_infall, z, ram_press);
				auto func_rvirdiv500 = ram_pressure_stripping_galaxy_gas(satellite_galaxy, satellite_subhalo.rvir_infall/500, z, ram_press);

				if(func_rvir > 0){
					r_rps = satellite_galaxy->r_rps;
				}
				else if (func_rvir < 0 && func_rvirdiv500 > 0){
					r_rps = process_ram_pressure_stripping_gas(central_subhalo, satellite_subhalo, z, ram_press, false, true);
				}
				else if (func_rvir < 0 && func_rvirdiv500 < 0){
					r_rps = 0;
				}

				// If the ram-pressure stripping radius has decreased from previous timesteps, then compute how much new gas is lost.
				if(r_rps < satellite_galaxy->r_rps && r_rps > 0){

					auto mism_removed_disk = satellite_galaxy->disk_gas.mass - satellite_galaxy->enclosed_mass_exponential(r_rps, satellite_galaxy->disk_gas.mass, satellite_galaxy->disk_gas.rscale);
					auto mism_removed_bulge = satellite_galaxy->bulge_gas.mass - satellite_galaxy->enclosed_mass_exponential(r_rps, satellite_galaxy->bulge_gas.mass, satellite_galaxy->bulge_gas.rscale);

					if(mism_removed_disk < 0){
						mism_removed_disk = 0;
					}
					if(mism_removed_bulge < 0){
						mism_removed_bulge = 0;
					}

					if(mism_removed_disk > satellite_galaxy->disk_gas.mass) {
						mism_removed_disk = satellite_galaxy->disk_gas.mass;
					}
					if(mism_removed_bulge > satellite_galaxy->bulge_gas.mass) {
						mism_removed_bulge = satellite_galaxy->bulge_gas.mass;
					}

					//stripped gas mass in proportion to the mass ratio of gas.
					auto metals_removed_bulge = remove_gas(satellite_galaxy->bulge_gas, mism_removed_bulge, 1);
					auto metals_removed_disk = remove_gas(satellite_galaxy->disk_gas, mism_removed_disk, 1);

					// Transfer mass to halo gas of central subhalo
					central_subhalo->hot_halo_gas.mass += mism_removed_disk + mism_removed_bulge;
					central_subhalo->hot_halo_gas.mass_metals += metals_removed_bulge + metals_removed_disk;

					// Keep track of mass loss and ram pressure stripping radii.
					satellite_galaxy->ram_pressure_stripped_gas.mass += mism_removed_disk + mism_removed_bulge;
					satellite_galaxy->ram_pressure_stripped_gas.mass_metals += metals_removed_bulge + metals_removed_disk;
					satellite_galaxy->r_rps = r_rps;
				}
				else if(r_rps == 0){
					// in this case all gas is transferred to the ISM gas to the central subhalo and restore baryon components.
					central_subhalo->hot_halo_gas.mass += mgal_gas;
					central_subhalo->hot_halo_gas.mass_metals += mgal_gas_metals;
					satellite_galaxy->disk_gas.restore_baryon();
					satellite_galaxy->bulge_gas.restore_baryon();

					// track mass loss
					satellite_galaxy->ram_pressure_stripped_gas.mass += mgal_gas;
					satellite_galaxy->ram_pressure_stripped_gas.mass_metals += mgal_gas_metals;
					satellite_galaxy->r_rps = satellite_subhalo.rvir_infall/500;
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

double Environment::process_ram_pressure_stripping_gas(const SubhaloPtr &primary,
		Subhalo &secondary,
		double z,
		double ram_press,
		bool halo_strip,
		bool ism_strip){

	double x_low = 0;

	if(halo_strip){
		x_low = secondary.rvir_infall/100.0;
	}
	else if(ism_strip){
		x_low = secondary.rvir_infall/500.0;
	}

	galaxy_properties_for_root_solver props = {
		z,
		ram_press,
		x_low,
		halo_strip,
		ism_strip,
		primary,
		secondary,
	};

	struct EnvironmentProcessAndProps {
		Environment *environment;
		galaxy_properties_for_root_solver *props;
	};

	auto f = [](double r, void *ctx) -> double {
		auto *env_and_props = static_cast<EnvironmentProcessAndProps *>(ctx);
		if(env_and_props->props->halo_strip){
			return env_and_props->environment->ram_pressure_stripping_hot_gas(env_and_props->props->primary,
						env_and_props->props->secondary,
						r,
						env_and_props->props->z,
						env_and_props->props->ram_pressure);
		}
		else if(env_and_props->props->ism_strip){
			return env_and_props->environment->ram_pressure_stripping_galaxy_gas(env_and_props->props->secondary.type1_galaxy(),
						r,
						env_and_props->props->z,
						env_and_props->props->ram_pressure);
		}
		else{
			std::ostringstream os;
			os << " ram pressure stripping calculation requires either halo_trip or ism_strip to be true";
			throw invalid_data(os.str());
		}
	};

	EnvironmentProcessAndProps env_and_props = {this, &props};
	double result = 0;
	try{
		result = root_solver.root_solver_function(f, &env_and_props, env_and_props.props->x_low, env_and_props.props->secondary.rvir_infall, 0, parameters.Accuracy_RPS);
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

	return result;

}

double Environment::ram_pressure_stripping_hot_gas(const SubhaloPtr &primary,
		const Subhalo &secondary,
		double r,
		double z,
		double ram_press){
	/* here we evaluate the function that needs to go to zero to find the best R_RPS solution for the hot gas
	this comes from the following equation:
	rho_hot_cen(Rsat) * vsat^2.0 = G * Msat(<r)*mhot_sat / (8 * Rvir_sat * r^3)
	Here, Rsat is the distance from the halo centre to the satellite galaxy; vsat is the relative velocity of the subhalo and r what we need to solve for.
	*/

	// Here use the sum of the current hot halo gas plus what has been stripped. This assumed that the density of gas is only affected by cooling
	// and not ram pressure stripping.
	auto enc_mass = darkmatterhalos->enclosed_total_mass(secondary, z, r);
	double func = parameters.alpha_rps_halo * shark::constants::G * enc_mass *
			(secondary.hot_halo_gas.mass + secondary.hot_halo_gas_stripped.mass) / (8 * secondary.rvir_infall * std::pow(r,3)) / 1e18 -
			ram_press;

	return func;
}

double Environment::ram_pressure_stripping_galaxy_gas(const GalaxyPtr &galaxy,
		double r,
		double z,
		double ram_press){
	/* here we evaluate the function that needs to go to zero to find the best R_RPS solution for the hot gas
	this comes from the following equation:
	rho_hot_cen(Rsat) * vsat^2.0 = 2 * PI * G * Sigma_gas(r) * (Sigma_gas(r) + Sigma_star(r))
	Here, Rsat is the distance from the halo centre to the satellite galaxy; vsat is the relative velocity of the subhalo and r what we need to solve for.
	*/

	// Here use the sum of the current hot halo gas plus what has been stripped. This assumed that the density of gas is only affected by cooling
	// and not ram pressure stripping.

	auto sigma_gas = galaxy->surface_density_gas(r) / 1e12; //In Msun/pc^2
	auto sigma_gal = (galaxy->surface_density_bulge(r) + galaxy->surface_density_disk(r)) / 1e12; //In Msun/pc^2
	double func = shark::constants::PI2 *  shark::constants::G * sigma_gas * sigma_gal -
			ram_press;

	return func;
}

double Environment::ram_pressure(const SubhaloPtr &primary,
		const Subhalo &secondary,
		double z)
{
	// Compute ram pressure in units of Msun/pc^3 * (km/s)^2

	// Find the objects physical position
	double conversion_factor = cosmo_params.Hubble_h * (1 +  z);
	auto pos_rel = (secondary.position - primary->position) / conversion_factor;
	auto rsat = pos_rel.norm();

	// Find the physical velocity + the hubble flow
	double hubble_flow = rsat * cosmology->hubble_parameter(z);
	auto vrel = secondary.velocity - primary->velocity + hubble_flow;
	auto vrel_norm = vrel.norm();

	auto rvir_prim = darkmatterhalos->halo_virial_radius(primary->host_halo, z);
	auto rho_cen = primary->hot_halo_gas.mass / (shark::constants::PI4 * std::pow(rvir_prim,2) * rsat) / 1e18 ; //in Msun/pc^3

	return rho_cen * std::pow(vrel_norm,2);
}

float Environment::remove_gas(BaryonBase &component, double m_removed, float f_gas){

	if(component.mass > 0){
		auto z = component.mass_metals / component.mass;
		auto metals_removed = z * m_removed * f_gas;
		component.mass_metals -= metals_removed;
		component.mass -= f_gas * m_removed;

		if(component.mass < 0){
			component.restore_baryon();
		}
		if(component.mass_metals < 0){
			component.mass_metals = 0;
		}

		return metals_removed;
	}
	else{
		return 0;
	}

}


} // namespace shark
