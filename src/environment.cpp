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

#include "environment.h"
#include "subhalo.h"
#include "galaxy.h"

namespace shark {

EnvironmentParameters::EnvironmentParameters(const Options &options)
{
	options.load("environment.stripping", stripping, true);
	options.load("environment.gradual_stripping", gradual_stripping);
	options.load("environment.tidal_stripping", tidal_stripping);
	options.load("environment.minimum_halo_mass_fraction", minimum_halo_mass_fraction);

}

Environment::Environment(const EnvironmentParameters &parameters):
	parameters(parameters){
	//no-opt
}

void Environment::process_satellite_subhalo_environment(Subhalo &satellite_subhalo, Subhalo &central_subhalo){

	// Ejected gas is moved to the budget of ejected gas of the central, as this gas escaped
	// the subhalo of the satellite. This is always the case.
	central_subhalo.ejected_galaxy_gas += satellite_subhalo.ejected_galaxy_gas;
	central_subhalo.lost_galaxy_gas += satellite_subhalo.lost_galaxy_gas;

	satellite_subhalo.ejected_galaxy_gas.restore_baryon();
	satellite_subhalo.lost_galaxy_gas.restore_baryon();


	// Remove hot halo gas only if stripping is applied
	if(parameters.stripping){

		if(parameters.gradual_stripping){
			//TODO
		}
		else{
			if(satellite_subhalo.hot_halo_gas.mass > 0 || satellite_subhalo.cold_halo_gas.mass > 0){

				central_subhalo.hot_halo_gas += satellite_subhalo.hot_halo_gas;
				central_subhalo.hot_halo_gas += satellite_subhalo.cold_halo_gas;

				satellite_subhalo.hot_halo_gas.restore_baryon();
				satellite_subhalo.cold_halo_gas.restore_baryon();
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
		float lost_stellar_mass = 0;
		float lost_stellar_mass_metals = 0;

		if(satellite_subhalo.central_galaxy()){
			// Apply here model of Errani et al. (2015).
			float ratio_mass = satellite_subhalo.Mvir / satellite_subhalo.Mvir_infall;
			// Apply a maximum stripping of 99% in halo mass.
			if(ratio_mass < parameters.minimum_halo_mass_fraction){
				ratio_mass = parameters.minimum_halo_mass_fraction;
			}

			// compute how much has been lost since galaxy infall
			float ratio_sm = std::pow(2, 3.43) * std::pow(ratio_mass, 1.86) / std::pow( 1 + ratio_mass, 3.43);
			lost_stellar_mass = (1 - ratio_sm) * satellite_subhalo.star_central_infall.mass;
			// now remove what has already been lost in previous snapshots to compute what should be subtracted from
			// the central galaxy in the satellite subhalo now.
			lost_stellar_mass = lost_stellar_mass - satellite_subhalo.central_galaxy()->stars_tidal_stripped.mass;

			if(lost_stellar_mass > 0){
				remove_tidal_stripped_stars(*satellite_subhalo.central_galaxy(), lost_stellar_mass, lost_stellar_mass_metals);

				// add the stripped material to central subhalo
				central_subhalo.stellar_halo.mass += lost_stellar_mass;
				central_subhalo.stellar_halo.mass_metals += lost_stellar_mass_metals;
			}
		}
		// check all type 2 galaxies in the central subhalo and strip only those that haven't been stripped yet.
		if(central_subhalo.type2_galaxies_count() > 0){
			float ratio_sm = std::pow(2.0, 3.43) * std::pow(parameters.minimum_halo_mass_fraction, 1.86) / std::pow( 1 + parameters.minimum_halo_mass_fraction, 3.43);
			// in this case loop over type II satellites and if they have not been stripped
			// then strip a fraction of the satellite mass.
			for (auto &satellite: central_subhalo.type2_galaxies() ){
				// If no stripping has occurred in this type II then assume 99% loss of DM.
				// is stars_tidal_stripped.mass>0 is because stripping should have occurred already.
				if(satellite.stars_tidal_stripped.mass == 0){
					// compute how much has been lost since galaxy infall
					lost_stellar_mass = (1 - ratio_sm) * satellite.stellar_mass();

					remove_tidal_stripped_stars(satellite, lost_stellar_mass, lost_stellar_mass_metals);

					// add the stripped material to central subhalo
					central_subhalo.stellar_halo.mass += lost_stellar_mass;
					central_subhalo.stellar_halo.mass_metals += lost_stellar_mass_metals;
				}
			}
		}
	}

}

void Environment::remove_tidal_stripped_stars(Galaxy &galaxy, float lost_stellar_mass, float lost_stellar_mass_metals){


	//check that lost mass does not exceed the total stellar mass of the galaxy to be stripped.
	if(lost_stellar_mass > galaxy.stellar_mass()){
		lost_stellar_mass = galaxy.stellar_mass();
		lost_stellar_mass_metals = galaxy.stellar_mass_metals();
		galaxy.disk_stars.restore_baryon();
		galaxy.bulge_stars.restore_baryon();
	}
	else{

		// strip first the disk of the galaxy and then the bulge:
		if(lost_stellar_mass < galaxy.disk_stars.mass){
			// in this case we strip material from the disk but not the bulge
			lost_stellar_mass_metals = lost_stellar_mass / galaxy.disk_stars.mass * galaxy.disk_stars.mass_metals;
			galaxy.disk_stars.mass -= lost_stellar_mass;
			galaxy.disk_stars.mass_metals -= lost_stellar_mass_metals;
		}
		else{
			// in this case we strip all the disk and remove a fraction of the bulge.
			float lost_bulge_mass = lost_stellar_mass - galaxy.disk_stars.mass;
			float lost_bulge_metals = lost_bulge_mass / galaxy.bulge_stars.mass * galaxy.bulge_stars.mass_metals;
			lost_stellar_mass_metals = galaxy.disk_stars.mass_metals + lost_bulge_metals;
			galaxy.disk_stars.restore_baryon();

			galaxy.bulge_stars.mass -= lost_bulge_mass;
			galaxy.bulge_stars.mass_metals -= lost_bulge_metals;
		}

		// sanity checks
		if(galaxy.disk_stars.mass < 0){
			galaxy.disk_stars.restore_baryon();
		}
		if(galaxy.bulge_stars.mass < 0){
			galaxy.bulge_stars.restore_baryon();
		}
	}

	// accummulate what has been lost to tidal stripping in this galaxy.
	galaxy.stars_tidal_stripped.mass += lost_stellar_mass;
	galaxy.stars_tidal_stripped.mass_metals += lost_stellar_mass_metals;

}

} // namespace shark
