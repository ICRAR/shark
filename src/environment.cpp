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

namespace shark {

EnvironmentParameters::EnvironmentParameters(const Options &options)
{
	options.load("environment.gradual_stripping", gradual_stripping);
	options.load("environment.stripping", stripping, true);
}

Environment::Environment(const EnvironmentParameters &parameters):
	parameters(parameters){
	//no-opt
}

void Environment::process_satellite_subhalo_environment(Subhalo &satellite_subhalo, Subhalo &central_subhalo){

	if(parameters.stripping){
		// Ejected gas is moved to the budget of ejected gas of the central, as this gas escaped
		// the subhalo of the satellite.
		central_subhalo.ejected_galaxy_gas += satellite_subhalo.ejected_galaxy_gas;
		central_subhalo.lost_galaxy_gas += satellite_subhalo.lost_galaxy_gas;

		satellite_subhalo.ejected_galaxy_gas.restore_baryon();
		satellite_subhalo.lost_galaxy_gas.restore_baryon();

		if(parameters.gradual_stripping){
			//TODO
		}
		else{
			if(satellite_subhalo.hot_halo_gas.mass > 0 || satellite_subhalo.ejected_galaxy_gas.mass > 0 || satellite_subhalo.cold_halo_gas.mass > 0){

				central_subhalo.hot_halo_gas += satellite_subhalo.hot_halo_gas;
				central_subhalo.hot_halo_gas += satellite_subhalo.cold_halo_gas;

				satellite_subhalo.hot_halo_gas.restore_baryon();
				satellite_subhalo.cold_halo_gas.restore_baryon();
			}
		}
	}

}

} // namespace shark
