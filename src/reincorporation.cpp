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
 */

#include <cmath>

#include "reincorporation.h"

namespace shark {

ReincorporationParameters::ReincorporationParameters(const Options &options)
{
	options.load("reincorporation.alpha_reheat",alpha_reheat);
	options.load("reincorporation.mhalo_norm",mhalo_norm);
	options.load("reincorporation.halo_mass_power",halo_mass_power);
}

Reincorporation::Reincorporation(const ReincorporationParameters &parameters, const DarkMatterHalosPtr &darkmatterhalo):
	parameters(parameters),
	darkmatterhalo(darkmatterhalo)
{
	//no-opt
}

double Reincorporation::reincorporated_mass(HaloPtr halo, double z, double delta_t){

	double mvir = halo->Mvir;
	double mgas = halo->central_subhalo->ejected_galaxy_gas.mass;

	double treinc = parameters.alpha_reheat * std::pow( (mvir / parameters.mhalo_norm), parameters.halo_mass_power);
	//double reincor_rate = mgas * parameters.alpha_reheat/tdyn * std::pow( (mvir / parameters.mhalo_norm), parameters.halo_mass_power);

	if(treinc == 0){
		return mgas;
	}

	return mgas / treinc * delta_t;

}

}
