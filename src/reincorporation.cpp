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

#include "halo.h"
#include "reincorporation.h"
#include "subhalo.h"

namespace shark {

ReincorporationParameters::ReincorporationParameters(const Options &options)
{
	options.load("reincorporation.tau_reinc",tau_reinc);
	options.load("reincorporation.mhalo_norm",mhalo_norm);
	options.load("reincorporation.halo_mass_power",halo_mass_power);
}

Reincorporation::Reincorporation(const ReincorporationParameters &parameters, DarkMatterHalosPtr darkmatterhalo):
	parameters(parameters),
	darkmatterhalo(std::move(darkmatterhalo))
{
	//no-opt
}

double Reincorporation::reincorporated_mass(Halo &halo, Subhalo &subhalo, double z, double delta_t){

	if(subhalo.subhalo_type == Subhalo::SATELLITE){
		return 0;
	}

	double mvir = halo.Mvir;
	double mgas = subhalo.ejected_galaxy_gas.mass;

	double treinc = parameters.tau_reinc * std::pow( (mvir / parameters.mhalo_norm), parameters.halo_mass_power);
	//double reincor_rate = mgas * parameters.tau_reinc/tdyn * std::pow( (mvir / parameters.mhalo_norm), parameters.halo_mass_power);

	if(treinc == 0){
		return mgas;
	}
	else if (treinc > 100){
		return 0;
	}

	return mgas / treinc * delta_t;

}

} // namespace shark
