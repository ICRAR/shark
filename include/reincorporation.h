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

/**
 * @file
 */

#ifndef INCLUDE_REINCORPORATION_H_
#define INCLUDE_REINCORPORATION_H_

#include <memory>
#include <utility>

#include "components.h"
#include "dark_matter_halos.h"
#include "options.h"

namespace shark {

class ReincorporationParameters{

public:
	ReincorporationParameters(const Options &options);

	double tau_reinc = 0;
	double mhalo_norm = 0;
	double halo_mass_power = 0;

};

class Reincorporation{

public:
	Reincorporation(const ReincorporationParameters &parameters, const DarkMatterHalosPtr &darkmatterhalo);

	double reincorporated_mass (Halo &halo, Subhalo &subhalo, double z, double delta_t);

private:

	ReincorporationParameters parameters;
	DarkMatterHalosPtr darkmatterhalo;
};

using ReincorporationPtr = std::shared_ptr<Reincorporation>;

template <typename ...Ts>
ReincorporationPtr make_reincorporation(Ts&&...ts)
{
	return std::make_shared<Reincorporation>(std::forward<Ts>(ts)...);
}

}//end namespace shark

#endif /* INCLUDE_REINCORPORATION_H_ */
