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

#ifndef SHARK_GALAXY_CREATOR_H_
#define SHARK_GALAXY_CREATOR_H_

#include "cosmology.h"
#include "components.h"
#include "dark_matter_halos.h"
#include "galaxy.h"
#include "gas_cooling.h"
#include "simulation.h"

namespace shark {

class GalaxyCreator {

public:
	GalaxyCreator(CosmologyPtr cosmology, GasCoolingParameters cool_params, SimulationParameters sim_params);
	void create_galaxies(const std::vector<MergerTreePtr> &merger_trees, TotalBaryon &AllBaryons);

private:
	bool create_galaxies(const HaloPtr &halo, double z, Galaxy::id_t ID);

	CosmologyPtr cosmology;
	GasCoolingParameters cool_params;
	SimulationParameters sim_params;
};

}  // namespace shark

#endif // SHARK_GALAXY_CREATOR_H_
