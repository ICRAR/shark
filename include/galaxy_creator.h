//
// Galaxy creator class declaration
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
#ifndef SHARK_GALAXY_CREATOR_H_
#define SHARK_GALAXY_CREATOR_H_

#include "cosmology.h"
#include "components.h"
#include "dark_matter_halos.h"
#include "gas_cooling.h"
#include "simulation.h"

namespace shark {

class GalaxyCreator {

public:
	GalaxyCreator(const CosmologyPtr &cosmology, GasCoolingParameters cool_params, SimulationParameters sim_params);
	void create_galaxies(const std::vector<MergerTreePtr> &merger_trees, TotalBaryon &AllBaryons);

private:
	bool create_galaxies(const HaloPtr &halo, double z, id_t ID);

	CosmologyPtr cosmology;
	GasCoolingParameters cool_params;
	SimulationParameters sim_params;
};

}  // namespace shark

#endif // SHARK_GALAXY_CREATOR_H_
