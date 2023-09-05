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
 * Halo evolution routines. Most of it now lives under the SharkRunner class,
 * but some things have not been moved yet.
 */

#ifndef SHARK_EVOLVE_HALOS_H_
#define SHARK_EVOLVE_HALOS_H_

#include <cmath>
#include <memory>

#include "components.h"
#include "cosmology.h"
#include "execution.h"
#include "physical_model.h"
#include "simulation.h"
#include "star_formation.h"

namespace shark {

/**
 * Transfers galaxies of the subhalos of this snapshot into the corresponding
 * subhalos of the next snapshot, and baryon components from subhalo to subhalo.
 *
 * @param halos The halos whose subhalos need to be transferred to the next snapshot
 * @param snapshot This snapshot
 * @param AllBaryons The TotalBaryon accummulation object
 */
void transfer_galaxies_to_next_snapshot(const std::vector<HaloPtr> &halos, int snapshot, TotalBaryon &AllBaryons);

void track_total_baryons(Cosmology &cosmology, ExecutionParameters execparams, SimulationParameters simulation_params, const std::vector<HaloPtr> &halos,
		TotalBaryon &AllBaryons, int snapshot, const molgas_per_galaxy &molgas, double deltat);

void reset_instantaneous_galaxy_properties(const std::vector<HaloPtr> &halos, int snapshot);

}  // namespace shark

#endif // SHARK_EVOLVE_HALOS_H_
