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
 * Galaxy creator class implementation
 */

#include <atomic>
#include <numeric>

#include "galaxy_creator.h"
#include "logging.h"
#include "omp_utils.h"
#include "timer.h"

namespace shark {

GalaxyCreator::GalaxyCreator(const CosmologyPtr &cosmology, GasCoolingParameters cool_params, SimulationParameters sim_params, unsigned int threads) :
	cosmology(cosmology),
	cool_params(cool_params),
	sim_params(sim_params),
	threads(threads)
{
	// no-op
}

void GalaxyCreator::create_galaxies(const std::vector<MergerTreePtr> &merger_trees, TotalBaryon &AllBaryons)
{
	std::atomic<Galaxy::id_t> galaxy_id(0);
	auto timer = Timer();

	std::vector<double> total_baryon(threads);

	for(int snapshot = sim_params.min_snapshot; snapshot <= sim_params.max_snapshot - 1; snapshot++) {
		auto z = sim_params.redshifts[snapshot];

		omp_static_for(merger_trees, threads, [&](const MergerTreePtr &merger_tree, int thread_num) {
			for(auto &halo: merger_tree->halos_at(snapshot)) {

				// Halo has a central subhalo with ascendants so ignore it, as it should already have galaxies in it.
				if(halo->central_subhalo->ascendants.size() > 0){
					continue;
				}

				create_galaxies(halo, z, galaxy_id++);
				total_baryon[thread_num] += halo->central_subhalo->hot_halo_gas.mass;
			}
		});

		// Keep track of the total amount of baryons integrated from the first snapshot to the current one.
		AllBaryons.baryon_total_created[snapshot] += std::accumulate(total_baryon.begin(), total_baryon.end(), 0.);
	}

	LOG(info) << "Created " << galaxy_id.load() << " initial galaxies in " << timer;
}

void GalaxyCreator::create_galaxies(const HaloPtr &halo, double z, Galaxy::id_t galaxy_id)
{
	// Central subhalo has a central galaxy (somehow!), ignore
	auto central_subhalo = halo->central_subhalo;
	if(central_subhalo->central_galaxy()) {
		std::ostringstream os;
		os << "Central Subhalo " << central_subhalo << " is first in merger tree but has central galaxy.";
		throw invalid_argument(os.str());
		//return false;
	}

	// Count how many galaxies this halo has.
	auto galaxy_count = central_subhalo->galaxy_count();
	if(galaxy_count > 0){
		std::ostringstream os;
		os << "Central Subhalo " << central_subhalo << " has no central galaxy but " << galaxy_count <<" satellites.";
		throw invalid_argument(os.str());
	}

	auto galaxy = std::make_shared<Galaxy>(galaxy_id);
	galaxy->vmax = central_subhalo->Vcirc;

	central_subhalo->hot_halo_gas.mass = halo->Mvir * cosmology->universal_baryon_fraction();

	// Assign metallicity to the minimum allowed.
	central_subhalo->hot_halo_gas.mass_metals = central_subhalo->hot_halo_gas.mass * cool_params.pre_enrich_z;

	central_subhalo->galaxies.emplace_back(std::move(galaxy));

	if (LOG_ENABLED(debug)) {
		LOG(debug) << "Added a central galaxy for subhalo " << central_subhalo;
	}

}

}
