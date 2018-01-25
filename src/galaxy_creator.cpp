//
// Galaxy creator class implementation
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

#include <chrono>

#include "galaxy_creator.h"
#include "logging.h"

namespace shark {

GalaxyCreator::GalaxyCreator(std::shared_ptr<Cosmology> cosmology, std::shared_ptr<DarkMatterHalos> darkmatterhalos, GasCoolingParameters cool_params, SimulationParameters sim_params) :
	cosmology(cosmology),
	darkmatterhalos(darkmatterhalos),
	cool_params(cool_params),
	sim_params(sim_params)
{
	// no-op
}

void GalaxyCreator::create_galaxies(const std::vector<MergerTreePtr> &merger_trees, TotalBaryon &AllBaryons)
{
	int galaxies_added = 0;
	double total_baryon = 0.0;

	auto start = std::chrono::steady_clock::now();
	for(int snapshot = sim_params.min_snapshot; snapshot <= sim_params.max_snapshot - 1; snapshot++) {
		for(auto &merger_tree: merger_trees) {
			for(auto &halo: merger_tree->halos[snapshot]) {
				if (create_galaxies(halo, sim_params.redshifts[snapshot])) {
					galaxies_added++;
					total_baryon += halo->central_subhalo->hot_halo_gas.mass;
				}
			}
		}
		// Keep track of the total amount of baryons integrated from the first snapshot to the current one.
		AllBaryons.baryon_total_created[snapshot] += total_baryon;
	}

	auto d = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
	LOG(info) << "Created " << galaxies_added << " initial galaxies in " << d << " [ms]";
}

bool GalaxyCreator::create_galaxies(const HaloPtr &halo, double z)
{

	// Halo has a central subhalo with ascendants so ignore it, as it should already have galaxies in it.
	if(halo->central_subhalo->ascendants.size() > 0){
		return false;
	}

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

	auto galaxy = std::make_shared<Galaxy>();
	galaxy->galaxy_type = Galaxy::CENTRAL;

	central_subhalo->galaxies.push_back(galaxy);
	LOG(debug) << "Added a central galaxy for subhalo " << central_subhalo;

	central_subhalo->hot_halo_gas.mass = halo->Mvir * cosmology->universal_baryon_fraction();

	// Assign metallicity to the minimum allowed.
	central_subhalo->hot_halo_gas.mass_metals = central_subhalo->hot_halo_gas.mass * cool_params.pre_enrich_z;

	//assign an ad-hoc half-mass radius and specific angular momentum to start with.
	galaxy->disk_gas.rscale = darkmatterhalos->disk_size_theory(*central_subhalo, z);

	darkmatterhalos->galaxy_velocity(*central_subhalo);

	return true;

}

}
