/*
 * SHArk_Evolve_Halos.cpp
 *
 *  Created on: 10Apr.,2017
 *      Author: clagos
 */

#include <cmath>
#include <memory>

#include "components.h"
#include "evolve_halos.h"
#include "logging.h"

using namespace std;

namespace shark {

static
void evolve_system( shared_ptr<BasicPhysicalModel> physicalmodel, SubhaloPtr &subhalo, int snapshot, double z, double delta_t){

	// Solve ODEs for this system
	for(auto &galaxy: subhalo->galaxies) {
		physicalmodel->evolve_galaxy(*subhalo, *galaxy, z, delta_t);
		//Solve_Systems();
	}

}

void populate_halos(shared_ptr<BasicPhysicalModel> physicalmodel, HaloPtr halo, int snapshot, double z, double delta_t) {


	for(auto &subhalo: halo->all_subhalos()) {
		evolve_system(physicalmodel, subhalo, snapshot, z, delta_t);
	}
}

void transfer_galaxies_to_next_snapshot(const std::vector<HaloPtr> &halos){

	/**
	 * This function transfer galaxies of the subhalos of this snapshot into the subhalos of the next snapshot, and baryon components from subhalo to subhalo.
	 */
	unsigned int subhalos_without_descendant = 0;
	for(auto &halo: halos){
		for(SubhaloPtr &subhalo: halo->all_subhalos()) {

			auto descendant_subhalo = subhalo->descendant;
			if (!descendant_subhalo) {
				subhalos_without_descendant++;
				continue;
			}

			// Check cases where the descendant subhalo will be a satellite, but the current is central. In that case
			// we modify the type of the central galaxy of this subhalo to type1. and the rest to type 2.

			if(subhalo->subhalo_type == Subhalo::CENTRAL && descendant_subhalo->subhalo_type == Subhalo::SATELLITE){
				for (auto &galaxy: subhalo->galaxies){
					if(galaxy->galaxy_type == Galaxy::CENTRAL){
						galaxy->galaxy_type = Galaxy::TYPE1;
					}
					else{
						galaxy->galaxy_type = Galaxy::TYPE2;
					}
				}
			}

			// Transfer galaxies.
			subhalo->transfer_galaxies_to(descendant_subhalo);

			// Transfer subhalo baryon components.
			descendant_subhalo->cold_halo_gas = subhalo->cold_halo_gas;
			descendant_subhalo->hot_halo_gas = subhalo->hot_halo_gas;
			descendant_subhalo->ejected_galaxy_gas = subhalo->ejected_galaxy_gas;
			descendant_subhalo->cooling_subhalo_tracking = subhalo->cooling_subhalo_tracking;

		}
	}

	if (subhalos_without_descendant) {
		LOG(warning) << "Found " << subhalos_without_descendant << " subhalos without descendant while transferring galaxies";
	}

}

}
