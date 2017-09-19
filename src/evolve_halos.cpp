/*
 * SHArk_Evolve_Halos.cpp
 *
 *  Created on: 10Apr.,2017
 *      Author: clagos
 */

#include <cmath>
#include <memory>

#include "evolve_halos.h"

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

void transfer_galaxies_to_next_snapshot(HaloPtr halo){

	/**
	 * This function transfer galaxies of the subhalos of this snapshot into the subhalos of the next snapshot, and baryon components from subhalo to subhalo.
	 */
	for(SubhaloPtr &subhalo: halo->all_subhalos()) {

		auto descendant_subhalo = subhalo->descendant;

		// Transfer galaxies.
		subhalo->copy_galaxies_to(descendant_subhalo);

		// Transfer subhalo baryon components.
		descendant_subhalo->cold_halo_gas = subhalo->cold_halo_gas;
		descendant_subhalo->hot_halo_gas = subhalo->hot_halo_gas;
		descendant_subhalo->ejected_galaxy_gas = subhalo->ejected_galaxy_gas;
		descendant_subhalo->cooling_subhalo_tracking = subhalo->cooling_subhalo_tracking;

	}


}

void destroy_galaxies_this_snapshot(vector<HaloPtr> all_halos_this_snapshot){

	for(auto &halo: all_halos_this_snapshot){

		for(auto &subhalo: halo->all_subhalos()) {
			subhalo->galaxies.clear();
		}

	}

}


}
