/*
 * SHArk_Evolve_Halos.cpp
 *
 *  Created on: 10Apr.,2017
 *      Author: clagos
 */

#include <cmath>
#include <memory>

#include "components.h"
#include "cosmology.h"
#include "numerical_constants.h"
#include "simulation.h"

using namespace std;

namespace shark {

static
void evolve_system(shared_ptr<Subhalo> &subhalo, int snapshot){

	//
	// Determine if galaxies in this subhalo merger on the current snapshot.
	//
	// Calculate_Merger_Galaxies();

	// Solve ODEs for this system
	for(shared_ptr<Galaxy> &galaxy: subhalo->galaxies) {
		//Solve_Systems();
	}

}

void populate_halos(shared_ptr<Halo> halo, int snapshot) {

	//
	// Determine if there is any exchange of galaxies between subhalos.
	// This could happen if galaxies go from type 0 to 2, and 1 to merging with the central in one snapshot.*/
	//
	//Calculates_Merger_Subhalos();

	for(shared_ptr<Subhalo> &subhalo: halo->subhalos) {
		evolve_system(subhalo, snapshot);
	}
}

}