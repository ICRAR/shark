//
// The main function for the shark executable
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

#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

#include "components.h"
#include "cosmology.h"
#include "evolve_halos.h"
#include "numerical_constants.h"
#include "parameters.h"
#include "physical_model.h"
#include "simulation.h"


using namespace shark;
using namespace std;

// TODO: All these are stub, dummy methods that will go away one by one
//       into their specific parts of the code
Parameters read_parameters(const string &name) {
	return Parameters();
}

SimulationParameters read_simulation_parameters(const string &name) {
	return SimulationParameters();
}

vector<MergerTree> load_merger_trees() {
	return vector<MergerTree>();
}

void physical_processes() {
	return;
}

void do_stuff_at_halo_level(vector<shared_ptr<Halo>> halos) {
	return;
}

void write_output(int snapshot, const vector<MergerTree> &merger_trees) {
	return;
}

/**
 * Main SHArk routine.
 *
 * Here we load the relevant information and do the basic loops to solve galaxy formation
 */
int main(int argc, char **argv) {

	if ( argc < 2 ) {
		cerr << "Usage: " << argv[0] << " <params-file>" << endl;
		return 1;
	}

	/* We read the parameters that have been given as input by the user.*/
	Parameters params = read_parameters(argv[1]);

	//GasCoolingParameters gas_cooling_params(argv[1]);
	GasCooling gas_cooling = GasCooling(GasCoolingParameters(argv[1]));
	StellarFeedback stellar_feedback = StellarFeedback(StellarFeedbackParameters(argv[1]));
	StarFormation star_formation = StarFormation(StarFormationParameters(argv[1]));
	RecyclingParameters recycling_parameters;
	BasicPhysicalModel basic_physicalmodel(1e-6, gas_cooling, stellar_feedback, star_formation, recycling_parameters);

	// We read the simulation parameters next. Note that by using a different
	// reader allows the user to put all the information in one parameter file.
	SimulationParameters sim_params = read_simulation_parameters(argv[1]);

	// Read the merger tree files.
	// Each merger tree will be a construction of halos and subhalos
	// with their growth history.
	vector<MergerTree> merger_trees = load_merger_trees();

	// This function should return the system of differential equations
	// to be solved at each snapshot and in each galaxy.
	// This set of ODEs apply on the ideal case of a central galaxy, with no
	// AGN feedback. These equations are modified later if galaxies are
	// satellites or have an AGN.
	physical_processes();

	// The way we solve for galaxy formation is snapshot by snapshot.
	// We first loop over those, and for a fixed snapshot,
	// we loop over merger trees.
	// Each merger trees has a set of halos at a given snapshot,
	// which in turn contain galaxies.
	for(int snapshot=sim_params.min_snapshots; snapshot <= sim_params.max_snapshots; snapshot++) {

		for(MergerTree &tree: merger_trees) {
			/*here loop over the halos this merger tree has at this time.*/
			for(shared_ptr<Halo> halo: tree.halos[snapshot]) {
				/*populate halos. This function should evolve the subhalos inside the halo.*/
				populate_halos(basic_physicalmodel, halo, snapshot);
			}
		}

		/*Here you shoud include the physics that allow halos to speak to each other. This could be useful e.g. during reionisation.*/
		vector<shared_ptr<Halo>> all_halos_for_this_snapshot;
		do_stuff_at_halo_level(all_halos_for_this_snapshot);

		/*write snapshots only if the user wants outputs at this time.*/
		if( std::find(params.writing_outputs.begin(), params.writing_outputs.end(), snapshot) != params.writing_outputs.end() ){
			write_output(snapshot, merger_trees);
		}
	}

	return 0;
}
