/*
 * simulation.cpp
 *
 *  Created on: 13Jun.,2017
 *      Author: clagos
 */

#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <tuple>

#include "components.h"
#include "exceptions.h"
#include "logging.h"
#include "numerical_constants.h"
#include "simulation.h"


namespace shark {

SimulationParameters::SimulationParameters(const Options &options)
{

	std::string redshift_file;

	options.load("simulation.volume", volume, true);
	options.load("simulation.particle_mass", particle_mass);
	options.load("simulation.lbox", lbox, true);
	options.load("simulation.tot_n_subvolumes", tot_nsubvols, true);
	options.load("simulation.min_snapshot", min_snapshot, true);
	options.load("simulation.max_snapshot", max_snapshot, true);
	options.load("simulation.sim_name", sim_name);
	options.load("simulation.tree_files_prefix", tree_files_prefix, true);
	options.load("simulation.redshift_file",redshift_file, true);

	load_simulation_tables(redshift_file);

}

void SimulationParameters::load_simulation_tables(const std::string &redshift_file)
{

	using namespace std;

	LOG(debug) << "Reading table " << redshift_file ;

	ifstream f = open_file(redshift_file);
	string line;
	while ( getline(f, line) ) {

		trim(line);
		if (empty_or_comment(line)) {
			continue;
		}

		int s;
		double r;

		istringstream iss(line);
		iss >> s >> r;

		redshifts[s]=r;

	}
	f.close();

	// Check that the redshift values descend when snapshots ascend
	// This also implies that there are no repeated values, which we have seen
	// which we have seen in the past
	auto prev_z = std::numeric_limits<double>::max();
	for(auto &pair: redshifts) {
		auto z = pair.second;
		if (z >= prev_z) {
			std::ostringstream os;
			os << "redshift at snapshot " << pair.first << " (" << z << ") is not lower than ";
			os << "redshift at snapshot " << pair.first - 1 << " (" << prev_z << ")";
			throw invalid_data(os.str());
		}
		prev_z = z;
	}

}

Simulation::Simulation(SimulationParameters parameters, const CosmologyPtr &cosmology) :
	parameters(parameters),
	cosmology(cosmology)
{
	// no-op
}

double Simulation::convert_snapshot_to_age(int s){

	return cosmology->convert_redshift_to_age(parameters.redshifts[s]);

}
}



