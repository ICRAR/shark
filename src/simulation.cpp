/*
 * simulation.cpp
 *
 *  Created on: 13Jun.,2017
 *      Author: clagos
 */

#include <cmath>
#include <fstream>
#include <map>
#include <tuple>

#include "simulation.h"
#include "logging.h"
#include "numerical_constants.h"
#include "components.h"


namespace shark {

SimulationParameters::SimulationParameters(const std::string &filename) :
	Options(filename),
	volume(0),
	particle_mass(0),
	min_snapshot(0),
	max_snapshot(0),
	redshifts()
	{

	std::string redshift_file;

	load("simulation.volume", volume);
	load("simulation.particle_mass", particle_mass);
	load("simulation.min_snapshot", min_snapshot);
	load("simulation.max_snapshot", max_snapshot);
	load("simulation.redshift_file",redshift_file);

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
		if (is_skipable(line)) {
			continue;
		}

		int s;
		double r;

		istringstream iss(line);
		iss >> s >> r;

		redshifts[s]=r;

	}
	f.close();

}

Simulation::Simulation(SimulationParameters parameters, std::shared_ptr<Cosmology> cosmology) :
	parameters(parameters),
	cosmology(cosmology)
{
	// no-op
}

double Simulation::convert_snapshot_to_age(int s){

	return cosmology->convert_redshift_to_age(parameters.redshifts[s]);

}
}



