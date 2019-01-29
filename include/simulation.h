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
 * Simulation parameters used as input for shark
 */

#ifndef SHARK_SIMULATION_H_
#define SHARK_SIMULATION_H_

#include <vector>
#include <string>

#include "cosmology.h"
#include "options.h"

namespace shark {

class SimulationParameters {

public:
	explicit SimulationParameters(const Options &options);

	float volume = 0;
	float particle_mass = 0;
	float lbox  =0;

	int min_snapshot = 0;
	int max_snapshot = 0;
	int tot_nsubvols = 0;

	std::string sim_name {};

	std::string tree_files_prefix {"tree."};

	std::map<int,double> redshifts {};


	void load_simulation_tables(const std::string &redshift_file);
};


class Simulation {

public:

	Simulation(SimulationParameters parameters, const CosmologyPtr &cosmology);

	double convert_snapshot_to_age(int s);

private:
	SimulationParameters parameters;
	CosmologyPtr cosmology;


};

}  // namespace shark

#endif // SHARK_SIMULATION_H_
