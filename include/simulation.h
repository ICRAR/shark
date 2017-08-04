//
// The simulation parameters used as input for shark
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

#ifndef SHARK_SIMULATION_H_
#define SHARK_SIMULATION_H_

#include <vector>
#include <memory>
#include <string>

#include "cosmology.h"
#include "options.h"

namespace shark {

class SimulationParameters : public Options {

public:
	SimulationParameters(const std::string &filename);

	float volume;
	float particle_mass;

	int min_snapshot;
	int max_snapshot;

	std::string tree_files_prefix;

	std::map<int,double> redshifts;


	void load_simulation_tables(const std::string &redshift_file);
};


class Simulation {

public:

	Simulation(SimulationParameters parameters, std::shared_ptr<Cosmology> cosmology);

	double convert_snapshot_to_age(int s);

private:
	SimulationParameters parameters;
	std::shared_ptr<Cosmology> cosmology;


};

}  // namespace shark

#endif // SHARK_SIMULATION_H_
