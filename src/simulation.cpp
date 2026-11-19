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
 */

#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <tuple>

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
	options.load("simulation.hydrorun", hydrorun, false);

	load_simulation_tables(redshift_file);

}

void SimulationParameters::load_simulation_tables(const std::string &redshift_file)
{
	LOG(debug) << "Reading table " << redshift_file ;

	std::ifstream f = open_file(redshift_file);
	std::string line;
	while ( std::getline(f, line) ) {

		trim(line);
		if (empty_or_comment(line)) {
			continue;
		}

		int s;
		double r;

		std::istringstream iss(line);
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

Simulation::Simulation(SimulationParameters parameters, CosmologyPtr cosmology) :
	parameters(std::move(parameters)),
	cosmology(std::move(cosmology))
{
	// no-op
}

double Simulation::convert_snapshot_to_age(int s){

	return cosmology->convert_redshift_to_age(parameters.redshifts[s]);

}

} // namespace shark
