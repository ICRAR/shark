//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2018
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

#ifndef SHARK_EXECUTION_H_
#define SHARK_EXECUTION_H_

#include <cassert>
#include <ctime>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "components/algorithms.h"
#include "options.h"

namespace shark {

class ExecutionParameters {

public:
	explicit ExecutionParameters(const Options &options);

	std::set<int> output_snapshots;
	Options::file_format_t output_format = Options::HDF5;
	std::string output_directory;
	std::string name_model;
	std::random_device::result_type seed = std::random_device()();
	std::vector<unsigned int> simulation_batches;
	std::time_t starting_time = std::time(nullptr);

	bool output_snapshot(int snapshot);
	int last_output_snapshot();

	template <typename Component>
	std::random_device::result_type get_seed(const Component &component)
	{
		auto id = get_id(component);
		assert(id >= 0);
		return seed + std::random_device::result_type(id);
	}

	bool skip_missing_descendants = true;
	bool warn_on_missing_descendants = true;
	bool ensure_mass_growth = true;
	bool ignore_late_massive_halos  = false;

	/**
	 * Parameters of sf histories:
	 * output_sf_histories: boolean parameter set to true if the user wants the star formation histories to be output.
	 * snapshots_sf_histories: vector of int with the snapshots the user wants the star formation histories output at.
	 */
	bool output_sf_histories = false;
	std::vector<int> snapshots_sf_histories;

	/**
	 * Parameters of BH histories:
	 * output_bh_histories: boolean parameter set to true if the user wants the black hole formation histories to be output.
	 * snapshots_bh_histories: vector of int with the snapshots the user wants the black hole formation histories output at.
	 */

	bool output_bh_histories = false;
	std::vector<int> snapshots_bh_histories;

	float ode_solver_precision = 0;
	int ignore_npart_threshold = 1000;
	float ignore_below_z = 1.0;
};

} // namespace shark

#endif // SHARK_EXECUTION_H_
