/*
 * execution.h
 *
 *  Created on: 14Jun.,2017
 *      Author: clagos
 */


/*
 * execution.h
 *
 *  Created on: 14Jun.,2017
 *      Author: clagos
 */

#ifndef SHARK_EXECUTION_H_
#define SHARK_EXECUTION_H_

#include <set>
#include <string>
#include <vector>

#include "options.h"

namespace shark {

class ExecutionParameters {

public:
	ExecutionParameters(const Options &options);

	std::vector<int> output_snapshots {};
	Options::file_format_t output_format = Options::HDF5;
	std::string output_directory {};
	std::string name_model {};
	std::vector<unsigned int> simulation_batches {};

	bool skip_missing_descendants = true;
	bool warn_on_missing_descendants = true;

	/**
	 * Parameters of sf histories:
	 * output_sf_histories: boolean parameter set to true if the user wants the star formation histories to be output.
	 * snapshots_sf_histories: vector of int with the snapshots the user wants the star formation histories output at.
	 */
	bool output_sf_histories = false;
	std::vector<int> snapshots_sf_histories {};

	float ode_solver_precision = 0;
};

} // namespace shark

#endif // SHARK_EXECUTION_H_
