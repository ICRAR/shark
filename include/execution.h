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

	std::vector<int> output_snapshots;
	std::string output_format;
	std::string output_directory;
	std::vector<int> simulation_batches;

	bool skip_missing_descendants;

	float ode_solver_precision;
};

} // namespace shark

#endif // SHARK_EXECUTION_H_
