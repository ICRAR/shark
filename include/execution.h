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

class ExecutionParameters : public Options {

public:
	ExecutionParameters(const std::string &filename);

	std::set<int> output_snapshots;
	std::string output_format;
	std::string output_directory;
	std::vector<int> simulation_batches;
};

} // namespace shark

#endif // SHARK_EXECUTION_H_