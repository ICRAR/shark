/*
 * execution.h
 *
 *  Created on: 14Jun.,2017
 *      Author: clagos
 */




#include <vector>
#include <string>

#include "options.h"

namespace shark {

class ExecutionParameters : public Options {

public:
	ExecutionParameters(const std::string &filename);

	std::vector<double> output_snapshots;
	char output_format;
	char output_directory;

};
}
