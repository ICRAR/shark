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



#include <vector>
#include <string>

#include "options.h"

namespace shark {

class ExecutionParameters : public Options {

public:
	ExecutionParameters(const std::string &filename);

	std::vector<double> output_snapshots;
	std::string output_format;
	std::string output_directory;

};
}


