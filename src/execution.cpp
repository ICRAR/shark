/*
 * execution.cpp
 *
 *  Created on: 14Jun.,2017
 *      Author: clagos
 */

#include <cmath>
#include <fstream>
#include <map>
#include <tuple>

#include "execution.h"
#include "logging.h"

namespace shark {


ExecutionParameters::ExecutionParameters(const std::string &filename) :
	Options(filename),
	output_snapshots(),
	output_format(),
	output_directory(),
	simulation_batches(),
	skip_missing_descendants(false)
{
	load("execution.output_snapshots", output_snapshots);
	load("execution.output_format", output_format);
	load("execution.output_directory", output_directory);
	load("execution.simulation_batches", simulation_batches, true);
	load("execution.skip_missing_descendants", skip_missing_descendants);
}

}
