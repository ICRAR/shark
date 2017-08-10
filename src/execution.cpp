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

namespace shark {


ExecutionParameters::ExecutionParameters(const Options &options) :
	output_snapshots(),
	output_format(),
	output_directory(),
	simulation_batches(),
	skip_missing_descendants(false)
{
	options.load("execution.output_snapshots", output_snapshots);
	options.load("execution.output_format", output_format);
	options.load("execution.output_directory", output_directory);
	options.load("execution.simulation_batches", simulation_batches, true);
	options.load("execution.skip_missing_descendants", skip_missing_descendants);
}

}
