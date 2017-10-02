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
	output_format(Options::HDF5),
	output_directory(),
	name_model(),
	simulation_batches(),
	skip_missing_descendants(true),
	ode_solver_precision()
{
	options.load("execution.output_snapshots", output_snapshots, true);
	options.load("execution.output_format", output_format, true);
	options.load("execution.output_directory", output_directory, true);
	options.load("execution.simulation_batches", simulation_batches, true);
	options.load("execution.skip_missing_descendants", skip_missing_descendants);
	options.load("execution.ode_solver_precision", ode_solver_precision, true);
	options.load("execution.name_model", name_model, true);

}

}
