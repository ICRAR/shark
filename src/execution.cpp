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


ExecutionParameters::ExecutionParameters(const Options &options)
{
	options.load("execution.output_snapshots", output_snapshots, true);
	options.load("execution.output_format", output_format, true);
	options.load("execution.output_directory", output_directory, true);
	options.load("execution.simulation_batches", simulation_batches, true);
	options.load("execution.skip_missing_descendants", skip_missing_descendants);
	options.load("execution.warn_on_missing_descendants", warn_on_missing_descendants);
        options.load("execution.ensure_mass_growth", ensure_mass_growth);

	options.load("execution.ode_solver_precision", ode_solver_precision, true);
	options.load("execution.name_model", name_model, true);

	options.load("execution.output_sf_histories", output_sf_histories);
	options.load("execution.snapshots_sf_histories", snapshots_sf_histories);
}

}
