//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2017
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

	options.load("execution.name_model", name_model, true);
	options.load("execution.seed", seed);

	options.load("execution.ode_solver_precision", ode_solver_precision, true);

	options.load("execution.ignore_late_massive_halos", ignore_late_massive_halos);
        options.load("execution.ignore_npart_threshold", ignore_npart_threshold);
        options.load("execution.ignore_below_z", ignore_below_z);

	options.load("execution.output_sf_histories", output_sf_histories);
	options.load("execution.snapshots_sf_histories", snapshots_sf_histories);

	options.load("execution.output_bh_histories", output_bh_histories);
	options.load("execution.snapshots_bh_histories", snapshots_bh_histories);
}

bool ExecutionParameters::output_snapshot(int snapshot)
{
	return output_snapshots.find(snapshot) != output_snapshots.end();
}

int ExecutionParameters::last_output_snapshot()
{
	return *output_snapshots.rbegin();
}

} // namespace shark
