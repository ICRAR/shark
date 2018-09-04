//
// Execution-related unit tests
//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2018
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

#include <cxxtest/TestSuite.h>

#include "execution.h"

using namespace shark;

class TestSubhalos : public CxxTest::TestSuite
{

private:

	void assert_output_snapshots(const std::string &snapshots, std::set<int> expected_snapshots, int expected_last_snapshot)
	{
		Options opts {};
		opts.add("execution.output_format = hdf5");
		opts.add("execution.output_directory = .");
		opts.add("execution.simulation_batches = 0");
		opts.add("execution.ode_solver_precision = 0.5");
		opts.add("execution.name_model = test");
		opts.add(std::string("execution.output_snapshots = ") + snapshots);
		ExecutionParameters params {opts};

		TS_ASSERT_EQUALS(params.output_snapshots, expected_snapshots);
		TS_ASSERT_EQUALS(params.last_output_snapshot(), expected_last_snapshot);
		for (int s: expected_snapshots) {
			TS_ASSERT(params.output_snapshot(s));
		}
		TS_ASSERT(!params.output_snapshot(params.last_output_snapshot() + 1));
	}

public:

	void test_output_snapshots()
	{
		assert_output_snapshots("199", {199}, 199);
		assert_output_snapshots("199 199", {199}, 199);
		assert_output_snapshots("199 198", {199, 198}, 199);
		assert_output_snapshots("199 0", {0, 199}, 199);
		assert_output_snapshots("199 0 199", {0, 199}, 199);
		assert_output_snapshots("0 199", {0, 199}, 199);
	}
};