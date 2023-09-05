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

/**
 * @file
 *
 * Main shark runner class definition
 */

#ifndef SHARK_SHARK_RUNNER_H
#define SHARK_SHARK_RUNNER_H

#include <memory>

namespace shark {

// Forward declaration to avoid including options.h
class Options;

/**
 * The main driver of a shark instance run.
 *
 * A shark runner is given a set of options, and a number of threads to run.
 * Using only this information it reads the initial data, creates merger trees,
 * and evolves the galaxies within.
 *
 * This particular implementation is separated using the pimpl pattern to hide
 * most of the internals to client applications (the main.cpp module in our case),
 * which need only to know about the Options class.
 */
class SharkRunner {

public:
	/**
	 * Constructor
	 *
	 * If there is any missing or invalid option this constructor will throw
	 * and exception.
	 *
	 * @param options The set of options used to run shark
	 * @param threads The number of threads used to run shark
	 */
	SharkRunner(const Options &options, unsigned int threads);
	~SharkRunner();

	/// Run shark until completion
	void run();

	/// Report total execution times
	void report_total_times();

private:
	class impl;
	std::unique_ptr<impl> pimpl;

};

} // namespace shark

#endif // SHARK_SHARK_RUNNER_H