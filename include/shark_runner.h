//
// Main shark runner class definition
//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2018
// Copyright by UWA (in the framework of the ICRAR)
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307  USA
//

#ifndef SHARK_SHARK_RUNNER_H
#define SHARK_SHARK_RUNNER_H

#include <memory>

#include "options.h"

namespace shark {

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
	SharkRunner(const Options &options, unsigned int threads);
	~SharkRunner();
	void run();

private:
	class impl;
	std::unique_ptr<impl> pimpl;

};

} // namespace shark

#endif // SHARK_SHARK_RUNNER_H