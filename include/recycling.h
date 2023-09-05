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
 *
 * Recycling used as input for shark
 */

#ifndef SHARK_RECYCLING_H_
#define SHARK_RECYCLING_H_

#include <vector>
#include <string>

#include "options.h"

namespace shark {

class RecyclingParameters {

public:
	explicit RecyclingParameters(const Options &options);

	double yield = 0;
	double recycle = 0;
	double zsun = 0;

	bool evolving_yield = false;
};

}  // namespace shark

#endif // SHARK_RECYCLING_H_
