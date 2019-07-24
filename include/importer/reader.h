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
 * Base importer reader class definition
 */

#ifndef SHARK_IMPORTER_READER
#define SHARK_IMPORTER_READER

#include <stdexcept>
#include <vector>

#include "components.h"

namespace shark {

namespace importer {

class Reader {

public:

	virtual ~Reader() = 0;

	/**
	 * Reads all the subhalos for a given snapshot
	 *
	 * @param snapshot The snapshot number
	 * @return A list of populated Subhalos
	 */
	virtual std::vector<Subhalo> read_subhalos(int snapshot) = 0;

};

}  // namespace importer

}  // namespace shark

#endif // SHARK_IMPORTER_READER