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
 */

#ifndef INCLUDE_MERGER_TREE_READER_H_
#define INCLUDE_MERGER_TREE_READER_H_

#include <map>
#include <memory>

#include "components.h"
#include "dark_matter_halos.h"
#include "simulation.h"

namespace shark {

class SURFSReader {

public:

	/**
	 * Constructor.
	 *
	 * @param trees_dir Directory where all tree files are located
	 */
	SURFSReader(const std::string &prefix, DarkMatterHalosPtr dark_matter_halos, SimulationParameters simulation_params, unsigned int threads);

	const std::vector<HaloPtr> read_halos(std::vector<unsigned int> batches);

private:
	std::string prefix;
	DarkMatterHalosPtr dark_matter_halos;
	SimulationParameters simulation_params;
	unsigned int threads;

	const std::vector<HaloPtr> read_halos(unsigned int batch);
	const std::vector<SubhaloPtr> read_subhalos(unsigned int batch);
	const std::string get_filename(unsigned int batch);

};

}  // namespace shark



#endif /* INCLUDE_MERGER_TREE_READER_H_ */
