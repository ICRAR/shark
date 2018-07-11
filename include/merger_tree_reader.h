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
	SURFSReader(const std::string &prefix, const DarkMatterHalosPtr &dark_matter_halos, const SimulationParameters &sim_params, unsigned int threads);

	const std::vector<HaloPtr> read_halos(std::vector<unsigned int> batches);

private:
	std::string prefix;
	DarkMatterHalosPtr dark_matter_halos;
	SimulationParameters simulation_params;
	unsigned int threads;

	const std::vector<HaloPtr> read_halos(unsigned int batch);
	const std::vector<SubhaloPtr> read_subhalos(unsigned int batch);
	const std::string get_filename(int batch);
	const std::vector<Subhalo> read_subhalos_batch(int batch);


};

}  // namespace shark



#endif /* INCLUDE_MERGER_TREE_READER_H_ */
