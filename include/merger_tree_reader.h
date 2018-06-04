/*
 * merger_tree_reader.h
 *
 *  Created on: 31Jul.,2017
 *      Author: clagos
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
