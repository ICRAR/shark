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

namespace shark {

class SURFSReader {

public:

	/**
	 * Constructor.
	 *
	 * @param trees_dir Directory where all tree files are located
	 */
	SURFSReader(const std::string &prefix);

	const std::vector<HaloPtr> read_halos(std::vector<int> batches);
	const std::vector<HaloPtr> read_halos(int batch);

private:

	std::string prefix;
	const std::string get_filename(int batch);
	const std::vector<Subhalo> read_subhalos_batch(int batch);
};

}  // namespace shark



#endif /* INCLUDE_MERGER_TREE_READER_H_ */
