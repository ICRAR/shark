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

#include "importer/reader.h"

namespace shark {

namespace importer {

class SURFSReader : public Reader {

public:

	/**
	 * Constructor.
	 *
	 * @param trees_dir Directory where all tree files are located
	 */
	SURFSReader(const std::string &filename);

	virtual const std::vector<Subhalo> read_subhalos(int snapshot) override;

private:

	std::string trees_dir;
	class MergerTreeReader: public Options{


	};

	const std::string get_filename(int snapshot, int batch);
	const std::vector<Subhalo> read_subhalos_batch(int batch);
};

}  // namespace importer

}  // namespace shark



#endif /* INCLUDE_MERGER_TREE_READER_H_ */
