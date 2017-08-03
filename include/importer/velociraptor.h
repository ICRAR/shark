//
// VELOCiraptor reader class definition
//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2017
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

#ifndef SHARK_IMPORTER_VELOCIRAPTOR
#define SHARK_IMPORTER_VELOCIRAPTOR

#include <map>
#include <memory>

#include "importer/descendants.h"
#include "importer/reader.h"

namespace shark {

namespace importer {

class VELOCIraptorReader : public Reader {

public:

	/**
	 * Constructor.
	 *
	 * @param trees_dir Directory where all tree files are located
	 */
	VELOCIraptorReader(std::shared_ptr<DescendantReader> &descendant_reader, const std::string &trees_dir);

	virtual const std::vector<Subhalo> read_subhalos(int snapshot) override;

private:

	// Indexed by halo_id
	std::map<long, descendants_data_t> descendants_data;

	std::string descendants_file;
	std::string trees_dir;

	const std::string get_filename(int snapshot, int batch);
	const std::vector<Subhalo> read_subhalos_batch(int snapshot, int batch);
};

}  // namespace importer

}  // namespace shark

#endif // SHARK_IMPORTER_VELOCIRAPTOR
