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
 * VELOCiraptor reader class definition
 */

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

	std::vector<Subhalo> read_subhalos(int snapshot) override;

private:

	// Indexed by halo_id
	std::map<long, descendants_data_t> descendants_data;

	std::string descendants_file;
	std::string trees_dir;

	const std::string get_filename(int snapshot, unsigned int batch);
	std::vector<Subhalo> read_subhalos_batch(int snapshot, unsigned int batch);
};

}  // namespace importer

}  // namespace shark

#endif // SHARK_IMPORTER_VELOCIRAPTOR
