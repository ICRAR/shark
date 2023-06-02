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
 * Descendants-related class implementations
 */

#include "exceptions.h"
#include "utils.h"
#include "hdf5/io/reader.h"
#include "importer/descendants.h"

namespace shark {

namespace importer {

//
// BaseDescendantReader methods follow
//
DescendantReader::DescendantReader(const std::string &filename) :
	filename(filename)
{
	if (filename.empty()) {
		throw invalid_argument("Descendants file has no value");
	}
}

DescendantReader::~DescendantReader() = default;

//
// AsciiDescendantReader methods follow
//
AsciiDescendantReader::AsciiDescendantReader(const std::string &filename) :
	DescendantReader(filename)
{
	// no-op
}

std::vector<descendants_data_t> AsciiDescendantReader::read_whole()
{
	std::string line;
	std::ifstream descendants_f = open_file(filename);

	// The first line tells us how many descendants there are
	unsigned int nhalos;
	std::getline(descendants_f, line);
	std::istringstream linestream(line);
	linestream >> nhalos;

	// continue reading the rest and check if we'll need a final sorting
	std::vector<descendants_data_t> descendants;
	descendants.reserve(nhalos);
	while ( getline(descendants_f, line) ) {
		descendants_data_t desc;
		std::istringstream linestream(line);
		linestream >> desc.halo_id >> desc.halo_snapshot >> desc.descendant_id >> desc.descendant_snapshot;
		descendants.push_back(desc);
	}

	return descendants;
}

//
// HDF5DescendantReader methods follow
//

HDF5DescendantReader::HDF5DescendantReader(const std::string &filename) :
	DescendantReader(filename)
{
	// no-op
}

std::vector<descendants_data_t> HDF5DescendantReader::read_whole()
{
	hdf5::Reader reader(filename);

	std::vector<long> halo_ids = reader.read_dataset_v<long>("Halo_IDs");
	std::vector<int> halo_snaps = reader.read_dataset_v<int>("Halo_Snapshots");
	std::vector<long> desc_ids = reader.read_dataset_v<long>("Descendant_IDs");
	std::vector<int> desc_snaps = reader.read_dataset_v<int>("Descendant_Snapshots");

	// Check that all sizes are the same
	auto size = halo_ids.size();
	if ( halo_snaps.size() != size ) {
		std::ostringstream os;
		os << "Halo_Snapshots length != Halo_IDs length: " << size << " != " << halo_snaps.size();
		throw invalid_data(os.str());
	}
	if ( desc_ids.size() != size ) {
		std::ostringstream os;
		os << "Descendant_IDs length != Halo_IDs length: " << size << " != " << desc_ids.size();
		throw invalid_data(os.str());
	}
	if ( desc_snaps.size() != size ) {
		std::ostringstream os;
		os << "Descendant_Snapshots length != Halo_IDs length: " << size << " != " << desc_snaps.size();
		throw invalid_data(os.str());
	}

	std::vector<descendants_data_t> descendants;
	descendants.reserve(size);
	for(unsigned int i=0; i!=size; i++) {
		descendants_data_t desc = {
			halo_ids[i],
			halo_snaps[i],
			desc_ids[i],
			desc_snaps[i]
		};
		descendants.push_back(desc);
	}

	return descendants;
}

}  // namespace importer

}  // namespace shark