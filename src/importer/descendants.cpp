//
// Descendants-related class implementations
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

#include "exceptions.h"
#include "utils.h"
#include "hdf5/reader.h"
#include "importer/descendants.h"

using namespace std;

namespace shark {

namespace importer {

//
// BaseDescendantReader methods follow
//
DescendantReader::DescendantReader(const string &filename) :
	filename(filename)
{
	if ( filename.size() == 0 ) {
		throw invalid_argument("Descendants file has no value");
	}
}

DescendantReader::~DescendantReader()
{
	// no-op
}

//
// AsciiDescendantReader methods follow
//
AsciiDescendantReader::AsciiDescendantReader(const string &filename) :
	DescendantReader(filename)
{
	// no-op
}

vector<descendants_data_t> AsciiDescendantReader::read_whole()
{
	string line;
	ifstream descendants_f = open_file(filename);

	// The first line tells us how many descendants there are
	unsigned int nhalos;
	getline(descendants_f, line);
	istringstream linestream(line);
	linestream >> nhalos;

	// continue reading the rest and check if we'll need a final sorting
	vector<descendants_data_t> descendants;
	descendants.reserve(nhalos);
	while ( getline(descendants_f, line) ) {
		descendants_data_t desc;
		istringstream linestream(line);
		linestream >> desc.halo_id >> desc.halo_snapshot >> desc.descendant_id >> desc.descendant_snapshot;
		descendants.push_back(move(desc));
	}

	return descendants;
}

//
// HDF5DescendantReader methods follow
//

HDF5DescendantReader::HDF5DescendantReader(const string &filename) :
	DescendantReader(filename)
{
	// no-op
}

vector<descendants_data_t> HDF5DescendantReader::read_whole()
{
	hdf5::Reader reader(filename);

	vector<long> halo_ids = reader.read_dataset_v<long>("Halo_IDs");
	vector<int> halo_snaps = reader.read_dataset_v<int>("Halo_Snapshots");
	vector<long> desc_ids = reader.read_dataset_v<long>("Descendant_IDs");
	vector<int> desc_snaps = reader.read_dataset_v<int>("Descendant_Snapshots");

	// Check that all sizes are the same
	auto size = halo_ids.size();
	if ( halo_snaps.size() != size ) {
		ostringstream os;
		os << "Halo_Snapshots length != Halo_IDs length: " << size << " != " << halo_snaps.size();
		throw invalid_data(os.str());
	}
	if ( desc_ids.size() != size ) {
		ostringstream os;
		os << "Descendant_IDs length != Halo_IDs length: " << size << " != " << desc_ids.size();
		throw invalid_data(os.str());
	}
	if ( desc_snaps.size() != size ) {
		ostringstream os;
		os << "Descendant_Snapshots length != Halo_IDs length: " << size << " != " << desc_snaps.size();
		throw invalid_data(os.str());
	}

	vector<descendants_data_t> descendants;
	descendants.reserve(size);
	for(unsigned int i=0; i!=size; i++) {
		descendants_data_t desc = {
			.halo_id = halo_ids[i],
			.halo_snapshot = halo_snaps[i],
			.descendant_id = desc_ids[i],
			.descendant_snapshot = desc_snaps[i]
		};
		descendants.push_back(move(desc));
	}

	return descendants;
}

}  // namespace trees

}  // namespace shark