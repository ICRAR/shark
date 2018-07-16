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
 * VELOCIraptor reader implementation
 */

#include <importer/velociraptor.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "exceptions.h"
#include "utils.h"
#include "hdf5/reader.h"

using namespace std;

namespace shark {

namespace importer {

VELOCIraptorReader::VELOCIraptorReader(shared_ptr<DescendantReader> &reader, const string &trees_dir) :
	Reader(),
	trees_dir(trees_dir)
{
	if ( trees_dir.size() == 0 ) {
		throw invalid_argument("Trees dir has no value");
	}

	// read all descendants info and put them into our internal map
	// for quick lookup
	auto descendants = reader->read_whole();

	for(auto &&descendant: descendants) {
		descendants_data[descendant.halo_id] = move(descendant);
	}
}

const string VELOCIraptorReader::get_filename(int snapshot, int batch)
{
	ostringstream os;
	os << trees_dir << "/snapshot_" << snapshot << ".VELOCIraptor.hdf.properties." << batch;
	return os.str();
}

const vector<Subhalo> VELOCIraptorReader::read_subhalos(int snapshot)
{
	unsigned int nbatches;

	// Give the hdf5 reader a scope so it gets destroyed quickly after use
	{
		hdf5::Reader batchfile_0(get_filename(snapshot, 0));
		nbatches = batchfile_0.read_dataset<unsigned int>("Num_of_files");
	}

	// read all batches and add them up to a single vector,
	// which we then return
	vector<Subhalo> subhalos;
	for(unsigned int batch=0; batch != nbatches; batch++) {
		auto batch_subhalos = read_subhalos_batch(snapshot, batch);
		subhalos.insert(subhalos.end(), batch_subhalos.begin(), batch_subhalos.end());
	}
	return subhalos;
}

const vector<Subhalo> VELOCIraptorReader::read_subhalos_batch(int snapshot, int batch)
{
	hdf5::Reader batch_file(get_filename(snapshot, batch));
	auto n_subhalos = batch_file.read_dataset<unsigned int>("Num_of_groups");
	if ( !n_subhalos ) {
		return {};
	}

	vector<double> inx = batch_file.read_dataset_v<double>("Xc");
	vector<double> iny = batch_file.read_dataset_v<double>("Yc");
	vector<double> inz = batch_file.read_dataset_v<double>("Zc");
	vector<double> invx = batch_file.read_dataset_v<double>("VXc");
	vector<double> invy = batch_file.read_dataset_v<double>("VYc");
	vector<double> invz = batch_file.read_dataset_v<double>("VZc");
	vector<double> inbmass = batch_file.read_dataset_v<double>("Mass_tot");
	vector<Subhalo::id_t> inhalo = batch_file.read_dataset_v<Subhalo::id_t>("ID");
	vector<Subhalo::id_t> hhalo = batch_file.read_dataset_v<Subhalo::id_t>("hostHaloID");
	vector<Subhalo::id_t> ID_mbp = batch_file.read_dataset_v<Subhalo::id_t>("ID_mbp");
	vector<double> invmax = batch_file.read_dataset_v<double>("Vmax");
	vector<double> r2 = batch_file.read_dataset_v<double>("R_HalfMass");
	vector<double> v_dispersion = batch_file.read_dataset_v<double>("sigV");
	vector<double> Lx = batch_file.read_dataset_v<double>("Lx");
	vector<double> Ly = batch_file.read_dataset_v<double>("Ly");
	vector<double> Lz = batch_file.read_dataset_v<double>("Lz");

	vector<Subhalo> subhalos;
	for(unsigned int i=0; i!=n_subhalos; i++) {

		Subhalo subhalo(inhalo[i], snapshot);

		//
		// TODO: here we assign properties, etc
		//

		// Find the corresponding descendant information
		auto it = descendants_data.find(subhalo.id);
		if( it == descendants_data.end() ) {
			ostringstream os;
			os << "No data could be found in the descendants file for subhalo id=" << subhalo.id;
			throw invalid_data(os.str());
		}

		// Fill the missing bits
		auto desc = it->second;
		subhalo.descendant_id = desc.descendant_id;
		subhalo.descendant_snapshot = desc.descendant_snapshot;

		// Done, save it now
		subhalos.push_back(move(subhalo));
	}

	return subhalos;
}


}  // namespace trees

}  // namespace shark
