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
#include <iterator>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "exceptions.h"
#include "subhalo.h"
#include "utils.h"
#include "hdf5/io/reader.h"

namespace shark {

namespace importer {

VELOCIraptorReader::VELOCIraptorReader(std::shared_ptr<DescendantReader> &reader, const std::string &trees_dir) :
	trees_dir(trees_dir)
{
	if (trees_dir.empty()) {
		throw invalid_argument("Trees dir has no value");
	}

	// read all descendants info and put them into our internal map
	// for quick lookup
	auto descendants = reader->read_whole();

	for(auto &&descendant: descendants) {
		descendants_data[descendant.halo_id] = descendant;
	}
}

const std::string VELOCIraptorReader::get_filename(int snapshot, unsigned int batch)
{
	std::ostringstream os;
	os << trees_dir << "/snapshot_" << snapshot << ".VELOCIraptor.hdf.properties." << batch;
	return os.str();
}

std::vector<Subhalo> VELOCIraptorReader::read_subhalos(int snapshot)
{
	unsigned int nbatches;

	// Give the hdf5 reader a scope so it gets destroyed quickly after use
	{
		hdf5::Reader batchfile_0(get_filename(snapshot, 0));
		nbatches = batchfile_0.read_dataset<unsigned int>("Num_of_files");
	}

	// read all batches and add them up to a single vector,
	// which we then return
	std::vector<Subhalo> subhalos;
	for(unsigned int batch=0; batch != nbatches; batch++) {
		auto batch_subhalos = read_subhalos_batch(snapshot, batch);
		std::move(batch_subhalos.begin(), batch_subhalos.end(), std::back_inserter(subhalos));
	}
	return subhalos;
}

std::vector<Subhalo> VELOCIraptorReader::read_subhalos_batch(int snapshot, unsigned int batch)
{
	hdf5::Reader batch_file(get_filename(snapshot, batch));
	auto n_subhalos = batch_file.read_dataset<unsigned int>("Num_of_groups");
	if (n_subhalos == 0) {
		return {};
	}

	std::vector<double> inx = batch_file.read_dataset_v<double>("Xc");
	std::vector<double> iny = batch_file.read_dataset_v<double>("Yc");
	std::vector<double> inz = batch_file.read_dataset_v<double>("Zc");
	std::vector<double> invx = batch_file.read_dataset_v<double>("VXc");
	std::vector<double> invy = batch_file.read_dataset_v<double>("VYc");
	std::vector<double> invz = batch_file.read_dataset_v<double>("VZc");
	std::vector<double> inbmass = batch_file.read_dataset_v<double>("Mass_tot");
	std::vector<Subhalo::id_t> inhalo = batch_file.read_dataset_v<Subhalo::id_t>("ID");
	std::vector<Subhalo::id_t> hhalo = batch_file.read_dataset_v<Subhalo::id_t>("hostHaloID");
	std::vector<Subhalo::id_t> ID_mbp = batch_file.read_dataset_v<Subhalo::id_t>("ID_mbp");
	std::vector<double> invmax = batch_file.read_dataset_v<double>("Vmax");
	std::vector<double> r2 = batch_file.read_dataset_v<double>("R_HalfMass");
	std::vector<double> v_dispersion = batch_file.read_dataset_v<double>("sigV");
	std::vector<double> Lx = batch_file.read_dataset_v<double>("Lx");
	std::vector<double> Ly = batch_file.read_dataset_v<double>("Ly");
	std::vector<double> Lz = batch_file.read_dataset_v<double>("Lz");

	std::vector<Subhalo> subhalos;
	for(unsigned int i=0; i!=n_subhalos; i++) {

		subhalos.emplace_back(inhalo[i], snapshot);
		auto &subhalo = subhalos.back();

		//
		// TODO: here we assign properties, etc
		//

		// Find the corresponding descendant information
		auto it = descendants_data.find(subhalo.id);
		if( it == descendants_data.end() ) {
			std::ostringstream os;
			os << "No data could be found in the descendants file for subhalo id=" << subhalo.id;
			throw invalid_data(os.str());
		}

		// Fill the missing bits
		auto desc = it->second;
		subhalo.descendant_id = desc.descendant_id;
		subhalo.descendant_snapshot = desc.descendant_snapshot;
	}

	return subhalos;
}


}  // namespace importer

}  // namespace shark
