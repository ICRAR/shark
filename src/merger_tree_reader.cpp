/*
 * merger_tree_reader.cpp
 *
 *  Created on: 31Jul.,2017
 *      Author: clagos
 */

#include <array>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <tuple>

#include "merger_tree_reader.h"
#include "logging.h"
#include "exceptions.h"
#include "utils.h"
#include "hdf5/reader.h"

using namespace std;

namespace shark {

SURFSReader::SURFSReader(const std::string &prefix) :
	prefix(prefix)
{
	if ( prefix.size() == 0 ) {
		throw invalid_argument("Trees dir has no value");
	}

}

const string SURFSReader::get_filename(int batch)
{
	ostringstream os;
	os << prefix << "." << batch << ".hdf5";
	return os.str();
}

const std::vector<std::shared_ptr<Halo>> SURFSReader::read_halos(std::vector<int> batches)
{

	// Check that batch numbers are within boundaries
	// (supposing that the file for batch 0 always exists)
	unsigned int nbatches;
	{
		auto batch0_fname = get_filename(0);
		LOG(debug) << "Opening " << batch0_fname << " for reading";
		hdf5::Reader batchfile_0(batch0_fname);
		nbatches = batchfile_0.read_attribute<unsigned int>("fileInfo/numberOfFiles");
	}

	for(auto batch: batches) {
		if (batch >= nbatches) {
			std::ostringstream os;
			os << "Batch is greater than numberOfFile specified in " << get_filename(0);
			os << ": " << batch << " > " << nbatches;
			throw invalid_argument(os.str());
		}
	}

	// Read halos for each batch, accumulate and return
	std::vector<std::shared_ptr<Halo>> all_halos;
	for(auto batch: batches) {
		LOG(info) << "Reading file for batch " << batch;
		auto halos_batch = read_halos(batch);
		all_halos.insert(all_halos.end(), halos_batch.begin(), halos_batch.end());
	}

	return all_halos;
}

const std::vector<std::shared_ptr<Halo>> SURFSReader::read_halos(int batch)
{
	const auto fname = get_filename(batch);
	hdf5::Reader batch_file(fname);

	//Read position and velocities first.
	vector<double> position = batch_file.read_dataset_v_2<double>("haloTrees/position");
	vector<double> velocity = batch_file.read_dataset_v_2<double>("haloTrees/velocity");

	//Read mass, circular velocity and angular momentum.
	vector<double> Mvir = batch_file.read_dataset_v<double>("haloTrees/nodeMass");
	vector<double> Vcirc = batch_file.read_dataset_v<double>("haloTrees/maximumCircularVelocity");
	vector<double> L = batch_file.read_dataset_v_2<double>("haloTrees/angularMomentum");

	//Read indices and the snapshot number at which the subhalo lives.
	vector<int> snap = batch_file.read_dataset_v<int>("haloTrees/snapshotNumber");
	vector<Subhalo::id_t> nodeIndex = batch_file.read_dataset_v<Subhalo::id_t>("haloTrees/nodeIndex");
	vector<Subhalo::id_t> descIndex = batch_file.read_dataset_v<Subhalo::id_t>("haloTrees/descendantIndex");
	vector<Halo::id_t> hostIndex = batch_file.read_dataset_v<Halo::id_t>("haloTrees/hostIndex");
	vector<Halo::id_t> descHost = batch_file.read_dataset_v<Halo::id_t>("haloTrees/descendantHost");

	//Read properties that characterise the position of the subhalo inside the halo.
	vector<int> IsMain = batch_file.read_dataset_v<int>("haloTrees/isMainProgenitor");
	vector<int> IsCentre = batch_file.read_dataset_v<int>("haloTrees/isDHaloCentre");


	int n_subhalos = Mvir.size();

	if ( !n_subhalos ) {
		return {};
	}

	vector<std::shared_ptr<Subhalo>> subhalos;
	for(unsigned int i=0; i!=n_subhalos; i++) {

		std::shared_ptr<Subhalo> subhalo = std::make_shared<Subhalo>();

		//Assign indices.
		subhalo->id = nodeIndex[i];
		subhalo->descendant_id = descIndex[i];
		subhalo->descendant_halo_id = descHost[i];
		subhalo->haloID = hostIndex[i];
		subhalo->snapshot = snap[i];

		//Determine if subhalo is centre of Dhalo. This is done using IsMainProgenitor, as this is the halo
		//that is found by dhalos to be the centre (although not necessarily the most massive).
		if(IsMain[i] == 1) {
			subhalo->subhalo_type = Subhalo::CENTRAL;
		}
		else {
			subhalo->subhalo_type = Subhalo::SATELLITE;
		}

		//Assign mass.
		subhalo->Mvir = Mvir[i];

		//Assign position
		subhalo->position.x = position[3 * i];
		subhalo->position.y = position[3 * i + 1];
		subhalo->position.z = position[3 * i + 2];

		//Assign velocity
		subhalo->velocity.x = velocity[3 * i];
		subhalo->velocity.y = velocity[3 * i + 1];
		subhalo->velocity.z = velocity[3 * i + 2];

		//Assign angular momentum
		subhalo->L.x = L[3 * i];
		subhalo->L.y = L[3 * i + 1];
		subhalo->L.z = L[3 * i + 2];

		subhalo->Vcirc = Vcirc[i];

		// Done, save it now
		subhalos.push_back(std::move(subhalo));
	}

	LOG(info) << "Read " << subhalos.size() << " Subhalos from " << fname;

	// Sort subhalos by host index (which intrinsically sorts them by snapshot
	// since host indices numbers are prefixed with the snapshot number)
	std::sort(subhalos.begin(), subhalos.end(), [](const std::shared_ptr<Subhalo> &lhs, const std::shared_ptr<Subhalo> &rhs) {
		return lhs->haloID < rhs->haloID;
	});

	//Assign properties to halos.
	std::shared_ptr<Halo> current_halo;
	std::vector<std::shared_ptr<Halo>> halos;
	Halo::id_t last_halo_id = -1;
	for(const auto &subhalo: subhalos) {

		auto halo_id = subhalo->haloID;
		if (halo_id != last_halo_id) {
			last_halo_id = halo_id;
			current_halo = std::make_shared<Halo>(halo_id, subhalo->snapshot);
			halos.push_back(current_halo);
		}

		current_halo->add_subhalo(subhalo);
		subhalo->host_halo = current_halo;
	}

	LOG(info) << "Created " << halos.size() << " Halos from these Subhalos";

	for(const auto &halo: halos) {
		// calculate vvir
	}

	return halos;
}

}  // namespace shark
