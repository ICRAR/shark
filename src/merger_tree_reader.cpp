/*
 * merger_tree_reader.cpp
 *
 *  Created on: 31Jul.,2017
 *      Author: clagos
 */

#include <array>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <tuple>

#include "config.h"
#include "dark_matter_halos.h"
#include "exceptions.h"
#include "logging.h"
#include "merger_tree_reader.h"
#include "simulation.h"
#include "timer.h"
#include "utils.h"
#include "hdf5/reader.h"

#ifdef SHARK_OPENMP
#include <omp.h>
#endif


using namespace std;

namespace shark {

SURFSReader::SURFSReader(const std::string &prefix, unsigned int threads) :
	prefix(prefix), threads(threads)
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

const std::vector<HaloPtr> SURFSReader::read_halos(std::vector<unsigned int> batches, DarkMatterHalos &darkmatterhalos, SimulationParameters &sim_params)
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
	std::vector<HaloPtr> all_halos;
	for(auto batch: batches) {
		LOG(info) << "Reading file for batch " << batch;
		auto halos_batch = read_halos(batch, darkmatterhalos, sim_params);
		all_halos.reserve(all_halos.size() + halos_batch.size());
		all_halos.insert(all_halos.end(), halos_batch.begin(), halos_batch.end());
	}

	return all_halos;
}

const std::vector<SubhaloPtr> SURFSReader::read_subhalos(unsigned int batch, DarkMatterHalos &darkmatterhalos, SimulationParameters &sim_params)
{
	Timer t;
	const auto fname = get_filename(batch);
	hdf5::Reader batch_file(fname);

	//Read position and velocities first.
	vector<float> position = batch_file.read_dataset_v_2<float>("haloTrees/position");
	vector<float> velocity = batch_file.read_dataset_v_2<float>("haloTrees/velocity");

	//Read mass, circular velocity and angular momentum.
	vector<float> Mvir = batch_file.read_dataset_v<float>("haloTrees/nodeMass");
	vector<float> Vcirc = batch_file.read_dataset_v<float>("haloTrees/maximumCircularVelocity");
	vector<float> L = batch_file.read_dataset_v_2<float>("haloTrees/angularMomentum");
	vector<float> cnfw = batch_file.read_dataset_v<float>("haloTrees/cnfw");
	vector<float> lambda = batch_file.read_dataset_v<float>("haloTrees/lambda");
	vector<float> vvir = batch_file.read_dataset_v<float>("haloTrees/Vvir");

	//Read indices and the snapshot number at which the subhalo lives.
	vector<int> snap = batch_file.read_dataset_v<int>("haloTrees/snapshotNumber");
	vector<Subhalo::id_t> nodeIndex = batch_file.read_dataset_v<Subhalo::id_t>("haloTrees/nodeIndex");
	vector<Subhalo::id_t> descIndex = batch_file.read_dataset_v<Subhalo::id_t>("haloTrees/descendantIndex");
	vector<Halo::id_t> hostIndex = batch_file.read_dataset_v<Halo::id_t>("haloTrees/hostIndex");
	vector<Halo::id_t> descHost = batch_file.read_dataset_v<Halo::id_t>("haloTrees/descendantHost");

	//Read properties that characterise the position of the subhalo inside the halo.
	vector<int> IsMain = batch_file.read_dataset_v<int>("haloTrees/isMainProgenitor");
	vector<int> IsCentre = batch_file.read_dataset_v<int>("haloTrees/isDHaloCentre");
	vector<int> IsInterpolated = batch_file.read_dataset_v<int>("haloTrees/isInterpolated");

	auto n_subhalos = Mvir.size();
	LOG(info) << "Read raw data of " << n_subhalos << " subhalos from " << fname << " in " << t;
	if ( !n_subhalos ) {
		return {};
	}

	std::ostringstream os;
	os << "File " << fname << " has " << n_subhalos << " subhalos. ";
	os << "After reading we should be using ~" << memory_amount(n_subhalos * sizeof(Subhalo)) << " of memory";
	LOG(info) << os.str();

	t = Timer();
	vector<vector<SubhaloPtr>> t_subhalos(threads);

#ifdef SHARK_OPENMP
	#pragma omp parallel for num_threads(threads) schedule(static, 10000)
#endif
	for(unsigned int i=0; i < n_subhalos; i++) {

		auto subhalo = std::make_shared<Subhalo>();

		// Subhalo and Halo index, snapshot
		subhalo->id = nodeIndex[i];
		subhalo->haloID = hostIndex[i];
		subhalo->snapshot = snap[i];

		// Descendant information. -1 means that the Subhalo has no descendant
		auto descendant_id = descIndex[i];
		if (descendant_id == -1) {
			subhalo->has_descendant = false;
		}
		else {
			subhalo->has_descendant = true;
			subhalo->descendant_id = descendant_id;
			subhalo->descendant_halo_id = descHost[i];
		}

		//Assign main progenitor flags.
		if(IsMain[i] == 1){
			subhalo->main_progenitor = true;
		}

		//Assign interpolated subhalo flags.
		if(IsInterpolated[i] == 1){
			subhalo->IsInterpolated = true;
		}

		//Make all subhalos satellite, because once we construct the merger tree we will find the main branch.
		subhalo->subhalo_type = Subhalo::SATELLITE;

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

		//Assign specific angular momentum
		subhalo->L.x = L[3 * i];
		subhalo->L.y = L[3 * i + 1];
		subhalo->L.z = L[3 * i + 2];

		subhalo->Vcirc = Vcirc[i];

		subhalo->concentration = darkmatterhalos.nfw_concentration(subhalo->Mvir,sim_params.redshifts[subhalo->snapshot]);

		if (subhalo->concentration < 1) {
			throw invalid_argument("concentration is <1, cannot continue. Please check input catalogue");
		}

		subhalo->lambda = darkmatterhalos.halo_lambda(lambda[i], sim_params.redshifts[subhalo->snapshot]);

		// Calculate virial velocity from the virial mass and redshift.
		subhalo->Vvir = darkmatterhalos.halo_virial_velocity(subhalo->Mvir, sim_params.redshifts[subhalo->snapshot]);

		// Done, save it now
		t_subhalos[omp_get_thread_num()].emplace_back(std::move(subhalo));
	}

	vector<SubhaloPtr> subhalos;
	if (threads == 0) {
		subhalos = std::move(t_subhalos[0]);
	}
	else {
		subhalos.reserve(n_subhalos);
		for (auto i = 0; i != threads; i++) {
			subhalos.insert(subhalos.end(), t_subhalos[i].begin(), t_subhalos[i].end());
		}
	}

	LOG(info) << "Created " << subhalos.size() << " Subhalos from " << fname << " in " << t;
	return subhalos;
}

const std::vector<HaloPtr> SURFSReader::read_halos(unsigned int batch, DarkMatterHalos &darkmatterhalos, SimulationParameters &sim_params)
{

	std::vector<SubhaloPtr> subhalos = read_subhalos(batch, darkmatterhalos, sim_params);

	// Sort subhalos by host index (which intrinsically sorts them by snapshot
	// since host indices numbers are prefixed with the snapshot number)
	std::sort(subhalos.begin(), subhalos.end(), [](const SubhaloPtr &lhs, const SubhaloPtr &rhs) {
		return lhs->haloID < rhs->haloID;
	});
	LOG(info) << "Sorted subhalos by haloID, creating Halos now";

	// Create and assign Halos
	HaloPtr halo;
	std::vector<HaloPtr> halos;
	Halo::id_t last_halo_id = -1;
	for(const auto &subhalo: subhalos) {

		auto halo_id = subhalo->haloID;
		if (halo_id != last_halo_id) {
			if (halo) {
				halos.emplace_back(std::move(halo));
			}
			last_halo_id = halo_id;
			halo = std::make_shared<Halo>(halo_id, subhalo->snapshot);
		}

		LOG(trace) << "Adding " << subhalo << " to " << halo;
		subhalo->host_halo = halo;
		halo->add_subhalo(std::move(subhalo));
	}
	subhalos.clear();

	std::ostringstream os;
	os << "Created " << halos.size() << " Halos from these Subhalos. ";
	os << "This should take another ~" << memory_amount(halos.size() * sizeof(Halo)) << " of memory";
	LOG(info) << os.str();

	// Calculate halos' vvir and concentration
	Timer t;
#ifdef SHARK_OPENMP
	#pragma omp parallel for num_threads(threads) schedule(dynamic, 10000)
#endif
	for(auto it = halos.begin(); it < halos.end(); it++) {
		const auto &halo = *it;
		halo->Vvir = darkmatterhalos.halo_virial_velocity(halo->Mvir, sim_params.redshifts[halo->snapshot]);
		halo->concentration = darkmatterhalos.nfw_concentration(halo->Mvir,sim_params.redshifts[halo->snapshot]);
	}
	LOG(info) << "Calculated Vvir and concentration for new Halos in " << t;

	return halos;
}

}  // namespace shark
