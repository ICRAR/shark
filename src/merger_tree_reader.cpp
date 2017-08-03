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

namespace importer {

SURFSReader::SURFSReader(const std::string &filename) :
	Options(filename),
	Reader(),
	trees_dir(trees_dir)
{
	load("simulation.tree_dir", tree_dir, true);

	if ( trees_dir.size() == 0 ) {
		throw invalid_argument("Trees dir has no value");
	}

}

const string SURFSReader::get_filename(int snapshot, int batch)
{
	ostringstream os;
	os << trees_dir << "." << batch;
	return os.str();
}

const vector<Subhalo> SURFSReader::read_subhalos(int snapshot)
{
	unsigned int nbatches;
	unsigned int total_subhalos;

	// Give the hdf5 reader a scope so it gets destroyed quickly after use
	{
		hdf5::Reader batchfile_0(get_filename(snapshot, 0));
		nbatches = batchfile_0.read_attribute<unsigned int>("fileInfo/numberOfFile");
		//total_subhalos = batchfile_0.read_dataset<unsigned int>("Total_num_of_groups");
	}

	// read all batches and add them up to a single vector,
	// which we then return
	vector<Subhalo> subhalos;
	for(unsigned int batch=0; batch != nbatches; batch++) {
		auto batch_subhalos = read_subhalos_batch(batch);
		subhalos.insert(subhalos.end(), batch_subhalos.begin(), batch_subhalos.end());
	}
	return subhalos;
}

const vector<Subhalo> SURFSReader::read_subhalos_batch(int batch)
{
	hdf5::Reader batch_file(get_filename(batch));


	//Read position and velocities first.
	vector<double> position = batch_file.read_dataset_v<double>("haloTrees/position");
	vector<double> velocity = batch_file.read_dataset_v<double>("haloTrees/velocity");

	//Read mass, circular velocity and angular momentum.
	vector<double> Mvir 	= batch_file.read_dataset_v<double>("haloTrees/nodeMass");
	vector<double> Vcirc 	= batch_file.read_dataset_v<double>("haloTrees/MaximumCircularVelocity");
	vector<double> L		= batch_file.read_dataset_v<double>("haloTrees/angularMomentum");

	//Read indices and the snapshot number at which the subhalo lives.
	vector<Subhalo::id_t> nodeIndex = batch_file.read_dataset_v<double>("haloTrees/nodeIndex");
	vector<Subhalo::id_t> hostIndex = batch_file.read_dataset_v<double>("haloTrees/hostIndex");
	vector<Subhalo::id_t> descIndex = batch_file.read_dataset_v<double>("haloTrees/descendantIndex");
	vector<double> snap 			= batch_file.read_dataset_v<double>("haloTrees/snapshotNumber");
	vector<Subhalo::id_t> descHost	= batch_file.read_dataset_v<double>("haloTrees/descendantHost");

	//Read properties that characterise the position of the subhalo inside the halo.
	vector<int> IsMain		= batch_file.read_dataset_v<int>("haloTrees/isMainProgenitor");
	vector<int> IsCentre	= batch_file.read_dataset_v<int>("haloTrees/isDHaloCentre");


	int n_subhalos = Mvir.size();

	if ( !n_subhalos ) {
		return {};
	}

	vector<Subhalo> subhalos;
	for(unsigned int i=0; i!=n_subhalos; i++) {

		Subhalo subhalo;

		//Assign indices.
		subhalo.id = nodeIndex[i];
		subhalo.descendant_id = descIndex[i];
		subhalo.haloID = hostIndex[i];
		subhalo.snapshot = snap[i];


		//Determine if subhalo is centre of Dhalo.
		if(IsCentre[i] == 1) subhalo.subhalo_type = 0;

		//Assign mass.
		subhalo.Mvir = Mvir[i];

		//Assign position
		subhalo.position.x = position[0,i];
		subhalo.position.y = position[1,i];
		subhalo.position.z = position[2,i];

		//Assign velocity
		subhalo.velocity.x = velocity[0,i];
		subhalo.velocity.y = velocity[1,i];
		subhalo.velocity.z = velocity[2,i];

		//Assign angular momentum
		subhalo.L[0] = L[0,i];
		subhalo.L[1] = L[1,i];
		subhalo.L[2] = L[2,i];

		subhalo.Vcirc = Vcirc[i];

		//TODO: this needs to be done outside the for, after reading all the subhalos.

		/* Find the corresponding descendant information
		auto it = descendants_data.find(subhalo.id);
		if( it == descendants_data.end() ) {
			ostringstream os;
			os << "No data could be found in the descendants file for subhalo id=" << subhalo.id;
			throw invalid_data(os.str());
		}

		// Fill the missing bits
		auto desc = it->second;
		subhalo.descendant_snapshot = desc.descendant_snapshot;*/

		// Done, save it now
		subhalos.push_back(move(subhalo));
	}



	//First count number of groups as the number of unique hostIndex values.
	std::set<double>unique_values(hostIndex.begin(),hostIndex.end());
	int n_halos = unique_values.size();

	//Assign properties to halos.
	for(unsigned int i=0; i!=n_halos; i++) {
		Halo halo;


	}
	//Calculate virial velocity for halos.


	/**
	 * Loop over subhalos to assign ascendants and descendant arrays, and to define halos.
	 */

	return subhalos;
}


}  // namespace trees

}  // namespace shark


