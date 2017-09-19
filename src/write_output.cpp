/*
 * write_output.cpp

 *
 *  Created on: 18Sep.,2017
 *      Author: clagos
 */

#include <iomanip>
#include <iterator>
#include <memory>
#include <numeric>
#include <vector>

#include "exceptions.h"
#include "logging.h"
#include "write_output.h"
#include "utils.h"
#include "hdf5/writer.h"

using namespace std;

namespace shark {

WriteOutput::WriteOutput(ExecutionParameters exec_params):
	exec_params(exec_params)
	{
		//no-opt
	}

void WriteOutput::write_galaxies(int snapshot, std::vector<HaloPtr> halos){

	string batch;

	if(exec_params.simulation_batches.size() == 1){
		batch = std::to_string(exec_params.simulation_batches[0]);
	}
	else{
		batch = "multiple_batches";
	}

	string fname = exec_params.output_directory + "/" + exec_params.name_model + "/" + std::to_string(snapshot) + batch + "/galaxies.hdf5";

	hdf5::Writer file(fname);

	//exec_params.simulation_batches.write_attribute<unsigned int>("fileInfo/numberOfFiles");
	file.write_dataset_v("runInfo/batches", exec_params.simulation_batches);


	// Loop over all halos and subhalos to write galaxy properties
	for (auto &halo: halos){

		// assign properties of host halo
		auto mhalo = halo->Mvir;
		auto halo_position = halo->position;
		auto halo_velocity = halo->velocity;
		auto vhalo = halo->Vvir;

		for (auto &subhalo: halo->all_subhalos()){

			// assign properties of host subhalo
			auto msubhalo = subhalo->Mvir;
			auto vsubhalo = subhalo->Vvir;
			auto cnfw_subhalo = subhalo->concentration;
			auto L_subhalo = subhalo->L;
			auto subhalo_position = subhalo->position;
			auto subhalo_velocity = subhalo->velocity;

			// Assign baryon properties of subhalo
			auto hot_subhalo = subhalo->hot_halo_gas;
			auto cold_subhalo = subhalo->cold_halo_gas;
			auto reheated_subhalo = subhalo->ejected_galaxy_gas;

			for (auto &galaxy: subhalo->galaxies){



			}
		}
	}


}

}// namespace shark

