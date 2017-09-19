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

#include "components.h"
#include "cosmology.h"
#include "exceptions.h"
#include "logging.h"
#include "write_output.h"
#include "utils.h"
#include "hdf5/writer.h"

using namespace std;

namespace shark {

WriteOutput::WriteOutput(ExecutionParameters exec_params, CosmologicalParameters cosmo_params, SimulationParameters sim_params):
	exec_params(exec_params),
	cosmo_params(cosmo_params),
	sim_params(sim_params)
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
	file.write_dataset("runInfo/model_name", exec_params.name_model);
	file.write_dataset("runInfo/ode_solver_precision", exec_params.ode_solver_precision);
	file.write_dataset("runInfo/skip_missing_descendants", exec_params.skip_missing_descendants);
	file.write_dataset("runInfo/snapshot", snapshot);
	file.write_dataset("runInfo/redshift", sim_params.redshifts[snapshot]);

	// Calculate effective volume of the run
	float volume = sim_params.volume * exec_params.simulation_batches.size();

	file.write_dataset("runInfo/EffectiveVolume", volume);
	file.write_dataset("runInfo/particle_mass", sim_params.particle_mass);

	// Write cosmological parameters
	file.write_dataset("Cosmology/OmegaM", cosmo_params.OmegaM);
	file.write_dataset("Cosmology/OmegaB", cosmo_params.OmegaB);
	file.write_dataset("Cosmology/OmegaL", cosmo_params.OmegaL);
	file.write_dataset("Cosmology/n_s", cosmo_params.n_s);
	file.write_dataset("Cosmology/sigma8", cosmo_params.sigma8);
	file.write_dataset("Cosmology/h", cosmo_params.Hubble_h);

	// Create all galaxies properties I want to write
	vector<float> mstars_disk;
	vector<float> mstars_bulge;
	vector<float> mgas_disk;
	vector<float> mgas_bulge;
	vector<float> mstars_metals_disk;
	vector<float> mstars_metals_bulge;
	vector<float> mgas_metals_disk;
	vector<float> mgas_metals_bulge;
	vector<float> mBH;

	vector<float> rdisk;
	vector<float> rbulge;

	vector<float> mhot;
	vector<float> mhot_metals;

	vector<float> mreheated;
	vector<float> mreheated_metals;

	vector<float> mvir_hosthalo;
	vector<float> mvir_subhalo;
	vector<float> vmax_subhalo;
	vector<float> vvir_hosthalo;

	vector<float> cnfw_subhalo;

	vector<float> position_x;
	vector<float> position_y;
	vector<float> position_z;

	vector<float> velocity_x;
	vector<float> velocity_y;
	vector<float> velocity_z;

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

void WriteOutput::write_cosmology(){




}

}// namespace shark

