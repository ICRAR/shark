//
// Galaxy writer classes implementations
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

#include <iomanip>
#include <iostream>
#include <iterator>
#include <memory>
#include <numeric>

#include <boost/filesystem.hpp>

#include "hdf5/writer.h"
#include "components.h"
#include "cosmology.h"
#include "exceptions.h"
#include "galaxy_writer.h"
#include "logging.h"
#include "star_formation.h"
#include "utils.h"


namespace shark {

GalaxyWriter::GalaxyWriter(ExecutionParameters exec_params, CosmologicalParameters cosmo_params, SimulationParameters sim_params, StarFormation starformation):
	exec_params(exec_params),
	sim_params(sim_params),
	cosmo_params(cosmo_params),
	starformation(starformation)
{
	//no-opt
}

std::string GalaxyWriter::get_output_directory(int snapshot)
{
	using namespace boost::filesystem;
	using std::string;

	string batch_dir = "multiple_batches";
	if (exec_params.simulation_batches.size() == 1) {
		batch_dir = std::to_string(exec_params.simulation_batches[0]);
	}

	string output_dir = exec_params.output_directory + "/" + sim_params.sim_name +
	                    "/" + exec_params.name_model + "/" + std::to_string(snapshot) +
	                    "/" + batch_dir;

	// Make sure the directory structure exists
	path dirname(output_dir);
	if (!exists(dirname)) {
		create_directories(dirname);
	}

	return output_dir;
}

void HDF5GalaxyWriter::write(int snapshot, const std::vector<HaloPtr> &halos, TotalBaryon AllBaryons){

	using std::string;
	using std::vector;

	hdf5::Writer file(get_output_directory(snapshot) + "/galaxies.hdf5");

	file.write_dataset("runInfo/batches", exec_params.simulation_batches);
	file.write_dataset("runInfo/ode_solver_precision", exec_params.ode_solver_precision);
	file.write_dataset("runInfo/skip_missing_descendants", exec_params.skip_missing_descendants);
	file.write_dataset("runInfo/snapshot", snapshot);
	file.write_dataset("runInfo/redshift", sim_params.redshifts[snapshot]);
	file.write_attribute("runInfo/model_name", exec_params.name_model);

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

	// Crate all subhalo properties to write.

	vector<long> descendant_id;
	vector<int> main;
	vector<long> id;

	// Create all galaxies properties to write
	vector<float> mstars_disk;
	vector<float> mstars_bulge;
	vector<float> mgas_disk;
	vector<float> mgas_bulge;
	vector<float> mstars_metals_disk;
	vector<float> mstars_metals_bulge;
	vector<float> mgas_metals_disk;
	vector<float> mgas_metals_bulge;
	vector<float> mmol_disk;
	vector<float> mmol_bulge;
	vector<float> matom_disk;
	vector<float> matom_bulge;

	vector<float> mBH;
	vector<float> mBH_acc_hh;
	vector<float> mBH_acc_sb;

	vector<float> sfr_disk;
	vector<float> sfr_burst;

	vector<float> rdisk;
	vector<float> rbulge;
	vector<float> sAM_disk;
	vector<float> sAM_bulge;

	vector<float> mhot;
	vector<float> mhot_metals;

	vector<float> mreheated;
	vector<float> mreheated_metals;

	vector<float> cooling_rate;

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

	vector<int> type;
	vector<Halo::id_t> id_halo;
	vector<Subhalo::id_t> id_subhalo;

	long j=1;
	// Loop over all halos and subhalos to write galaxy properties
	for (auto &halo: halos){

		// assign properties of host halo
		auto mhalo = halo->Mvir;
		auto vhalo = halo->Vvir;
		long i=1;

		for (auto &subhalo: halo->all_subhalos()){

			// assign properties of host subhalo
			auto msubhalo = subhalo->Mvir;
			auto vsubhalo = subhalo->Vcirc;
			auto cnfw = subhalo->concentration;
			auto subhalo_position = subhalo->position;
			auto subhalo_velocity = subhalo->velocity;

			// Assign baryon properties of subhalo
			auto hot_subhalo = subhalo->hot_halo_gas;
			auto cold_subhalo = subhalo->cold_halo_gas;
			auto reheated_subhalo = subhalo->ejected_galaxy_gas;

			descendant_id.push_back(subhalo->descendant_id);
			int m = 0;
			if(subhalo->main_progenitor){
				m = 1;
			}
			main.push_back(m);
			id.push_back(subhalo->id);

			for (auto &galaxy: subhalo->galaxies){

				//Calculate molecular gass mass of disk and bulge:
				double m_mol;
				double m_atom;
				double m_mol_b;
				double m_atom_b;
				starformation.get_molecular_gas(galaxy, sim_params.redshifts[snapshot], &m_mol, &m_atom, &m_mol_b, &m_atom_b);

				mmol_disk.push_back(m_mol);
				mmol_bulge.push_back(m_mol_b);
				matom_disk.push_back(m_atom);
				matom_bulge.push_back(m_atom_b);

				mstars_disk.push_back(galaxy->disk_stars.mass);
				mstars_bulge.push_back(galaxy->bulge_stars.mass);
				mgas_disk.push_back(galaxy->disk_gas.mass);
				mgas_bulge.push_back(galaxy->bulge_gas.mass);
				mstars_metals_disk.push_back(galaxy->disk_stars.mass_metals);
				mstars_metals_bulge.push_back(galaxy->bulge_stars.mass_metals);
				mgas_metals_disk.push_back(galaxy->disk_gas.mass_metals);
				mgas_metals_bulge.push_back(galaxy->bulge_gas.mass_metals);
				sfr_disk.push_back(galaxy->sfr_disk);
				sfr_burst.push_back(galaxy->sfr_bulge);
				mBH.push_back(galaxy->smbh.mass);
				mBH_acc_hh.push_back(galaxy->smbh.macc_hh);
				mBH_acc_sb.push_back(galaxy->smbh.macc_sb);

				rdisk.push_back(galaxy->disk_stars.rscale);
				rbulge.push_back(galaxy->bulge_stars.rscale);

				sAM_disk.push_back(galaxy->disk_stars.sAM);
				sAM_bulge.push_back(galaxy->bulge_stars.sAM);

				double mhot_gal = 0;
				double mzhot_gal = 0;
				double mreheat = 0;
				double mzreheat = 0;
				double rcool = 0;
				int t = galaxy->galaxy_type;
				if(galaxy->galaxy_type == Galaxy::CENTRAL){
					mhot_gal = hot_subhalo.mass + cold_subhalo.mass;
					mzhot_gal = hot_subhalo.mass_metals + cold_subhalo.mass_metals;
					mreheat = reheated_subhalo.mass;
					mzreheat = reheated_subhalo.mass_metals;
					rcool = halo->cooling_rate;
				}

				cooling_rate.push_back(rcool);

				mhot.push_back(mhot_gal);
				mhot_metals.push_back(mzhot_gal);
				mreheated.push_back(mreheat);
				mreheated_metals.push_back(mzreheat);

				mvir_hosthalo.push_back(mhalo);

				double mvir_gal = 0 ;
				double vmax_sub = 0;
				double c_sub = 0;
				xyz<float> pos;
				xyz<float> vel;

				if(galaxy->galaxy_type == Galaxy::CENTRAL || galaxy->galaxy_type == Galaxy::TYPE1){
					mvir_gal = msubhalo;
					vmax_sub = vsubhalo;
					c_sub = cnfw;
					pos = subhalo_position;
					vel = subhalo_velocity;
				}
				else{
					// In case of type 2 galaxies assign negative positions and velocities.
					pos.x = -1;
					pos.y = -1;
					pos.z = -1;
					vel.x = -1;
					vel.y = -1;
					vel.z = -1;
				}

				mvir_subhalo.push_back(mvir_gal);
				vmax_subhalo.push_back(vmax_sub);
				vvir_hosthalo.push_back(vhalo);
				cnfw_subhalo.push_back(c_sub);

				position_x.push_back(pos.x);
				position_y.push_back(pos.y);
				position_z.push_back(pos.z);

				velocity_x.push_back(vel.x);
				velocity_y.push_back(vel.y);
				velocity_z.push_back(vel.z);

				type.push_back(t);

				id_halo.push_back(j);
				id_subhalo.push_back(i);
			}
			i++;
		}
		j++;
	}


	file.write_dataset("Subhalo/id", id);
	file.write_dataset("Subhalo/main_progenitor", main);
	file.write_dataset("Subhalo/descendant_id", descendant_id);

	file.write_dataset("Galaxies/mstars_disk", mstars_disk);
	file.write_dataset("Galaxies/mstars_bulge", mstars_bulge);
	file.write_dataset("Galaxies/mgas_disk", mgas_disk);
	file.write_dataset("Galaxies/mgas_bulge", mgas_bulge);
	file.write_dataset("Galaxies/mstars_metals_disk",mstars_metals_disk);
	file.write_dataset("Galaxies/mstars_metals_bulge", mstars_metals_bulge);
	file.write_dataset("Galaxies/mgas_metals_disk", mgas_metals_disk);
	file.write_dataset("Galaxies/mgas_metals_bulge", mgas_metals_bulge);
	file.write_dataset("Galaxies/mmol_disk",mmol_disk);
	file.write_dataset("Galaxies/mmol_bulge",mmol_bulge);
	file.write_dataset("Galaxies/matom_disk",matom_disk);
	file.write_dataset("Galaxies/matom_bulge",matom_bulge);

	file.write_dataset("Galaxies/sfr_disk", sfr_disk);
	file.write_dataset("Galaxies/sfr_burst", sfr_burst);

	file.write_dataset("Galaxies/mBH", mBH);
	file.write_dataset("Galaxies/BH_accretion_rate_hh", mBH_acc_hh);
	file.write_dataset("Galaxies/BH_accretion_rate_sb", mBH_acc_sb);

	file.write_dataset("Galaxies/rdisk", rdisk);
	file.write_dataset("Galaxies/rbulge", rbulge);

	file.write_dataset("Galaxies/specific_angular_momentum_disk", sAM_disk);
	file.write_dataset("Galaxies/specific_angular_momentum_bulge", sAM_bulge);

	file.write_dataset("Galaxies/mhot", mhot);
	file.write_dataset("Galaxies/mhot_metals", mhot_metals);
	file.write_dataset("Galaxies/mreheated", mreheated);
	file.write_dataset("Galaxies/mreheated_metals", mreheated_metals);
	file.write_dataset("Galaxies/cooling_rate", cooling_rate);

	file.write_dataset("Galaxies/mvir_hosthalo", mvir_hosthalo);
	file.write_dataset("Galaxies/mvir_subhalo", mvir_subhalo);
	file.write_dataset("Galaxies/vmax_subhalo", vmax_subhalo);
	file.write_dataset("Galaxies/vvir_hosthalo", vvir_hosthalo);
	file.write_dataset("Galaxies/cnfw_subhalo", cnfw_subhalo);

	file.write_dataset("Galaxies/position_x", position_x);
	file.write_dataset("Galaxies/position_y", position_y);
	file.write_dataset("Galaxies/position_z", position_z);

	file.write_dataset("Galaxies/velocity_x", velocity_x);
	file.write_dataset("Galaxies/velocity_y", velocity_y);
	file.write_dataset("Galaxies/velocity_z", velocity_z);

	file.write_dataset("Galaxies/type", type);

	file.write_dataset("Galaxies/id_subhalo", id_subhalo);
	file.write_dataset("Galaxies/id_halo", id_halo);

	vector<float> redshifts;

	for (int i=sim_params.min_snapshot; i <= snapshot; i++){
		redshifts.push_back(sim_params.redshifts[i]);
	}

	file.write_dataset("Global/redshifts", redshifts);
	file.write_dataset("Global/mcold",AllBaryons.get_masses(AllBaryons.mcold));
	file.write_dataset("Global/mcold_metals",AllBaryons.get_metals(AllBaryons.mcold));
	file.write_dataset("Global/mstars",AllBaryons.get_masses(AllBaryons.mstars));
	file.write_dataset("Global/mstars_metals",AllBaryons.get_metals(AllBaryons.mstars));
	file.write_dataset("Global/mHI",AllBaryons.get_masses(AllBaryons.mHI));
	file.write_dataset("Global/mH2",AllBaryons.get_masses(AllBaryons.mH2));
	file.write_dataset("Global/mBH",AllBaryons.get_masses(AllBaryons.mBH));
	file.write_dataset("Global/SFR",AllBaryons.SFR);

	file.write_dataset("Global/mhot_halo",AllBaryons.get_masses(AllBaryons.mhot_halo));
	file.write_dataset("Global/mhot_metals",AllBaryons.get_metals(AllBaryons.mhot_halo));
	file.write_dataset("Global/mcold_halo",AllBaryons.get_masses(AllBaryons.mcold_halo));
	file.write_dataset("Global/mcold_halo_metals",AllBaryons.get_metals(AllBaryons.mcold_halo));
	file.write_dataset("Global/mejected_halo",AllBaryons.get_masses(AllBaryons.mejected_halo));
	file.write_dataset("Global/mejected_halo_metals",AllBaryons.get_masses(AllBaryons.mejected_halo));

	file.write_dataset("Global/mDM",AllBaryons.get_masses(AllBaryons.mDM));

}

void ASCIIGalaxyWriter::write(int snapshot, const std::vector<HaloPtr> &halos, TotalBaryon AllBaryons)
{

	using std::vector;
	using std::string;

	std::ofstream output(get_output_directory(snapshot) + "/galaxies.dat");

	// TODO: Write a header?

	// Each galaxy corresponds to one line
	for (const auto &halo: halos) {
		for(const auto &subhalo: halo->all_subhalos()) {
			for(const auto &galaxy: subhalo->galaxies) {
				write_galaxy(galaxy, subhalo, snapshot, output);
			}
		}
	}

	output.close();
}

void ASCIIGalaxyWriter::write_galaxy(const GalaxyPtr &galaxy, const SubhaloPtr &subhalo, int snapshot, std::ofstream &f)
{
	auto mstars_disk = galaxy->disk_stars.mass;
	auto mstars_bulge = galaxy->bulge_stars.mass;
	auto mgas_disk = galaxy->disk_gas.mass;
	auto mgas_metals_disk = galaxy->disk_gas.mass_metals;
	auto mBH = galaxy->smbh.mass;
	auto rdisk = galaxy->disk_stars.rscale;
	auto rbulge = galaxy->bulge_stars.rscale;
	double m_mol;
	double m_atom;
	double m_mol_b;
	double m_atom_b;

	starformation.get_molecular_gas(galaxy, sim_params.redshifts[snapshot], &m_mol, &m_atom, &m_mol_b, &m_atom_b);

	f << mstars_disk << " " << mstars_bulge << " " <<  m_atom + m_atom_b
	  << " " << mBH << " " << mgas_metals_disk / mgas_disk << " "
	  << mstars_disk + mstars_bulge << " " << rdisk << " " << rbulge << " "
	  << subhalo->id << " " << subhalo->host_halo->id << "\n";

}


}// namespace shark

