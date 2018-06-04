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
#include "timer.h"
#include "utils.h"


namespace shark {

GalaxyWriter::GalaxyWriter(ExecutionParameters exec_params, CosmologicalParameters cosmo_params,  std::shared_ptr<Cosmology> cosmology, std::shared_ptr<DarkMatterHalos> darkmatterhalo, SimulationParameters sim_params):
	exec_params(exec_params),
	cosmo_params(cosmo_params),
	cosmology(cosmology),
	darkmatterhalo(darkmatterhalo),
	sim_params(sim_params){
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

void HDF5GalaxyWriter::write(int snapshot, const std::vector<HaloPtr> &halos, TotalBaryon &AllBaryons, const molgas_per_galaxy &molgas_per_gal){

	using std::string;
	using std::vector;

	//Write output with the number of the coming snapshot. This is because we evolved galaxies to the end of the current snapshot.
	int snap_to_write = snapshot + 1;

	string comment;

	hdf5::Writer file(get_output_directory(snap_to_write) + "/galaxies.hdf5");

	//Write header
	write_header(file, snap_to_write);

	//Write galaxies
	write_galaxies(file, snap_to_write, halos, molgas_per_gal);

	//Write total baryon components
	write_global_properties(file, snap_to_write, AllBaryons);

	// Write star formation histories.
	write_histories(snap_to_write, halos);

}

void HDF5GalaxyWriter::write_header(hdf5::Writer &file, int snapshot){

	std::string comment;

	comment = "number of batches analysed";
	file.write_dataset("runInfo/batches", exec_params.simulation_batches, comment);

	comment = "accuracy applied when solving the ODE system of the physical model.";
	file.write_dataset("runInfo/ode_solver_precision", exec_params.ode_solver_precision, comment);

	comment = "boolean parameter that sets whether the code ignores subhalos that have no descendants.";
	file.write_dataset("runInfo/skip_missing_descendants", exec_params.skip_missing_descendants, comment);

	comment = "output snapshot";
	file.write_dataset("runInfo/snapshot", snapshot, comment);

	comment = "output redshift";
	file.write_dataset("runInfo/redshift", sim_params.redshifts[snapshot], comment);

	file.write_attribute("runInfo/model_name", exec_params.name_model);

	// Calculate effective volume of the run
	float volume = sim_params.volume * exec_params.simulation_batches.size();

	comment = "effective volume of this run [cMpc/h]";
	file.write_dataset("runInfo/EffectiveVolume", volume, comment);

	comment = "dark matter particle mass of this simulation [Msun/h]";
	file.write_dataset("runInfo/particle_mass", sim_params.particle_mass, comment);

	comment = "Box side size of the full simulated volume [Mpc/h]";
	file.write_dataset("runInfo/lbox", sim_params.lbox, comment);

	comment = "Total number of subvolumes in which the simulated box was divided into";
	file.write_dataset("runInfo/tot_n_subvolumes", sim_params.tot_nsubvols, comment);

	// Write cosmological parameters

	comment = "omega matter assumed in simulation";
	file.write_dataset("Cosmology/OmegaM", cosmo_params.OmegaM, comment);

	comment = "omega baryon assumed in simulation";
	file.write_dataset("Cosmology/OmegaB", cosmo_params.OmegaB, comment);

	comment = "omega lambda assumed in simulation";
	file.write_dataset("Cosmology/OmegaL", cosmo_params.OmegaL, comment);

	comment = "scalar spectral index assumed in simulation";
	file.write_dataset("Cosmology/n_s", cosmo_params.n_s, comment);

	comment = "fluctuation amplitude at 8 Mpc/h";
	file.write_dataset("Cosmology/sigma8", cosmo_params.sigma8, comment);

	comment = "normalization of hubble parameter H0 = h * 100 Mpc * km/s";
	file.write_dataset("Cosmology/h", cosmo_params.Hubble_h, comment);
}

template<typename T>
static inline
std::size_t report_vsize(const std::vector<T> &v, std::ostringstream &os, const char *name) {
	const std::size_t amount = sizeof(T) * v.size();
	os << " " << name << ": " << memory_amount(amount);
	return amount;
};

void HDF5GalaxyWriter::write_galaxies(hdf5::Writer &file, int snapshot, const std::vector<HaloPtr> &halos, const molgas_per_galaxy &molgas_per_gal){

	Timer t;


	using std::string;
	using std::vector;

	string comment;

	// Crate all subhalo properties to write.

	vector<long> descendant_id;
	vector<int> main;
	vector<long> id;
	vector<long> host_id;
	vector<long> id_galaxy;

	// Create all galaxies properties to write
	vector<float> mstars_disk;
	vector<float> mstars_bulge;
	vector<float> mstars_burst;
	vector<float> mgas_disk;
	vector<float> mgas_bulge;
	vector<float> mstars_metals_disk;
	vector<float> mstars_metals_bulge;
	vector<float> mstars_metals_burst;
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

	vector<float> rdisk_gas;
	vector<float> rbulge_gas;
	vector<float> sAM_disk_gas;
	vector<float> sAM_disk_gas_atom;
	vector<float> sAM_disk_gas_mol;
	vector<float> sAM_bulge_gas;

	vector<float> rdisk_star;
	vector<float> rbulge_star;
	vector<float> sAM_disk_star;
	vector<float> sAM_bulge_star;

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
	vector<float> lambda_subhalo;

	vector<float> position_x;
	vector<float> position_y;
	vector<float> position_z;

	vector<float> velocity_x;
	vector<float> velocity_y;
	vector<float> velocity_z;

	vector<float> L_x;
	vector<float> L_y;
	vector<float> L_z;

	vector<int> type;
	vector<Halo::id_t> id_halo;
	vector<Halo::id_t> id_halo_tree;
	vector<Subhalo::id_t> id_subhalo;
	vector<Subhalo::id_t> id_subhalo_tree;



	long j = 1;
	long gal_id = 1;
	// Loop over all halos and subhalos to write galaxy properties
	for (auto &halo: halos){

		// assign properties of host halo
		auto mhalo = halo->Mvir;
		auto vhalo = halo->Vvir;
		long i=1;

		for (auto &subhalo: halo->all_subhalos()){

			host_id.push_back(halo->id);

			// assign properties of host subhalo
			auto msubhalo = subhalo->Mvir;
			auto vsubhalo = subhalo->Vcirc;
			auto cnfw = subhalo->concentration;
			auto lambda = subhalo->lambda;
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

			for (const auto &galaxy: subhalo->galaxies){

				id_halo_tree.push_back(halo->id);
				id_subhalo_tree.push_back(subhalo->id);

				//Calculate molecular gas mass of disk and bulge, and specific angular momentum in atomic/molecular disk.
				auto &molecular_gas = molgas_per_gal.at(galaxy);

				// Gas components separated into HI and H2.
				mmol_disk.push_back(molecular_gas.m_mol);
				mmol_bulge.push_back(molecular_gas.m_mol_b);
				matom_disk.push_back(molecular_gas.m_atom);
				matom_bulge.push_back(molecular_gas.m_atom_b);

				// Stellar components
				mstars_disk.push_back(galaxy->disk_stars.mass);
				mstars_bulge.push_back(galaxy->bulge_stars.mass);
				mstars_burst.push_back(galaxy->burst_stars.mass);

				// Gas components
				mgas_disk.push_back(galaxy->disk_gas.mass);
				mgas_bulge.push_back(galaxy->bulge_gas.mass);

				// Metals of the stellar components.
				mstars_metals_disk.push_back(galaxy->disk_stars.mass_metals);
				mstars_metals_bulge.push_back(galaxy->bulge_stars.mass_metals);
				mstars_metals_burst.push_back(galaxy->burst_stars.mass_metals);

				// Metals of the gas components.
				mgas_metals_disk.push_back(galaxy->disk_gas.mass_metals);
				mgas_metals_bulge.push_back(galaxy->bulge_gas.mass_metals);

				// SFRs in disks and bulges.
				sfr_disk.push_back(galaxy->sfr_disk);
				sfr_burst.push_back(galaxy->sfr_bulge);

				// Black hole properties.
				mBH.push_back(galaxy->smbh.mass);
				mBH_acc_hh.push_back(galaxy->smbh.macc_hh);
				mBH_acc_sb.push_back(galaxy->smbh.macc_sb);

				// Sizes and specific angular momentum of disks and bulges.

				rdisk_gas.push_back(galaxy->disk_gas.rscale);
				rbulge_gas.push_back(galaxy->bulge_gas.rscale);
				sAM_disk_gas.push_back(galaxy->disk_gas.sAM);
				sAM_disk_gas_atom.push_back(molecular_gas.j_atom);
				sAM_disk_gas_mol.push_back(molecular_gas.j_mol);
				sAM_bulge_gas.push_back(galaxy->bulge_gas.sAM);

				rdisk_star.push_back(galaxy->disk_stars.rscale);
				rbulge_star.push_back(galaxy->bulge_stars.rscale);
				sAM_disk_star.push_back(galaxy->disk_stars.sAM);
				sAM_bulge_star.push_back(galaxy->bulge_stars.sAM);

				// Halo properties below.
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
				double l_sub = 0;
				xyz<float> pos;
				xyz<float> vel;
				xyz<float> L;

				if(galaxy->galaxy_type == Galaxy::CENTRAL || galaxy->galaxy_type == Galaxy::TYPE1){
					mvir_gal = msubhalo;
					vmax_sub = vsubhalo;
					c_sub    = cnfw;
					l_sub    = lambda;
					pos      = subhalo_position;
					vel      = subhalo_velocity;
					L.x      = (subhalo->L.x/subhalo->L.norm()) * galaxy->angular_momentum();
					L.y      = (subhalo->L.y/subhalo->L.norm()) * galaxy->angular_momentum();
					L.z      = (subhalo->L.z/subhalo->L.norm()) * galaxy->angular_momentum();
				}
				else{
					// In case of type 2 galaxies assign negative positions, velocities and angular momentum.
					darkmatterhalo->generate_random_orbits(pos, vel, L, galaxy->angular_momentum(), halo);
				}

				mvir_subhalo.push_back(mvir_gal);
				vmax_subhalo.push_back(vmax_sub);
				vvir_hosthalo.push_back(vhalo);
				cnfw_subhalo.push_back(c_sub);
				lambda_subhalo.push_back(l_sub);

				// Galaxy position and velocity.
				position_x.push_back(pos.x);
				position_y.push_back(pos.y);
				position_z.push_back(pos.z);

				velocity_x.push_back(vel.x);
				velocity_y.push_back(vel.y);
				velocity_z.push_back(vel.z);

				L_x.push_back(cosmology->comoving_to_physical_angularmomentum(L.x,sim_params.redshifts[snapshot]));
				L_y.push_back(cosmology->comoving_to_physical_angularmomentum(L.y,sim_params.redshifts[snapshot]));
				L_z.push_back(cosmology->comoving_to_physical_angularmomentum(L.z,sim_params.redshifts[snapshot]));

				type.push_back(t);

				id_halo.push_back(j);
				id_subhalo.push_back(i);

				id_galaxy.push_back(gal_id);

				gal_id ++;
			}
			i++;
		}
		j++;
	}

	std::ostringstream os;
	std::size_t total = 0;
#define REPORT(x) total += report_vsize(x, os, #x)
	REPORT(descendant_id);
	REPORT(main);
	REPORT(id);
	REPORT(id_galaxy);
	REPORT(host_id);
	REPORT(mstars_disk);
	REPORT(mstars_bulge);
	REPORT(mstars_burst);
	REPORT(mgas_disk);
	REPORT(mgas_bulge);
	REPORT(mstars_metals_disk);
	REPORT(mstars_metals_bulge);
	REPORT(mstars_metals_burst);
	REPORT(mgas_metals_disk);
	REPORT(mgas_metals_bulge);
	REPORT(mmol_disk);
	REPORT(mmol_bulge);
	REPORT(matom_disk);
	REPORT(matom_bulge);
	REPORT(mBH);
	REPORT(mBH_acc_hh);
	REPORT(mBH_acc_sb);
	REPORT(sfr_disk);
	REPORT(sfr_burst);
	REPORT(rdisk_gas);
	REPORT(rbulge_gas);
	REPORT(sAM_disk_gas);
	REPORT(sAM_disk_gas_atom);
	REPORT(sAM_disk_gas_mol);
	REPORT(sAM_bulge_gas);
	REPORT(rdisk_star);
	REPORT(rbulge_star);
	REPORT(sAM_disk_star);
	REPORT(sAM_bulge_star);
	REPORT(mhot);
	REPORT(mhot_metals);
	REPORT(mreheated);
	REPORT(mreheated_metals);
	REPORT(cooling_rate);
	REPORT(mvir_hosthalo);
	REPORT(mvir_subhalo);
	REPORT(vmax_subhalo);
	REPORT(vvir_hosthalo);
	REPORT(cnfw_subhalo);
	REPORT(lambda_subhalo);
	REPORT(position_x);
	REPORT(position_y);
	REPORT(position_z);
	REPORT(velocity_x);
	REPORT(velocity_y);
	REPORT(velocity_z);
	REPORT(L_x);
	REPORT(L_y);
	REPORT(L_z);
	REPORT(type);
	REPORT(id_halo);
	REPORT(id_subhalo);

	LOG(info) << "Total amount of memory used by the writing process: " << memory_amount(total);
	LOG(debug) << "Detailed amounts follow: " << os.str();

	LOG(info) << "Galaxies pivoted and memory reported in " << t;

	t = Timer();

	//Write subhalo properties.
	comment = "Subhalo id";
	file.write_dataset("Subhalo/id", id, comment);

	comment = "=1 if subhalo is the main progenitor' =0 otherwise.";
	file.write_dataset("Subhalo/main_progenitor", main, comment);

	comment = "id of the subhalo that is the descendant of this subhalo";
	file.write_dataset("Subhalo/descendant_id", descendant_id, comment);

	comment = "id of the host halo of this subhalo";
	file.write_dataset("Subhalo/host_id", host_id, comment);

	//Write galaxy properties.
	comment = "stellar mass in the disk [Msun/h]";
	file.write_dataset("Galaxies/mstars_disk", mstars_disk, comment);

	comment = "stellar mass in the bulge [Msun/h]";
	file.write_dataset("Galaxies/mstars_bulge", mstars_bulge, comment);

	comment = "stellar mass formed via starbursts [Msun/h]";
	file.write_dataset("Galaxies/mstars_burst", mstars_burst, comment);

	comment = "total gas mass in the disk [Msun/h]";
	file.write_dataset("Galaxies/mgas_disk", mgas_disk, comment);

	comment = "gas mass in the bulge [Msun/h]";
	file.write_dataset("Galaxies/mgas_bulge", mgas_bulge, comment);

	comment = "mass of metals locked in stars in the disk [Msun/h]";
	file.write_dataset("Galaxies/mstars_metals_disk",mstars_metals_disk, comment);

	comment = "mass of metals locked in stars in the bulge [Msun/h]";
	file.write_dataset("Galaxies/mstars_metals_bulge", mstars_metals_bulge, comment);

	comment = "mass of metals locked in stars that formed via starbursts [Msun/h]";
	file.write_dataset("Galaxies/mstars_metals_burst", mstars_metals_burst, comment);

	comment = "mass of metals locked in the gas of the disk [Msun/h]";
	file.write_dataset("Galaxies/mgas_metals_disk", mgas_metals_disk, comment);

	comment = "mass of metals locked in the gas of the bulge [Msun/h]";
	file.write_dataset("Galaxies/mgas_metals_bulge", mgas_metals_bulge, comment);

	comment = "molecular gas mass (helium plus hydrogen) in the disk [Msun/h]";
	file.write_dataset("Galaxies/mmol_disk",mmol_disk, comment);

	comment ="molecular gas mass (helium plus hydrogen) in the bulge [Msun/h]";
	file.write_dataset("Galaxies/mmol_bulge",mmol_bulge, comment);

	comment = "atomic gas mass (helium plus hydrogen) in the disk [Msun/h]";
	file.write_dataset("Galaxies/matom_disk",matom_disk, comment);

	comment ="atomic gas mass (helium plus hydrogen) in the bulge [Msun/h]";
	file.write_dataset("Galaxies/matom_bulge",matom_bulge, comment);

	comment = "star formation rate in the disk [Msun/Gyr/h]";
	file.write_dataset("Galaxies/sfr_disk", sfr_disk, comment);

	comment = "star formation rate in the bulge [Msun/Gyr/h]";
	file.write_dataset("Galaxies/sfr_burst", sfr_burst, comment);

	comment = "black hole mass [Msun/h]";
	file.write_dataset("Galaxies/mBH", mBH, comment);

	comment = "accretion rate onto the black hole during the hot halo mode [Msun/Gyr/h]";
	file.write_dataset("Galaxies/BH_accretion_rate_hh", mBH_acc_hh, comment);

	comment = "accretion rate onto the black hole during the starburst mode [Msun/Gyr/h]";
	file.write_dataset("Galaxies/BH_accretion_rate_sb", mBH_acc_sb, comment);

	comment = "half-mass radius of the stellar disk [cMpc/h]";
	file.write_dataset("Galaxies/rdisk_star", rdisk_star, comment);

	comment = "half-mass radius of the stellar bulge [cMpc/h]";
	file.write_dataset("Galaxies/rbulge_star", rbulge_star, comment);

	comment = "specific angular momentum of the stellar disk [km/s * cMpc/h]";
	file.write_dataset("Galaxies/specific_angular_momentum_disk_star", sAM_disk_star, comment);

	comment = "specific angular momentum of the stellar bulge [km/s * cMpc/h]";
	file.write_dataset("Galaxies/specific_angular_momentum_bulge_star", sAM_bulge_star, comment);

	comment = "half-mass radius of the gas disk [cMpc/h]";
	file.write_dataset("Galaxies/rdisk_gas", rdisk_gas, comment);

	comment = "half-mass radius of the gas bulge [cMpc/h]";
	file.write_dataset("Galaxies/rbulge_gas", rbulge_gas, comment);

	comment = "specific angular momentum of the gas disk [km/s * cMpc/h]";
	file.write_dataset("Galaxies/specific_angular_momentum_disk_gas", sAM_disk_gas, comment);

	comment = "specific angular momentum of the atomic gas disk [km/s * cMpc/h]";
	file.write_dataset("Galaxies/specific_angular_momentum_disk_gas_atom", sAM_disk_gas_atom, comment);

	comment = "specific angular momentum of the molecular gas disk [km/s * cMpc/h]";
	file.write_dataset("Galaxies/specific_angular_momentum_disk_gas_mol", sAM_disk_gas_mol, comment);

	comment = "specific angular momentum of the gas bulge [km/s * cMpc/h]";
	file.write_dataset("Galaxies/specific_angular_momentum_bulge_gas", sAM_bulge_gas, comment);

	comment = "hot gas mass in the halo [Msun/h]";
	file.write_dataset("Galaxies/mhot", mhot, comment);

	comment = "mass of metals locked in the hot halo gas [Msun/h]";
	file.write_dataset("Galaxies/mhot_metals", mhot_metals, comment);

	comment = "gas mass in the ejected gas component [Msun/h]";
	file.write_dataset("Galaxies/mreheated", mreheated, comment);

	comment = "mass of metals locked in the ejected gas component [Msun/h]";
	file.write_dataset("Galaxies/mreheated_metals", mreheated_metals, comment);

	comment = "cooling rate of the hot halo component [Msun/Gyr/h].";
	file.write_dataset("Galaxies/cooling_rate", cooling_rate, comment);

	comment = "Dark matter mass of the host halo in which this galaxy resides [Msun/h]";
	file.write_dataset("Galaxies/mvir_hosthalo", mvir_hosthalo, comment);

	comment = "Dark matter mass of the subhalo in which this galaxy resides [Msun/h]";
	file.write_dataset("Galaxies/mvir_subhalo", mvir_subhalo, comment);

	comment = "Maximum circular velocity of the dark matter subhalo in which this galaxy resides [km/s]";
	file.write_dataset("Galaxies/vmax_subhalo", vmax_subhalo, comment);

	comment = "Virial velocity of the dark matter halo in which this galaxy resides [km/s]";
	file.write_dataset("Galaxies/vvir_hosthalo", vvir_hosthalo, comment);

	comment = "NFW concentration parameter of the dark matter subhalo in which this galaxy resides [dimensionless]";
	file.write_dataset("Galaxies/cnfw_subhalo", cnfw_subhalo, comment);

	comment = "Spin parameter of the dark matter subhalo in which this galaxy resides [dimensionless]";
	file.write_dataset("Galaxies/lambda_subhalo", lambda_subhalo, comment);

	//Galaxy position
	comment = "position component x of galaxy [cMpc/h]. In the case of type 2 galaxies, the positions are generated to randomly sample an NFW halo with the concentration of the halo the galaxy lives in.";
	file.write_dataset("Galaxies/position_x", position_x, comment);
	comment = "position component y of galaxy [cMpc/h]. In the case of type 2 galaxies, the positions are generated to randomly sample an NFW halo with the concentration of the halo the galaxy lives in.";
	file.write_dataset("Galaxies/position_y", position_y, comment);
	comment = "position component z of galaxy [cMpc/h]. In the case of type 2 galaxies, the positions are generated to randomly sample an NFW halo with the concentration of the halo the galaxy lives in.";
	file.write_dataset("Galaxies/position_z", position_z, comment);

	//Galaxy velocity
	comment = "peculiar velocity component x of galaxy [km/s]. In the case of type 2 galaxies, the velocity is generated to randomly sample the velocity dispersion of a NFW halo with the concentration of the halo the galaxy lives in.";
	file.write_dataset("Galaxies/velocity_x", velocity_x, comment);
	comment = "peculiar velocity component y of galaxy [km/s]. In the case of type 2 galaxies, the velocity is generated to randomly sample the velocity dispersion of a NFW halo with the concentration of the halo the galaxy lives in.";
	file.write_dataset("Galaxies/velocity_y", velocity_y, comment);
	comment = "peculiar velocity component z of galaxy [km/s]. In the case of type 2 galaxies, the velocity is generated to randomly sample the velocity dispersion of a NFW halo with the concentration of the halo the galaxy lives in.";
	file.write_dataset("Galaxies/velocity_z", velocity_z, comment);

	//Galaxy AM vector
	comment = "total angular momentum component x of galaxy [Msun pMpc km/s]. In the case of type 2 galaxies, the AM vector is randomly oriented.";
	file.write_dataset("Galaxies/L_x", L_x,  comment);
	comment = "total angular momentum component y of galaxy [Msun pMpc km/s]. In the case of type 2 galaxies, the AM vector is randomly oriented.";
	file.write_dataset("Galaxies/L_y", L_y, comment);
	comment = "total angular momentum component z of galaxy [Msun pMpc km/s]. In the case of type 2 galaxies, the AM vector is randomly oriented.";
	file.write_dataset("Galaxies/L_z", L_z, comment);

	//Galaxy type.
	comment = "galaxy type; =0 for centrals; =1 for satellites that reside in well identified subhalos; =2 for orphan satellites";
	file.write_dataset("Galaxies/type", type, comment);

	//Galaxy IDs.
	comment = "subhalo ID. Unique to this snapshot.";
	file.write_dataset("Galaxies/id_subhalo", id_subhalo, comment);
	comment = "halo ID. Unique to this snapshot.";
	file.write_dataset("Galaxies/id_halo", id_halo, comment);
	comment = "galaxy ID. Unique to this snapshot.";
	file.write_dataset("Galaxies/id_galaxy", id_galaxy, comment);
	comment = "subhalo id in the tree (unique to entire halo catalogue).";
	file.write_dataset("Galaxies/id_subhalo_tree", id_subhalo_tree, comment);
	comment = "halo id in the tree (unique to entire halo catalogue).";
	file.write_dataset("Galaxies/id_halo_tree", id_halo_tree, comment);

	LOG(info) << "Galaxies data written in " << t;

}

void HDF5GalaxyWriter::write_global_properties (hdf5::Writer &file, int snapshot, TotalBaryon &AllBaryons){

	using std::string;
	using std::vector;

	string comment;

	vector<float> redshifts;
	vector<double> baryons_ever_created;
	vector<double> baryons_ever_lost;

	double baryons_lost = 0;

	for (int i=sim_params.min_snapshot+1; i <= snapshot; i++){
		redshifts.push_back(sim_params.redshifts[i]);
		baryons_ever_created.push_back(AllBaryons.baryon_total_created[i]);

		// Accummulate baryons lost.
		baryons_lost += AllBaryons.baryon_total_lost[i];
		baryons_ever_lost.push_back(baryons_lost);
	}

	comment = "redshifts of the global outputs.";
	file.write_dataset("Global/redshifts", redshifts, comment);

	comment = "total cold gas mass (interstellar medium) in the simulated box [Msun/h]";
	file.write_dataset("Global/mcold",AllBaryons.get_masses(AllBaryons.mcold), comment);

	comment = "total mass of metals locked in cold gas in the simulated box [Msun/h]";
	file.write_dataset("Global/mcold_metals",AllBaryons.get_metals(AllBaryons.mcold), comment);

	comment = "total stellar mass in the simulated box [Msun/h]";
	file.write_dataset("Global/mstars",AllBaryons.get_masses(AllBaryons.mstars), comment);

	comment = "total mass of metals locked in stars in the simulated box [Msun/h]";
	file.write_dataset("Global/mstars_metals",AllBaryons.get_metals(AllBaryons.mstars), comment);

	comment = "total stellar mass formed via starbursts in the simulated box [Msun/h]";
	file.write_dataset("Global/mstars_bursts",AllBaryons.get_masses(AllBaryons.mstars_burst), comment);

	comment = "total mass of metals locked in stars that formed via starbursts in the simulated box [Msun/h]";
	file.write_dataset("Global/mstars_metals_bursts",AllBaryons.get_metals(AllBaryons.mstars_burst), comment);

	comment = "total atomic gas mass in the simulated box [Msun/h]";
	file.write_dataset("Global/mHI",AllBaryons.get_masses(AllBaryons.mHI), comment);

	comment = "total molecular gas mass in the simulated box [Msun/h]";
	file.write_dataset("Global/mH2",AllBaryons.get_masses(AllBaryons.mH2), comment);

	comment = "total mass locked up in black holes in the simulated box [Msun/h]";
	file.write_dataset("Global/mBH",AllBaryons.get_masses(AllBaryons.mBH), comment);

	comment = "total star formation rate taking place in disks in the simulated box [Msun/Gyr/h]";
	file.write_dataset("Global/SFR_quiescent",AllBaryons.SFR_disk, comment);

	comment = "total star formation rate taking place in bulges in the simulated box [Msun/Gyr/h]";
	file.write_dataset("Global/SFR_burst",AllBaryons.SFR_bulge, comment);

	comment = "total hot gas mass in halos in the simulated box [Msun/h]";
	file.write_dataset("Global/mhot_halo",AllBaryons.get_masses(AllBaryons.mhot_halo),comment);
	comment = "total mass of metals in the hot gas mass in halos in the simulated box [Msun/h]";
	file.write_dataset("Global/mhot_metals",AllBaryons.get_metals(AllBaryons.mhot_halo), comment);

	comment = "total halo cold gas in the simulated box [Msun/h]";
	file.write_dataset("Global/mcold_halo",AllBaryons.get_masses(AllBaryons.mcold_halo), comment);
	comment = "total mass of metals in the halo cold gas mass in the simulated box [Msun/h]";
	file.write_dataset("Global/mcold_halo_metals",AllBaryons.get_metals(AllBaryons.mcold_halo), comment);

	comment = "total gas mass ejected from halos (and that has not yet been reincorporated) in the simulated box [Msun/h]";
	file.write_dataset("Global/mejected_halo",AllBaryons.get_masses(AllBaryons.mejected_halo), comment);
	comment = "total mass of metals in the ejected gas reservoir in the simulated box [Msun/h]";
	file.write_dataset("Global/mejected_halo_metals",AllBaryons.get_metals(AllBaryons.mejected_halo), comment);

	comment = "total dark matter mass locked up in halos in the simulated box [Msun/h].";
	file.write_dataset("Global/mDM",AllBaryons.get_masses(AllBaryons.mDM), comment);

	comment = "total baryon mass in the simulated box [Msun/h]";
	file.write_dataset("Global/mbar_created",baryons_ever_created, comment);
	comment = "total baryons lost in the simulated box [Msun/h] (ideally this should be =0)";
	file.write_dataset("Global/mbar_lost", baryons_ever_lost, comment);
}

void HDF5GalaxyWriter::write_histories (int snapshot, const std::vector<HaloPtr> &halos){


	using std::string;
	using std::vector;

	string comment;

	if(exec_params.output_sf_histories){
		if(std::find(exec_params.snapshots_sf_histories.begin(), exec_params.snapshots_sf_histories.end(), snapshot) != exec_params.snapshots_sf_histories.end()){
			hdf5::Writer file_sfh(get_output_directory(snapshot) + "/star_formation_histories.hdf5");

			//Create the vectors that will save the information of the galaxies
			vector<vector<float>> sfhs_disk;
			vector<vector<float>> stellar_mass_disk;
			vector<vector<float>> stellar_metals_disk;
			vector<vector<float>> gas_hs_disk;
			vector<vector<float>> gas_metals_hs_disk;

			vector<vector<float>> sfhs_bulge;
			vector<vector<float>> stellar_mass_bulge;
			vector<vector<float>> stellar_metals_bulge;
			vector<vector<float>> gas_hs_bulge;
			vector<vector<float>> gas_metals_hs_bulge;

			vector<long> id_galaxy;
			long gal_id = 1;

			float defl_value = -1;

			for (auto &halo: halos){
				for (auto &subhalo: halo->all_subhalos()){
					for (auto &galaxy: subhalo->galaxies){

						vector<float> sfh_gal_disk;
						vector<float> star_metals_gal_disk;
						vector<float> sfh_gal_bulge;
						vector<float> star_metals_gal_bulge;

						bool star_gal_bulge_exists = false;
						for(int s=sim_params.min_snapshot+1; s <= snapshot; s++) {

							auto it = std::find_if(galaxy->history.begin(), galaxy->history.end(), [s](const HistoryItem &hitem) {
								//information in snapshot corresponds to the end of it, so effectively, when writing, we need to
								//compare to s-1.
								return hitem.snapshot == s-1;
							});

							if (it == galaxy->history.end()) {
								if (star_gal_bulge_exists) {
									std::ostringstream os;
									os << "The history of the StellarMass of the bulge of " << galaxy << " ceased to exist (temporarily). ";
									os << "These are the snapshots for which there is a history item: ";
									vector<int> hsnaps(galaxy->history.size());
									std::transform(galaxy->history.begin(), galaxy->history.end(), hsnaps.begin(), [](const HistoryItem &hitem) {
										return hitem.snapshot;
									});
									std::copy(hsnaps.begin(), hsnaps.end(), std::ostream_iterator<int>(os, " "));
									LOG(warning) << os.str();
									for (auto &item: galaxy->history){
										std::ostringstream os;
										os << "snap history: "<< item.snapshot;
										LOG(warning) << os.str();
									}
								}
								sfh_gal_disk.push_back(defl_value);
								star_metals_gal_disk.push_back(defl_value);

								sfh_gal_bulge.push_back(defl_value);
								star_metals_gal_bulge.push_back(defl_value);
							}
							else {
								star_gal_bulge_exists = true;
								auto item = *it;
								// assign disk properties
								sfh_gal_disk.push_back(item.sfr_disk/constants::GIGA);
								if(item.sfr_disk > 0){
									star_metals_gal_disk.push_back(item.sfr_z_disk/item.sfr_disk);
								}
								else{
									star_metals_gal_disk.push_back(0);
								}

								// assign bulge properties
								sfh_gal_bulge.push_back(item.sfr_bulge/constants::GIGA);
								if(item.sfr_bulge > 0){
									star_metals_gal_bulge.push_back(item.sfr_z_bulge/item.sfr_bulge);
								}
								else{
									star_metals_gal_bulge.push_back(0);
								}
							}
						}

						// save galaxies only if they have a stellar mass >0 by the output snapshot.
						if(galaxy->stellar_mass() > 0){
							sfhs_disk.emplace_back(std::move(sfh_gal_disk));
							stellar_metals_disk.emplace_back(std::move(star_metals_gal_disk));

							sfhs_bulge.emplace_back(std::move(sfh_gal_bulge));
							stellar_metals_bulge.emplace_back(std::move(star_metals_gal_bulge));

							id_galaxy.push_back(gal_id);
						}

						gal_id ++;
					}
				}
			}

			vector<float> redshifts;
			vector<float> age_mean;
			vector<float> delta_t;

			double age_uni = std::abs(cosmology->convert_redshift_to_age(0));
			for (int i=sim_params.min_snapshot+1; i <= snapshot; i++){
				redshifts.push_back(sim_params.redshifts[i]);
				double delta = std::abs(cosmology->convert_redshift_to_age(sim_params.redshifts[i]) - cosmology->convert_redshift_to_age(sim_params.redshifts[i-1]));
				double age = age_uni - 0.5 * (std::abs(cosmology->convert_redshift_to_age(sim_params.redshifts[i]) + cosmology->convert_redshift_to_age(sim_params.redshifts[i-1])));
				delta_t.push_back(delta);
				age_mean.push_back(age);
			}

			//Write header
			write_header(file_sfh, snapshot);

			comment = "galaxy ID. Unique to this snapshot.";
			file_sfh.write_dataset("Galaxies/id_galaxy", id_galaxy, comment);

			//Write disk component history.
			comment = "Star formation history of stars formed that by this output time end up in the disk [Msun/yr/h]";
			file_sfh.write_dataset("Disks/StarFormationRateHistories", sfhs_disk, comment);

			comment = "Stellar metallicity of the stars formed in a timestep that by this output time ends up in the disk";
			file_sfh.write_dataset("Disks/MetallicityHistories", stellar_metals_disk, comment);

			//Write bulge component history.
			comment = "Star formation history of stars formed that by this output time end up in the bulge [Msun/yr/h]";
			file_sfh.write_dataset("Bulges/StarFormationRateHistories", sfhs_bulge, comment);

			comment = "Stellar metallicity of the stars formed in a timestep that by this output time ends up in the disk";
			file_sfh.write_dataset("Bulges/MetallicityHistories", stellar_metals_bulge, comment);

			comment = "Redshifts of the history outputs";
			file_sfh.write_dataset("Redshifts", redshifts, comment);

			comment = "Look back time to mean time between snapshots [Gyr]";
			file_sfh.write_dataset("age_mean", age_mean, comment);

			comment = "Time interval covered between snapshots [Gyr]";
			file_sfh.write_dataset("delta_t", delta_t, comment);

		}

	}
}

void ASCIIGalaxyWriter::write(int snapshot, const std::vector<HaloPtr> &halos, TotalBaryon &AllBaryons, const molgas_per_galaxy &molgas_per_gal)
{

	using std::vector;
	using std::string;

	std::ofstream output(get_output_directory(snapshot) + "/galaxies.dat");

	// TODO: Write a header?

	// Each galaxy corresponds to one line
	for (const auto &halo: halos) {
		for(const auto &subhalo: halo->all_subhalos()) {
			for(const auto &galaxy: subhalo->galaxies) {
				write_galaxy(galaxy, subhalo, snapshot, output, molgas_per_gal);
			}
		}
	}

	output.close();
}

void ASCIIGalaxyWriter::write_galaxy(const GalaxyPtr &galaxy, const SubhaloPtr &subhalo, int snapshot, std::ofstream &f, const molgas_per_galaxy &molgas_per_gal)
{
	auto mstars_disk = galaxy->disk_stars.mass;
	auto mstars_bulge = galaxy->bulge_stars.mass;
	auto mgas_disk = galaxy->disk_gas.mass;
	auto mgas_metals_disk = galaxy->disk_gas.mass_metals;
	auto mBH = galaxy->smbh.mass;
	auto rdisk = galaxy->disk_stars.rscale;
	auto rbulge = galaxy->bulge_stars.rscale;
	auto &molecular_gas = molgas_per_gal.at(galaxy);

	f << mstars_disk << " " << mstars_bulge << " " <<  molecular_gas.m_atom + molecular_gas.m_atom_b
	  << " " << mBH << " " << mgas_metals_disk / mgas_disk << " "
	  << mstars_disk + mstars_bulge << " " << rdisk << " " << rbulge << " "
	  << subhalo->id << " " << subhalo->host_halo->id << "\n";

}

}// namespace shark

