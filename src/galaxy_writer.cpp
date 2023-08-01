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
 * Galaxy writer classes implementations
 */

#include <ctime>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <memory>
#include <numeric>

#include <boost/filesystem.hpp>

#include "hdf5/io/writer.h"
#include "config.h"
#include "cosmology.h"
#include "exceptions.h"
#include "galaxy.h"
#include "galaxy_writer.h"
#include "halo.h"
#include "git_revision.h"
#include "logging.h"
#include "star_formation.h"
#include "subhalo.h"
#include "timer.h"
#include "total_baryon.h"
#include "utils.h"


namespace shark {

GalaxyWriter::GalaxyWriter(ExecutionParameters exec_params, CosmologicalParameters cosmo_params,  CosmologyPtr cosmology, DarkMatterHalosPtr darkmatterhalo, SimulationParameters sim_params, AGNFeedbackParameters agn_params):
	exec_params(std::move(exec_params)),
	cosmo_params(std::move(cosmo_params)),
	cosmology(std::move(cosmology)),
	darkmatterhalo(std::move(darkmatterhalo)),
	sim_params(std::move(sim_params)),
	agn_params(std::move(agn_params))
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

void HDF5GalaxyWriter::write(int snapshot, const std::vector<HaloPtr> &halos, TotalBaryon &AllBaryons, const molgas_per_galaxy &molgas_per_gal)
{
	hdf5::Writer file(get_output_directory(snapshot) + "/galaxies.hdf5");
	write_header(file, snapshot);
	write_galaxies(file, snapshot, halos, molgas_per_gal);
	write_global_properties(file, snapshot, AllBaryons);
	write_sf_histories(snapshot, halos);
	write_bh_histories(snapshot, halos);

}

void HDF5GalaxyWriter::write_header(hdf5::Writer &file, int snapshot){

	std::string comment;

	comment = "the shark version";
	file.write_dataset("run_info/shark_version", std::string(SHARK_VERSION), comment);
	comment = "the git revision of shark used to produce this data";
	file.write_dataset("run_info/shark_git_revision", git_sha1(), comment);
	comment = "whether this shark instance had uncommitted local changes";
	file.write_dataset("run_info/shark_git_has_local_changes", git_has_local_changes(), comment);

	comment = "number of batches analysed";
	file.write_dataset("run_info/batches", exec_params.simulation_batches, comment);

	comment = "accuracy applied when solving the ODE system of the physical model.";
	file.write_dataset("run_info/ode_solver_precision", exec_params.ode_solver_precision, comment);

	comment = "boolean parameter that sets whether the code ignores subhalos that have no descendants.";
	file.write_dataset("run_info/skip_missing_descendants", exec_params.skip_missing_descendants, comment);

	comment = "output snapshot";
	file.write_dataset("run_info/snapshot", snapshot, comment);

	comment = "output redshift";
	file.write_dataset("run_info/redshift", sim_params.redshifts[snapshot], comment);

	comment = "time at which this shark execution started";
	char time_str[20];
	std::strftime(time_str, sizeof(time_str), "%Y-%m-%dT%H:%M:%S", std::gmtime(&exec_params.starting_time));
	file.write_dataset("run_info/timestamp", std::string(time_str), comment);

	file.write_attribute("run_info/model_name", exec_params.name_model);

	comment = "The seed value used in the random number engines";
	file.write_dataset("run_info/seed", exec_params.seed, comment);

	// Calculate effective volume of the run
	float volume = sim_params.volume * exec_params.simulation_batches.size();

	comment = "effective volume of this run [(cMpc/h)^3]";
	file.write_dataset("run_info/effective_volume", volume, comment);

	comment = "dark matter particle mass of this simulation [Msun/h]";
	file.write_dataset("run_info/particle_mass", sim_params.particle_mass, comment);

	comment = "Box side size of the full simulated volume [Mpc/h]";
	file.write_dataset("run_info/lbox", sim_params.lbox, comment);

	comment = "Total number of subvolumes in which the simulated box was divided into";
	file.write_dataset("run_info/tot_n_subvolumes", sim_params.tot_nsubvols, comment);

	// Write cosmological parameters

	comment = "omega matter assumed in simulation";
	file.write_dataset("cosmology/omega_m", cosmo_params.OmegaM, comment);

	comment = "omega baryon assumed in simulation";
	file.write_dataset("cosmology/omega_b", cosmo_params.OmegaB, comment);

	comment = "omega lambda assumed in simulation";
	file.write_dataset("cosmology/omega_l", cosmo_params.OmegaL, comment);

	comment = "scalar spectral index assumed in simulation";
	file.write_dataset("cosmology/n_s", cosmo_params.n_s, comment);

	comment = "fluctuation amplitude at 8 Mpc/h";
	file.write_dataset("cosmology/sigma8", cosmo_params.sigma8, comment);

	comment = "normalization of hubble parameter H0 = h * 100 (km/s)/Mpc";
	file.write_dataset("cosmology/h", cosmo_params.Hubble_h, comment);
}

template<typename T>
static inline
std::size_t report_vsize(const std::vector<T> &v, std::ostringstream &os, const char *name) {
	const std::size_t amount = sizeof(T) * v.size();
	os << " " << name << ": " << memory_amount(amount);
	return amount;
}

void HDF5GalaxyWriter::write_galaxies(hdf5::Writer &file, int snapshot, const std::vector<HaloPtr> &halos, const molgas_per_galaxy &molgas_per_gal){

	Timer t;


	using std::string;
	using std::vector;

	string comment;

	// compute universe age at this redshift:
	double age_uni = std::abs(cosmology->convert_redshift_to_age(sim_params.redshifts[snapshot]));

	// Crate all subhalo properties to write.

	vector<Subhalo::id_t> descendant_id;
	vector<int> main;
	vector<Subhalo::id_t> id;
	vector<Halo::id_t> host_id;
	vector<Galaxy::id_t> id_galaxy;
	vector<Galaxy::id_t> descendant_id_galaxy;
	vector<Halo::id_t> halo_id;

	// Create all galaxies properties to write
	vector<float> mstars_disk;
	vector<float> mstars_bulge;
	vector<float> mstars_burst_mergers;
	vector<float> mstars_burst_diskinstabilities;
	vector<float> mstars_bulge_mergers_assembly;
	vector<float> mstars_bulge_diskins_assembly;
	vector<float> mstars_stripped;
	vector<float> mgas_disk;
	vector<float> mgas_bulge;
	vector<float> mgas_stripped;

	vector<float> mstars_metals_disk;
	vector<float> mstars_metals_bulge;
	vector<float> mstars_metals_burst_mergers;
	vector<float> mstars_metals_burst_diskinstabilities;
	vector<float> mstars_metals_bulge_mergers_assembly;
	vector<float> mstars_metals_bulge_diskins_assembly;
	vector<float> mstars_metals_stripped;
	vector<float> mgas_metals_disk;
	vector<float> mgas_metals_bulge;
	vector<float> mgas_stripped_metals;

	vector<float> mmol_disk;
	vector<float> mmol_bulge;
	vector<float> matom_disk;
	vector<float> matom_bulge;

	vector<float> mBH;
	vector<float> mBH_assembly;
	vector<float> mBH_acc_hh;
	vector<float> mBH_acc_sb;
	vector<float> bh_spin;

	vector<float> sfr_disk;
	vector<float> sfr_burst;
	vector<float> sfr_burst_mergers;
	vector<float> sfr_burst_diskins;
	vector<float> mean_stellar_age;

	vector<float> rdisk_gas;
	vector<float> rbulge_gas;
	vector<float> r_stripped_ism;
	vector<float> sAM_disk_gas;
	vector<float> sAM_disk_gas_atom;
	vector<float> sAM_disk_gas_mol;
	vector<float> sAM_bulge_gas;

	vector<float> rdisk_star;
	vector<float> rbulge_star;
	vector<float> sAM_disk_star;
	vector<float> sAM_bulge_star;

	vector<float> redshift_of_merger;

	vector<float> mhot;
	vector<float> mhot_metals;

	vector<float> mreheated;
	vector<float> mreheated_metals;

	vector<float> mhot_stripped;
	vector<float> mhot_stripped_metals;
	vector<float> r_stripped;

	vector<float> stellar_halo;
	vector<float> stellar_halo_metals;
	vector<float> mean_stellar_mass_galaxies_ihsc;

	vector<float> mlost;
	vector<float> mlost_metals;

	vector<float> cooling_rate;
	vector<int> on_hydrostatic_eq;

	vector<float> mvir_hosthalo;
	vector<float> mvir_subhalo;
	vector<float> vmax_subhalo;
	vector<float> vvir_hosthalo;
	vector<float> vvir_subhalo;
	vector<float> mvir_infall_subhalo;

	vector<float> cnfw_subhalo;
	vector<float> lambda_subhalo;

	vector<float> halo_m;
	vector<float> halo_v;
	vector<float> halo_lambda;
	vector<float> halo_concentration;
	vector<float> age_80_halo;
	vector<float> age_50_halo;
	vector<float> halo_final_m;

	vector<float> infall_time_subhalo;

	vector<float> position_x;
	vector<float> position_y;
	vector<float> position_z;

	vector<float> velocity_x;
	vector<float> velocity_y;
	vector<float> velocity_z;

	vector<float> L_x;
	vector<float> L_y;
	vector<float> L_z;

	vector<float> L_x_subhalo;
	vector<float> L_y_subhalo;
	vector<float> L_z_subhalo;

	vector<int> type;

	vector<Halo::id_t> id_halo;
	vector<Halo::id_t> id_halo_tree;
	vector<Subhalo::id_t> id_subhalo;
	vector<Subhalo::id_t> id_subhalo_tree;

	Halo::id_t j = 1;
	Subhalo::id_t i = 1;

	// Loop over all halos and subhalos to write galaxy properties
	for (auto &halo: halos){


		// assign properties of host halo
		auto mhalo = halo->Mvir;
		auto vhalo = halo->Vvir;

		halo_m.push_back(mhalo);
		halo_v.push_back(vhalo);
		halo_lambda.push_back(halo->lambda);
		halo_concentration.push_back(halo->concentration);
		age_80_halo.push_back(halo->age_80);
		age_50_halo.push_back(halo->age_50);
		halo_id.push_back(halo->id);
		halo_final_m.push_back(halo->final_halo()->Mvir);

		for (auto &subhalo: halo->all_subhalos()){

			host_id.push_back(halo->id);

			// assign properties of host subhalo (note that if these subhalos have descendants, then we assign those properties)
			auto msubhalo = subhalo->Mvir;
			auto cnfw     = subhalo->concentration;
			auto lambda   = subhalo->lambda;
			auto vvir_sh  = subhalo->Vvir;

			// Assign baryon properties of subhalo (note that here we use the subhalo as galaxies and baryons have not yet been transferred to the descendant)
			auto hot_subhalo = subhalo->hot_halo_gas;
			auto cold_subhalo = subhalo->cold_halo_gas;
			auto reheated_subhalo = subhalo->ejected_galaxy_gas;
			auto lost_subhalo = subhalo->lost_galaxy_gas;
			auto stellarhalo = subhalo->stellar_halo;
			auto msub_infall = subhalo->Mvir_infall;
			auto mmeanstellarhalo = subhalo->mean_galaxy_making_stellar_halo;
			auto stripped_subhalo = subhalo->hot_halo_gas_stripped;
			auto r_rps_halo = subhalo->hot_halo_gas_r_rps;

			descendant_id.push_back(subhalo->descendant_id);
			infall_time_subhalo.push_back(subhalo->infall_t);

			int m = 0;
			if(subhalo->main_progenitor){
				m = 1;
			}
			main.push_back(m);
			id.push_back(subhalo->id);

			L_x_subhalo.push_back(subhalo->L.x);
			L_y_subhalo.push_back(subhalo->L.y);
			L_z_subhalo.push_back(subhalo->L.z);

			for (auto &galaxy: subhalo->galaxies){
				//ignore this galaxy if it will appear for the first time in the coming snapshot.
				if(galaxy.birth_snapshot == snapshot) continue;

				if(halo->hydrostatic_eq){
					on_hydrostatic_eq.push_back(1);
				}
				else{
					on_hydrostatic_eq.push_back(0);
				}


				id_halo_tree.push_back(halo->id);
				id_subhalo_tree.push_back(subhalo->id);

				//Calculate molecular gas mass of disk and bulge, and specific angular momentum in atomic/molecular disk.
				auto &molecular_gas = molgas_per_gal.at(galaxy.id);
				// Gas components separated into HI and H2.
				mmol_disk.push_back(molecular_gas.m_mol);
				mmol_bulge.push_back(molecular_gas.m_mol_b);
				matom_disk.push_back(molecular_gas.m_atom);
				matom_bulge.push_back(molecular_gas.m_atom_b);

				// Stellar components
				mstars_disk.push_back(galaxy.disk_stars.mass);
				mstars_bulge.push_back(galaxy.bulge_stars.mass);
				mstars_burst_mergers.push_back(galaxy.galaxymergers_burst_stars.mass);
				mstars_bulge_mergers_assembly.push_back(galaxy.galaxymergers_assembly_stars.mass);
				mstars_burst_diskinstabilities.push_back(galaxy.diskinstabilities_burst_stars.mass);
				mstars_bulge_diskins_assembly.push_back(galaxy.diskinstabilities_assembly_stars.mass);
				mstars_stripped.push_back(galaxy.stars_tidal_stripped.mass);
				auto age = 0;
				if(galaxy.total_stellar_mass_ever_formed > 0){
					age = age_uni - galaxy.mean_stellar_age / galaxy.total_stellar_mass_ever_formed;
				}
				mean_stellar_age.push_back(age);

				// Gas components
				mgas_disk.push_back(galaxy.disk_gas.mass);
				mgas_bulge.push_back(galaxy.bulge_gas.mass);
				mgas_stripped.push_back(galaxy.ram_pressure_stripped_gas.mass);

				// Metals of the stellar components.
				mstars_metals_disk.push_back(galaxy.disk_stars.mass_metals);
				mstars_metals_bulge.push_back(galaxy.bulge_stars.mass_metals);
				mstars_metals_burst_mergers.push_back(galaxy.galaxymergers_burst_stars.mass_metals);
				mstars_metals_bulge_mergers_assembly.push_back(galaxy.galaxymergers_assembly_stars.mass_metals);
				mstars_metals_burst_diskinstabilities.push_back(galaxy.diskinstabilities_burst_stars.mass);
				mstars_metals_bulge_diskins_assembly.push_back(galaxy.diskinstabilities_burst_stars.mass_metals);
				mstars_metals_stripped.push_back(galaxy.stars_tidal_stripped.mass_metals);

				// Metals of the gas components.
				mgas_metals_disk.push_back(galaxy.disk_gas.mass_metals);
				mgas_metals_bulge.push_back(galaxy.bulge_gas.mass_metals);
				mgas_stripped_metals.push_back(galaxy.ram_pressure_stripped_gas.mass_metals);

				// SFRs in disks and bulges.
				sfr_disk.push_back(galaxy.sfr_disk);
				sfr_burst.push_back(galaxy.sfr_bulge_mergers + galaxy.sfr_bulge_diskins);
				sfr_burst_mergers.push_back(galaxy.sfr_bulge_mergers);
				sfr_burst_diskins.push_back(galaxy.sfr_bulge_diskins);

				// Black hole properties.
				mBH.push_back(galaxy.smbh.mass);
				mBH_assembly.push_back(galaxy.smbh.massembly);
				mBH_acc_hh.push_back(galaxy.smbh.macc_hh);
				mBH_acc_sb.push_back(galaxy.smbh.macc_sb);
				bh_spin.push_back(galaxy.smbh.spin);

				// Sizes and specific angular momentum of disks and bulges.
				rdisk_gas.push_back(galaxy.disk_gas.rscale);
				rbulge_gas.push_back(galaxy.bulge_gas.rscale);
				r_stripped_ism.push_back(galaxy.r_rps);
				sAM_disk_gas.push_back(galaxy.disk_gas.sAM);
				sAM_disk_gas_atom.push_back(molecular_gas.j_atom);
				sAM_disk_gas_mol.push_back(molecular_gas.j_mol);
				sAM_bulge_gas.push_back(galaxy.bulge_gas.sAM);

				rdisk_star.push_back(galaxy.disk_stars.rscale);
				rbulge_star.push_back(galaxy.bulge_stars.rscale);
				sAM_disk_star.push_back(galaxy.disk_stars.sAM);
				sAM_bulge_star.push_back(galaxy.bulge_stars.sAM);

				// Halo properties below.
				double mhot_gal = 0;
				double mzhot_gal = 0;
				double mreheat = 0;
				double mzreheat = 0;
				double lostm = 0;
				double lostzm = 0;
				double mstellarhalo = 0;
				double mzstellarhalo = 0;
				double ms_mean_stellarhalo = 0;

				double rcool = 0;
				double mhalo_stripped = 0;
				double mhalo_stripped_metals = 0;
				double r_rps_subhalo = 0;
				int t = galaxy.galaxy_type;

				// Define cooling rate, halo gas and related properties only if subhalo is not a type 2.
				if(galaxy.galaxy_type != Galaxy::TYPE2){
					rcool		= subhalo->cooling_rate;
					mhot_gal	= hot_subhalo.mass + cold_subhalo.mass;
					mzhot_gal	= hot_subhalo.mass_metals + cold_subhalo.mass_metals;
					mreheat		= reheated_subhalo.mass;
					mzreheat 	= reheated_subhalo.mass_metals;
					lostm 		= lost_subhalo.mass;
					lostzm 		= lost_subhalo.mass_metals;
				}

				// the properties below can only be > 0 for type 1 galaxies.
				if(galaxy.galaxy_type == Galaxy::TYPE1){
					mhalo_stripped			= stripped_subhalo.mass;
					mhalo_stripped_metals 	= stripped_subhalo.mass_metals;
					r_rps_subhalo 			= r_rps_halo;
				}

				//stellar halo is >0 only in central galaxies.
				if(galaxy.galaxy_type == Galaxy::CENTRAL){
					mstellarhalo = stellarhalo.mass;
					mzstellarhalo = stellarhalo.mass_metals;
					if(stellarhalo.mass > 0){
						ms_mean_stellarhalo = mmeanstellarhalo / stellarhalo.mass;
					}
					else{
						ms_mean_stellarhalo = 0;
					}
				}

				cooling_rate.push_back(rcool);
				mhot_stripped.push_back(mhalo_stripped);
				mhot_stripped_metals.push_back(mhalo_stripped_metals);
				r_stripped.push_back(r_rps_subhalo);

				mhot.push_back(mhot_gal);
				mhot_metals.push_back(mzhot_gal);
				mreheated.push_back(mreheat);
				mreheated_metals.push_back(mzreheat);
				mlost.push_back(lostm);
				mlost_metals.push_back(lostzm);

				stellar_halo.push_back(mstellarhalo);
				stellar_halo_metals.push_back(mzstellarhalo);
				mean_stellar_mass_galaxies_ihsc.push_back(ms_mean_stellarhalo);

				mvir_hosthalo.push_back(mhalo);
				vvir_hosthalo.push_back(vhalo);

				double mvir_gal = 0 ;
				double c_sub = 0;
				double l_sub = 0;
				double m_infall = 0;
				xyz<float> pos;
				xyz<float> vel;
				xyz<float> L;

				if(galaxy.galaxy_type == Galaxy::CENTRAL || galaxy.galaxy_type == Galaxy::TYPE1){
					mvir_gal = msubhalo;
					c_sub    = cnfw;
					l_sub    = lambda;
					pos      = subhalo->position;
					vel      = subhalo->velocity;
					L        = subhalo->L.unit() * galaxy.angular_momentum();
					vvir_subhalo.push_back(vvir_sh);
					mvir_subhalo.push_back(mvir_gal);
					cnfw_subhalo.push_back(c_sub);
					lambda_subhalo.push_back(l_sub);
					redshift_of_merger.push_back(-1);
					if(galaxy.descendant_id < 0 && snapshot < sim_params.max_snapshot){
						galaxy.descendant_id = galaxy.id;
					}
					if(galaxy.galaxy_type == Galaxy::TYPE1){
						m_infall = msub_infall;
					}
				}
				else{
					// In case of type 2 galaxies assign negative positions, velocities and angular momentum.
					darkmatterhalo->generate_random_orbits(pos, vel, L, galaxy.angular_momentum(), halo, galaxy);
					mvir_subhalo.push_back(galaxy.msubhalo_type2);
					cnfw_subhalo.push_back(galaxy.concentration_type2);
					lambda_subhalo.push_back(galaxy.lambda_type2);
					vvir_subhalo.push_back(galaxy.vvir_type2);

					// calculate the age of the universe by the time this galaxy will merge.
					double tmerge  = cosmology->convert_redshift_to_age(sim_params.redshifts[snapshot-1]) + galaxy.tmerge;
					double redshift_merger = cosmology->convert_age_to_redshift_lcdm(tmerge);
					redshift_of_merger.push_back(redshift_merger);

					if(galaxy.descendant_id < 0 ){
						galaxy.descendant_id = galaxy.id;
					}

				}
				mvir_infall_subhalo.push_back(m_infall);

				//force the descendant Id to be = -1 if this is the last snapshot. If not, check that all descendant_ids are positive.
				if(snapshot == sim_params.max_snapshot){
					galaxy.descendant_id = -1;
				}
				else if (galaxy.descendant_id < 0){
					std::ostringstream os;
					os << "Descendant_id of galaxy to be written is negative";
					throw invalid_argument(os.str());
				}

				id_galaxy.push_back(galaxy.id);
				descendant_id_galaxy.push_back(galaxy.descendant_id);

				vmax_subhalo.push_back(galaxy.vmax);

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
	REPORT(halo_id);
	REPORT(id_galaxy);
	REPORT(descendant_id_galaxy);
	REPORT(host_id);
	REPORT(halo_m);
	REPORT(halo_v);
	REPORT(halo_lambda);
	REPORT(halo_concentration);
	REPORT(age_50_halo);
	REPORT(age_80_halo);
	REPORT(halo_final_m);
	REPORT(infall_time_subhalo);
	REPORT(mstars_disk);
	REPORT(mstars_bulge);
	REPORT(mstars_burst_mergers);
	REPORT(mstars_burst_diskinstabilities);
	REPORT(mstars_bulge_mergers_assembly);
	REPORT(mstars_bulge_diskins_assembly);
	REPORT(mstars_stripped);
	REPORT(mgas_disk);
	REPORT(mgas_bulge);
	REPORT(mgas_stripped);
	REPORT(mstars_metals_disk);
	REPORT(mstars_metals_bulge);
	REPORT(mstars_metals_burst_mergers);
	REPORT(mstars_metals_burst_diskinstabilities);
	REPORT(mstars_metals_bulge_mergers_assembly);
	REPORT(mstars_metals_bulge_diskins_assembly);
	REPORT(mstars_metals_stripped);
	REPORT(mean_stellar_age);
	REPORT(mgas_metals_disk);
	REPORT(mgas_metals_bulge);
	REPORT(mgas_stripped_metals);
	REPORT(mmol_disk);
	REPORT(mmol_bulge);
	REPORT(matom_disk);
	REPORT(matom_bulge);
	REPORT(mBH);
	REPORT(mBH_assembly);
	REPORT(mBH_acc_hh);
	REPORT(mBH_acc_sb);
	REPORT(bh_spin);
	REPORT(sfr_disk);
	REPORT(sfr_burst);
	REPORT(sfr_burst_mergers);
	REPORT(sfr_burst_diskins);
	REPORT(rdisk_gas);
	REPORT(rbulge_gas);
	REPORT(r_stripped_ism);
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
	REPORT(mlost);
	REPORT(mlost_metals);
	REPORT(stellar_halo);
	REPORT(stellar_halo_metals);
	REPORT(mean_stellar_mass_galaxies_ihsc);
	REPORT(cooling_rate);
	REPORT(on_hydrostatic_eq);
	REPORT(mhot_stripped);
	REPORT(mhot_stripped_metals);
	REPORT(mvir_hosthalo);
	REPORT(mvir_subhalo);
	REPORT(vmax_subhalo);
	REPORT(vvir_hosthalo);
	REPORT(vvir_subhalo);
	REPORT(cnfw_subhalo);
	REPORT(r_stripped);
	REPORT(lambda_subhalo);
	REPORT(mvir_infall_subhalo);
	REPORT(L_x_subhalo);
	REPORT(L_y_subhalo);
	REPORT(L_z_subhalo);
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
	if (LOG_ENABLED(debug)) {
		LOG(debug) << "Detailed amounts follow: " << os.str();
	}

	LOG(info) << "Galaxies pivoted and memory reported in " << t;

	t = Timer();

	//Write halo properties.

	comment = "halo id in the tree (unique to entire halo catalogue)";
	file.write_dataset("halo/halo_id", halo_id, comment);

	comment = "virial mass of halo [Msun/h]";
	file.write_dataset("halo/mvir", halo_m, comment);

	comment = "virial velocity of halo [km/s]";
	file.write_dataset("halo/vvir", halo_m, comment);

	comment = "halo concentration";
	file.write_dataset("halo/concentration", halo_concentration, comment);

	comment = "halo spin";
	file.write_dataset("halo/lambda", halo_lambda, comment);

	comment = "redshift at which the halo had 80% of its current mass";
	file.write_dataset("halo/age_80", age_80_halo, comment);

	comment = "redshift at which the halo had 50% of its current mass";
	file.write_dataset("halo/age_50", age_50_halo, comment);

	comment = "virial mass of the halo in which this halo will end up in by z=0 [Msun/h]";
	file.write_dataset("halo/final_z0_mvir", halo_final_m, comment);

	//Write subhalo properties.
	comment = "Subhalo id";
	file.write_dataset("subhalo/id", id, comment);

	comment = "=1 if subhalo is the main progenitor' =0 otherwise.";
	file.write_dataset("subhalo/main_progenitor", main, comment);

	comment = "id of the subhalo that is the descendant of this subhalo";
	file.write_dataset("subhalo/descendant_id", descendant_id, comment);

	comment = "id of the host halo of this subhalo";
	file.write_dataset("subhalo/host_id", host_id, comment);

	comment = "redshift at which the subhalo became a SATELLITE (only well defined for satellite subhalos)";
	file.write_dataset("subhalo/infall_time_subhalo", infall_time_subhalo, comment);

	//Subhalo AM vector
	comment = "total angular momentum component x of subhalo [Msun pMpc km/s]. From VELOCIraptor.";
	file.write_dataset("subhalo/l_x", L_x_subhalo,  comment);
	comment = "total angular momentum component y of galaxy [Msun pMpc km/s]. From VELOCIraptor.";
	file.write_dataset("subhalo/l_y", L_y_subhalo, comment);
	comment = "total angular momentum component z of galaxy [Msun pMpc km/s]. From VELOCIraptor.";
	file.write_dataset("subhalo/l_z", L_z_subhalo, comment);

	//Write galaxy properties.

	comment = "stellar mass in the disk [Msun/h]";
	file.write_dataset("galaxies/mstars_disk", mstars_disk, comment);

	comment = "stellar mass in the bulge [Msun/h]";
	file.write_dataset("galaxies/mstars_bulge", mstars_bulge, comment);

	comment = "stellar mass formed via starbursts driven by galaxy mergers [Msun/h]";
	file.write_dataset("galaxies/mstars_burst_mergers", mstars_burst_mergers, comment);

	comment = "stellar mass formed via starbursts driven by disk instabilities [Msun/h]";
	file.write_dataset("galaxies/mstars_burst_diskinstabilities", mstars_burst_diskinstabilities, comment);

	comment = "stellar mass in the bulge brought via galaxy mergers (but that formed in disks) [Msun/h]";
	file.write_dataset("galaxies/mstars_bulge_mergers_assembly", mstars_bulge_mergers_assembly, comment);

	comment = "stellar mass in the bulge brought via disk instabilities from the disk [Msun/h]";
	file.write_dataset("galaxies/mstars_bulge_diskins_assembly", mstars_bulge_diskins_assembly, comment);

	comment = "stellar mass that was tidally stripped from this galaxy [Msun/h]";
	file.write_dataset("galaxies/mstars_tidally_stripped", mstars_stripped, comment);

	comment = "total gas mass in the disk [Msun/h]";
	file.write_dataset("galaxies/mgas_disk", mgas_disk, comment);

	comment = "gas mass in the bulge [Msun/h]";
	file.write_dataset("galaxies/mgas_bulge", mgas_bulge, comment);

	comment = "gas mass that has been stripped out of the ISM due to ram pressure stripping [Msun/h]";
	file.write_dataset("galaxies/mism_stripped", mgas_stripped, comment);

	comment = "mass of metals locked in stars in the disk [Msun/h]";
	file.write_dataset("galaxies/mstars_metals_disk",mstars_metals_disk, comment);

	comment = "mass of metals locked in stars in the bulge [Msun/h]";
	file.write_dataset("galaxies/mstars_metals_bulge", mstars_metals_bulge, comment);

	comment = "mass of metals locked in stars that formed via starbursts driven by galaxy mergers [Msun/h]";
	file.write_dataset("galaxies/mstars_metals_burst_mergers", mstars_metals_burst_mergers, comment);

	comment = "mass of metals locked in stars that formed via starbursts driven by disk instabilities [Msun/h]";
	file.write_dataset("galaxies/mstars_metals_burst_diskinstabilities", mstars_metals_burst_diskinstabilities, comment);

	comment = "mass of metals locked in stars in the bulge that was brought via galaxy mergers (but that formed in disks) [Msun/h]";
	file.write_dataset("galaxies/mstars_metals_bulge_mergers_assembly", mstars_metals_bulge_mergers_assembly, comment);

	comment = "mass of metals locked in stars in the bulge that was brought via disk instabilities from the disk [Msun/h]";
	file.write_dataset("galaxies/mstars_metals_bulge_diskins_assembly", mstars_metals_bulge_diskins_assembly, comment);

	comment = "mass of metals locked in stars that was tidally stripped from this galaxy [Msun/h]";
	file.write_dataset("galaxies/mstars_metals_tidally_stripped", mstars_metals_stripped, comment);

	comment = "stellar mass-weighted stellar age [Gyr]";
	file.write_dataset("galaxies/mean_stellar_age", mean_stellar_age, comment);

	comment = "mass of metals locked in the gas of the disk [Msun/h]";
	file.write_dataset("galaxies/mgas_metals_disk", mgas_metals_disk, comment);

	comment = "mass of metals locked in the gas of the bulge [Msun/h]";
	file.write_dataset("galaxies/mgas_metals_bulge", mgas_metals_bulge, comment);

	comment = "mass of metals that has been stripped out of the ISM due to ram pressure stripping [Msun/h]";
	file.write_dataset("galaxies/mism_metals_stripped", mgas_stripped_metals, comment);

	comment = "molecular gas mass (helium plus hydrogen) in the disk [Msun/h]";
	file.write_dataset("galaxies/mmol_disk",mmol_disk, comment);

	comment ="molecular gas mass (helium plus hydrogen) in the bulge [Msun/h]";
	file.write_dataset("galaxies/mmol_bulge",mmol_bulge, comment);

	comment = "atomic gas mass (helium plus hydrogen) in the disk [Msun/h]";
	file.write_dataset("galaxies/matom_disk",matom_disk, comment);

	comment ="atomic gas mass (helium plus hydrogen) in the bulge [Msun/h]";
	file.write_dataset("galaxies/matom_bulge",matom_bulge, comment);

	comment = "star formation rate in the disk [Msun/Gyr/h]";
	file.write_dataset("galaxies/sfr_disk", sfr_disk, comment);

	comment = "star formation rate in the bulge [Msun/Gyr/h]";
	file.write_dataset("galaxies/sfr_burst", sfr_burst, comment);

	comment = "star formation rate in the bulge driven by galaxy mergers [Msun/Gyr/h]";
	file.write_dataset("galaxies/sfr_burst_mergers", sfr_burst_mergers, comment);

	comment = "star formation rate in the bulge driven by disk instabilities [Msun/Gyr/h]";
	file.write_dataset("galaxies/sfr_burst_diskins", sfr_burst_diskins, comment);

	comment = "black hole mass [Msun/h]";
	file.write_dataset("galaxies/m_bh", mBH, comment);

	comment = "black hole mass that comes from assembly (BH-BH mergers) [Msun/h]";
	file.write_dataset("galaxies/m_bh_assembly", mBH_assembly, comment);

	comment = "accretion rate onto the black hole during the hot halo mode [Msun/Gyr/h]";
	file.write_dataset("galaxies/bh_accretion_rate_hh", mBH_acc_hh, comment);

	comment = "accretion rate onto the black hole during the starburst mode [Msun/Gyr/h]";
	file.write_dataset("galaxies/bh_accretion_rate_sb", mBH_acc_sb, comment);

	comment = "black hole spin [dimensionless]";
	file.write_dataset("galaxies/bh_spin", bh_spin, comment);

	comment = "half-mass radius of the stellar disk [cMpc/h]";
	file.write_dataset("galaxies/rstar_disk", rdisk_star, comment);

	comment = "half-mass radius of the stellar bulge [cMpc/h]";
	file.write_dataset("galaxies/rstar_bulge", rbulge_star, comment);

	comment = "specific angular momentum of the stellar disk [km/s * cMpc/h]";
	file.write_dataset("galaxies/specific_angular_momentum_disk_star", sAM_disk_star, comment);

	comment = "specific angular momentum of the stellar bulge [km/s * cMpc/h]";
	file.write_dataset("galaxies/specific_angular_momentum_bulge_star", sAM_bulge_star, comment);

	comment = "half-mass radius of the gas disk [cMpc/h]";
	file.write_dataset("galaxies/rgas_disk", rdisk_gas, comment);

	comment = "half-mass radius of the gas bulge [cMpc/h]";
	file.write_dataset("galaxies/rgas_bulge", rbulge_gas, comment);

	comment = "ram pressure stripping radius of the ISM [cMpc/h]";
	file.write_dataset("galaxies/r_ism_stripped", r_stripped_ism, comment);

	comment = "specific angular momentum of the gas disk [km/s * cMpc/h]";
	file.write_dataset("galaxies/specific_angular_momentum_disk_gas", sAM_disk_gas, comment);

	comment = "specific angular momentum of the atomic gas disk [km/s * cMpc/h]";
	file.write_dataset("galaxies/specific_angular_momentum_disk_gas_atom", sAM_disk_gas_atom, comment);

	comment = "specific angular momentum of the molecular gas disk [km/s * cMpc/h]";
	file.write_dataset("galaxies/specific_angular_momentum_disk_gas_mol", sAM_disk_gas_mol, comment);

	comment = "specific angular momentum of the gas bulge [km/s * cMpc/h]";
	file.write_dataset("galaxies/specific_angular_momentum_bulge_gas", sAM_bulge_gas, comment);

	comment = "redshift at which this galaxy will merge onto a central galaxy (only relevant for type 2 galaxies)";
	file.write_dataset("galaxies/redshift_merger", redshift_of_merger, comment);

	comment = "hot gas mass in the halo [Msun/h]";
	file.write_dataset("galaxies/mhot", mhot, comment);

	comment = "mass of metals locked in the hot halo gas [Msun/h]";
	file.write_dataset("galaxies/mhot_metals", mhot_metals, comment);

	comment = "gas mass in the ejected gas component [Msun/h]";
	file.write_dataset("galaxies/mreheated", mreheated, comment);

	comment = "mass of metals locked in the ejected gas component [Msun/h]";
	file.write_dataset("galaxies/mreheated_metals", mreheated_metals, comment);

	comment = "gas mass in the lost gas component - due to QSO feedback [Msun/h]";
	file.write_dataset("galaxies/mlost", mlost, comment);

	comment = "mass of metals locked in the lost gas component - due to QSO feedback [Msun/h]";
	file.write_dataset("galaxies/mlost_metals", mlost_metals, comment);

	comment = "stellar mass in the halo built by tidal stripping [Msun/h]";
	file.write_dataset("galaxies/mstellar_halo", stellar_halo, comment);

	comment = "mass of metals locked up in the stellar halo built by tidal stripping [Msun/h]";
	file.write_dataset("galaxies/mstellar_halo_metals", stellar_halo_metals, comment);

	comment = "mass weighted stellar mass of the galaxies that contributed to building the stellar halo [Msun/h]";
	file.write_dataset("galaxies/mean_mstellar_galaxies_stellarhalo", mean_stellar_mass_galaxies_ihsc, comment);

	comment = "cooling rate of the hot halo component [Msun/Gyr/h].";
	file.write_dataset("galaxies/cooling_rate", cooling_rate, comment);

	comment = "is halo on quasi hydrostatic equilibrium (=1 for true, =0 for false).";
	file.write_dataset("galaxies/on_hydrostatic_eq", on_hydrostatic_eq, comment);

	comment = "gas mass that has been stripped out of this subhalo due to ram pressure stripping [Msun/h].";
	file.write_dataset("galaxies/mhot_stripped", mhot_stripped, comment);

	comment = "mass of metals that has been stripped out of this subhalo due to ram pressure stripping [Msun/h].";
	file.write_dataset("galaxies/mhot_metals_stripped", mhot_stripped_metals, comment);

	comment = "Dark matter mass of the host halo in which this galaxy resides [Msun/h]";
	file.write_dataset("galaxies/mvir_hosthalo", mvir_hosthalo, comment);

	comment = "Dark matter mass of the subhalo in which this galaxy resides [Msun/h]. In the case of type 2 satellites, this corresponds to the mass its subhalo had before disappearing from the subhalo catalogs.";
	file.write_dataset("galaxies/mvir_subhalo", mvir_subhalo, comment);

	comment = "Maximum circular velocity of this galaxy [km/s]";
	file.write_dataset("galaxies/vmax_subhalo", vmax_subhalo, comment);

	comment = "Virial velocity of the dark matter subhalo in which this galaxy resides [km/s]. In the case of type 2 satellites, this corresponds to the virial velocity its subhalo had before disappearing from the subhalo catalogs.";
	file.write_dataset("galaxies/vvir_subhalo", vvir_subhalo, comment);

	comment = "Virial velocity of the dark matter host halo in which this galaxy resides [km/s].";
	file.write_dataset("galaxies/vvir_hosthalo", vvir_hosthalo, comment);

	comment = "ram pressure stripping radius of the halo gas [cMpc/h]";
	file.write_dataset("galaxies/r_halo_stripped", r_stripped, comment);

	comment = "NFW concentration parameter of the dark matter subhalo in which this galaxy resides [dimensionless]. In the case of type 2 satellites, this corresponds to the concentration its subhalo had before disappearing from the subhalo catalogs.";
	file.write_dataset("galaxies/cnfw_subhalo", cnfw_subhalo, comment);

	comment = "Spin parameter of the dark matter subhalo in which this galaxy resides [dimensionless].  In the case of type 2 satellites, this corresponds to the lambda its subhalo had before disappearing from the subhalo catalogs.";
	file.write_dataset("galaxies/lambda_subhalo", lambda_subhalo, comment);

	comment = "Dark matter mass at infall of the host halo in which this galaxy reside when it was last central [Msun/h]";
	file.write_dataset("galaxies/mvir_infall_subhalo", mvir_infall_subhalo, comment);

	//Galaxy position
	comment = "position component x of galaxy [cMpc/h]. In the case of type 2 galaxies, the positions are generated to randomly sample an NFW halo with the concentration of the halo the galaxy lives in.";
	file.write_dataset("galaxies/position_x", position_x, comment);
	comment = "position component y of galaxy [cMpc/h]. In the case of type 2 galaxies, the positions are generated to randomly sample an NFW halo with the concentration of the halo the galaxy lives in.";
	file.write_dataset("galaxies/position_y", position_y, comment);
	comment = "position component z of galaxy [cMpc/h]. In the case of type 2 galaxies, the positions are generated to randomly sample an NFW halo with the concentration of the halo the galaxy lives in.";
	file.write_dataset("galaxies/position_z", position_z, comment);

	//Galaxy velocity
	comment = "peculiar velocity component x of galaxy [km/s]. In the case of type 2 galaxies, the velocity is generated to randomly sample the velocity dispersion of a NFW halo with the concentration of the halo the galaxy lives in.";
	file.write_dataset("galaxies/velocity_x", velocity_x, comment);
	comment = "peculiar velocity component y of galaxy [km/s]. In the case of type 2 galaxies, the velocity is generated to randomly sample the velocity dispersion of a NFW halo with the concentration of the halo the galaxy lives in.";
	file.write_dataset("galaxies/velocity_y", velocity_y, comment);
	comment = "peculiar velocity component z of galaxy [km/s]. In the case of type 2 galaxies, the velocity is generated to randomly sample the velocity dispersion of a NFW halo with the concentration of the halo the galaxy lives in.";
	file.write_dataset("galaxies/velocity_z", velocity_z, comment);

	//Galaxy AM vector
	comment = "total angular momentum component x of galaxy [Msun pMpc km/s]. In the case of type 2 galaxies, the AM vector is randomly oriented.";
	file.write_dataset("galaxies/l_x", L_x,  comment);
	comment = "total angular momentum component y of galaxy [Msun pMpc km/s]. In the case of type 2 galaxies, the AM vector is randomly oriented.";
	file.write_dataset("galaxies/l_y", L_y, comment);
	comment = "total angular momentum component z of galaxy [Msun pMpc km/s]. In the case of type 2 galaxies, the AM vector is randomly oriented.";
	file.write_dataset("galaxies/l_z", L_z, comment);

	//Galaxy type.
	comment = "galaxy type; =0 for centrals; =1 for satellites that reside in well identified subhalos; =2 for orphan satellites";
	file.write_dataset("galaxies/type", type, comment);

	//Galaxy IDs.
	comment = "subhalo ID. Unique to this snapshot.";
	file.write_dataset("galaxies/id_subhalo", id_subhalo, comment);

	comment = "halo ID. Unique to this snapshot.";
	file.write_dataset("galaxies/id_halo", id_halo, comment);

	comment = "galaxy ID. Unique to this galaxy throughout time. If this galaxy never mergers onto a central, then its ID is always the same.";
	file.write_dataset("galaxies/id_galaxy", id_galaxy, comment);

	comment = "descendant galaxy ID. Different to galaxy id only if galaxy is type 2 and merges on the next snapshot.";
	file.write_dataset("galaxies/descendant_id_galaxy", descendant_id_galaxy, comment);

	comment = "subhalo id in the tree (unique to entire halo catalogue).";
	file.write_dataset("galaxies/id_subhalo_tree", id_subhalo_tree, comment);

	comment = "halo id in the tree (unique to entire halo catalogue).";
	file.write_dataset("galaxies/id_halo_tree", id_halo_tree, comment);

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
	file.write_dataset("global/redshifts", redshifts, comment);

	comment = "total cold gas mass (interstellar medium) in the simulated box [Msun/h]";
	file.write_dataset("global/mcold",AllBaryons.get_masses(AllBaryons.mcold), comment);

	comment = "total mass of metals locked in cold gas in the simulated box [Msun/h]";
	file.write_dataset("global/mcold_metals",AllBaryons.get_metals(AllBaryons.mcold), comment);

	comment = "total stellar mass in the simulated box [Msun/h]";
	file.write_dataset("global/mstars",AllBaryons.get_masses(AllBaryons.mstars), comment);

	comment = "total mass of metals locked in stars in the simulated box [Msun/h]";
	file.write_dataset("global/mstars_metals",AllBaryons.get_metals(AllBaryons.mstars), comment);

	comment = "total stellar mass formed via starbursts triggered by galaxy mergers in the simulated box [Msun/h]";
	file.write_dataset("global/mstars_bursts_mergers",AllBaryons.get_masses(AllBaryons.mstars_burst_galaxymergers), comment);

	comment = "total mass of metals locked in stars that formed via starbursts triggered by galaxy mergers in the simulated box [Msun/h]";
	file.write_dataset("global/mstars_metals_bursts_mergers",AllBaryons.get_metals(AllBaryons.mstars_burst_galaxymergers), comment);

	comment = "total stellar mass formed via starbursts triggered by disk instabilities in the simulated box [Msun/h]";
	file.write_dataset("global/mstars_bursts_diskinstabilities",AllBaryons.get_masses(AllBaryons.mstars_burst_diskinstabilities), comment);

	comment = "total mass of metals locked in stars that formed via starbursts triggered by disk instabilities in the simulated box [Msun/h]";
	file.write_dataset("global/mstars_metals_bursts_diskinstabilities",AllBaryons.get_metals(AllBaryons.mstars_burst_diskinstabilities), comment);

	comment = "total atomic gas mass in the simulated box [Msun/h]";
	file.write_dataset("global/m_hi",AllBaryons.get_masses(AllBaryons.mHI), comment);

	comment = "total molecular gas mass in the simulated box [Msun/h]";
	file.write_dataset("global/m_h2",AllBaryons.get_masses(AllBaryons.mH2), comment);

	comment = "total mass locked up in black holes in the simulated box [Msun/h]";
	file.write_dataset("global/m_bh",AllBaryons.get_masses(AllBaryons.mBH), comment);

	comment = "total star formation rate taking place in disks in the simulated box [Msun/Gyr/h]";
	file.write_dataset("global/sfr_quiescent",AllBaryons.SFR_disk, comment);

	comment = "total star formation rate taking place in bulges in the simulated box [Msun/Gyr/h]";
	file.write_dataset("global/sfr_burst",AllBaryons.SFR_bulge, comment);

	comment = "Maximum mass of the SMBHs in the simulated box [Msun/h]";
	file.write_dataset("global/smbh_maximum",AllBaryons.max_BH, comment);

	comment = "number of major mergers taking place in the simulated box at this snapshot.";
	file.write_dataset("global/number_major_mergers", AllBaryons.major_mergers, comment);

	comment = "number of minor mergers taking place in the simulated box at this snapshot.";
	file.write_dataset("global/number_minor_mergers", AllBaryons.minor_mergers, comment);

	comment = "number of disk instability episodes taking place in the simulated box at this snapshot.";
	file.write_dataset("global/number_disk_instabilities", AllBaryons.disk_instabil, comment);

	comment = "total hot gas mass in halos in the simulated box [Msun/h]";
	file.write_dataset("global/mhot_halo",AllBaryons.get_masses(AllBaryons.mhot_halo),comment);

	comment = "total mass of metals in the hot gas mass in halos in the simulated box [Msun/h]";
	file.write_dataset("global/mhot_metals",AllBaryons.get_metals(AllBaryons.mhot_halo), comment);

	comment = "total halo cold gas in the simulated box [Msun/h]";
	file.write_dataset("global/mcold_halo",AllBaryons.get_masses(AllBaryons.mcold_halo), comment);

	comment = "total mass of metals in the halo cold gas mass in the simulated box [Msun/h]";
	file.write_dataset("global/mcold_halo_metals",AllBaryons.get_metals(AllBaryons.mcold_halo), comment);

	comment = "total gas mass ejected from halos due to stellar feedback (and that has not yet been reincorporated) in the simulated box [Msun/h]";
	file.write_dataset("global/mejected_halo",AllBaryons.get_masses(AllBaryons.mejected_halo), comment);

	comment = "total mass of metals in the ejected gas reservoir due to stellar feedback in the simulated box [Msun/h]";
	file.write_dataset("global/mejected_halo_metals",AllBaryons.get_metals(AllBaryons.mejected_halo), comment);

	comment = "total gas mass ejected from halos due to QSO feedback in the simulated box [Msun/h]";
	file.write_dataset("global/mlost_halo",AllBaryons.get_masses(AllBaryons.mlost_halo), comment);

	comment = "total mass of metals in the ejected gas reservoir due to QSO feedback in the simulated box [Msun/h]";
	file.write_dataset("global/mlost_halo_metals",AllBaryons.get_metals(AllBaryons.mlost_halo), comment);

	comment = "total dark matter mass locked up in halos in the simulated box [Msun/h].";
	file.write_dataset("global/m_dm",AllBaryons.get_masses(AllBaryons.mDM), comment);

	comment = "total baryon mass in the simulated box [Msun/h]";
	file.write_dataset("global/mbar_created",baryons_ever_created, comment);

	comment = "total baryons lost in the simulated box [Msun/h] (ideally this should be =0)";
	file.write_dataset("global/mbar_lost", baryons_ever_lost, comment);
}

void HDF5GalaxyWriter::write_sf_histories (int snapshot, const std::vector<HaloPtr> &halos){


	using std::string;
	using std::vector;

	string comment;

	if(exec_params.output_sf_histories){
		if(std::find(exec_params.snapshots_sf_histories.begin(), exec_params.snapshots_sf_histories.end(), snapshot) != exec_params.snapshots_sf_histories.end()){

			Timer sfh_writer_timer;
			hdf5::Writer file_sfh(get_output_directory(snapshot) + "/star_formation_histories.hdf5");

			//Create the vectors that will save the information of the galaxies
			vector<vector<float>> sfhs_disk;
			vector<vector<float>> stellar_mass_disk;
			vector<vector<float>> stellar_metals_disk;

			vector<vector<float>> sfhs_bulge_mergers;
			vector<vector<float>> stellar_mass_bulge_mergers;
			vector<vector<float>> stellar_metals_bulge_mergers;

			vector<vector<float>> sfhs_bulge_diskins;
			vector<vector<float>> stellar_mass_bulge_diskins;
			vector<vector<float>> stellar_metals_bulge_diskins;

			vector<Galaxy::id_t> id_galaxy;

			float defl_value = 0;

			for (auto &halo: halos){
				for (auto &subhalo: halo->all_subhalos()){
					for (auto &galaxy: subhalo->galaxies){
						//ignore this galaxy if it will appear for the first time in the coming snapshot.
						if(galaxy.birth_snapshot == snapshot) continue;

						vector<float> sfh_gal_disk;
						vector<float> star_metals_gal_disk;
						vector<float> sfh_gal_bulge_mergers;
						vector<float> star_metals_gal_bulge_mergers;
						vector<float> sfh_gal_bulge_diskins;
						vector<float> star_metals_gal_bulge_diskins;

						bool star_gal_bulge_exists = false;
						for(int s=sim_params.min_snapshot+1; s <= snapshot; s++) {

							auto it = std::find_if(galaxy.history.begin(), galaxy.history.end(), [s](const HistoryItem &hitem) {
								//information in snapshot corresponds to the end of it, so effectively, when writing, we need to
								//compare to s-1.
								return hitem.snapshot == s-1;
							});

							if (it == galaxy.history.end()) {
								if (star_gal_bulge_exists) {
									std::ostringstream os;
									os << "The history of the StellarMass of the bulge of " << galaxy << " ceased to exist (temporarily). ";
									os << "These are the snapshots for which there is a history item: ";
									vector<int> hsnaps(galaxy.history.size());
									std::transform(galaxy.history.begin(), galaxy.history.end(), hsnaps.begin(), [](const HistoryItem &hitem) {
										return hitem.snapshot;
									});
									std::copy(hsnaps.begin(), hsnaps.end(), std::ostream_iterator<int>(os, " "));
									LOG(warning) << os.str();
									for (auto &item: galaxy.history){
										std::ostringstream os;
										os << "snap history: "<< item.snapshot;
										LOG(warning) << os.str();
									}
								}
								sfh_gal_disk.push_back(defl_value);
								star_metals_gal_disk.push_back(defl_value);

								sfh_gal_bulge_mergers.push_back(defl_value);
								star_metals_gal_bulge_mergers.push_back(defl_value);

								sfh_gal_bulge_diskins.push_back(defl_value);
								star_metals_gal_bulge_diskins.push_back(defl_value);
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

								// assign bulge properties driven by mergers
								sfh_gal_bulge_mergers.push_back(item.sfr_bulge_mergers/constants::GIGA);
								if(item.sfr_bulge_mergers > 0){
									star_metals_gal_bulge_mergers.push_back(item.sfr_z_bulge_mergers/item.sfr_bulge_mergers);
								}
								else{
									star_metals_gal_bulge_mergers.push_back(0);
								}

								// assign bulge properties driven by disk instabilities
								sfh_gal_bulge_diskins.push_back(item.sfr_bulge_diskins/constants::GIGA);
								if(item.sfr_bulge_diskins > 0){
									star_metals_gal_bulge_diskins.push_back(item.sfr_z_bulge_diskins/item.sfr_bulge_diskins);
								}
								else{
									star_metals_gal_bulge_diskins.push_back(0);
								}
							}
						}

						// save galaxies only if they have a stellar mass >0 by the output snapshot.
						if(galaxy.stellar_mass() > 0){
							sfhs_disk.emplace_back(std::move(sfh_gal_disk));
							stellar_metals_disk.emplace_back(std::move(star_metals_gal_disk));

							sfhs_bulge_mergers.emplace_back(std::move(sfh_gal_bulge_mergers));
							stellar_metals_bulge_mergers.emplace_back(std::move(star_metals_gal_bulge_mergers));

							sfhs_bulge_diskins.emplace_back(std::move(sfh_gal_bulge_diskins));
							stellar_metals_bulge_diskins.emplace_back(std::move(star_metals_gal_bulge_diskins));

							id_galaxy.push_back(galaxy.id);
						}
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

			comment = "galaxy ID. Unique to this galaxy throughout time. If this galaxy never mergers onto a central, then its ID is always the same.";
			file_sfh.write_dataset("galaxies/id_galaxy", id_galaxy, comment);

			//Write disk component history.
			comment = "Star formation history of stars formed that by this output time end up in the disk [Msun/yr/h]";
			file_sfh.write_dataset("disks/star_formation_rate_histories", sfhs_disk, comment);

			comment = "Stellar metallicity of the stars formed in a timestep that by this output time ends up in the disk";
			file_sfh.write_dataset("disks/metallicity_histories", stellar_metals_disk, comment);

			//Write bulge component history, for the mass build up due to galaxy mergers.
			comment = "Star formation history of stars formed that by this output time end up in the bulge formed via galaxy mergers [Msun/yr/h]";
			file_sfh.write_dataset("bulges_mergers/star_formation_rate_histories", sfhs_bulge_mergers, comment);

			comment = "Stellar metallicity of the stars formed in a timestep that by this output time ends up in the bulge formed via galaxy mergers";
			file_sfh.write_dataset("bulges_mergers/metallicity_histories", stellar_metals_bulge_mergers, comment);

			//Write bulge component history.
			comment = "Star formation history of stars formed that by this output time end up in the bulge formed via disk instabilities [Msun/yr/h]";
			file_sfh.write_dataset("bulges_diskins/star_formation_rate_histories", sfhs_bulge_diskins, comment);

			comment = "Stellar metallicity of the stars formed in a timestep that by this output time ends up in the bulge formed via disk instabilities";
			file_sfh.write_dataset("bulges_diskins/metallicity_histories", stellar_metals_bulge_diskins, comment);

			comment = "Redshifts of the history outputs";
			file_sfh.write_dataset("redshifts", redshifts, comment);

			comment = "Look back time to mean time between snapshots [Gyr]";
			file_sfh.write_dataset("lbt_mean", age_mean, comment);

			comment = "Time interval covered between snapshots [Gyr]";
			file_sfh.write_dataset("delta_t", delta_t, comment);

			LOG(info) << "Galaxies SFH data written in " << sfh_writer_timer;
		}

	}
}

void HDF5GalaxyWriter::write_bh_histories (int snapshot, const std::vector<HaloPtr> &halos){


	using std::string;
	using std::vector;

	string comment;

	if(exec_params.output_bh_histories){
		if(std::find(exec_params.snapshots_bh_histories.begin(), exec_params.snapshots_bh_histories.end(), snapshot) != exec_params.snapshots_bh_histories.end()){
			hdf5::Writer file_bh(get_output_directory(snapshot) + "/black_hole_histories.hdf5");

			//Create the vectors that will save the information of the galaxies
			vector<vector<float>> bh_mass;
			vector<vector<float>> bh_spin;
			vector<vector<float>> bh_assembly;

			vector<vector<float>> macc_hh;
			vector<vector<float>> macc_sb;

			vector<Galaxy::id_t> id_galaxy;

			float defl_value = 0;

			for (auto &halo: halos){
				for (auto &subhalo: halo->all_subhalos()){
					for (auto &galaxy: subhalo->galaxies){
						//ignore this galaxy if it will appear for the first time in the coming snapshot.
						if(galaxy.birth_snapshot == snapshot) continue;

						vector<float> bh_mass_gal;
						vector<float> bh_spin_gal;
						vector<float> bh_assembly_gal;
						vector<float> macc_hh_gal;
						vector<float> macc_sb_gal;

						for(int s=sim_params.min_snapshot+1; s <= snapshot; s++) {

							auto it = std::find_if(galaxy.bh_history.begin(), galaxy.bh_history.end(), [s](const BHHistoryItem &hitem) {
								//information in snapshot corresponds to the end of it, so effectively, when writing, we need to
								//compare to s-1.
								return hitem.snapshot == s-1;
							});

							if (it == galaxy.bh_history.end()) {
								bh_mass_gal.push_back(defl_value);
								bh_spin_gal.push_back(defl_value);
								bh_assembly_gal.push_back(defl_value);

								macc_hh_gal.push_back(defl_value);
								macc_sb_gal.push_back(defl_value);
							}
							else {
								auto item = *it;
								// assign disk properties
								macc_hh_gal.push_back(item.macc_hh/constants::GIGA);
								macc_sb_gal.push_back(item.macc_sb/constants::GIGA);
								bh_mass_gal.push_back(item.mbh);
								bh_spin_gal.push_back(item.spin);
								bh_assembly_gal.push_back(item.massembly);	
							}
						}

						// save galaxies only if they have a bh mass >seed BH by the output snapshot.
						if(galaxy.smbh.mass > agn_params.mseed){
							macc_hh.emplace_back(std::move(macc_hh_gal));
							macc_sb.emplace_back(std::move(macc_sb_gal));

							bh_mass.emplace_back(std::move(bh_mass_gal));
							bh_spin.emplace_back(std::move(bh_spin_gal));
							bh_assembly.emplace_back(std::move(bh_assembly_gal));

							id_galaxy.push_back(galaxy.id);
						}
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
			write_header(file_bh, snapshot);

			comment = "galaxy ID. Unique to this galaxy throughout time. If this galaxy never mergers onto a central, then its ID is always the same.";
			file_bh.write_dataset("galaxies/id_galaxy", id_galaxy, comment);

			//Write accreion rates
			comment = "Black hole accretion rate due to hot halo cooling [Msun/yr/h].";
			file_bh.write_dataset("galaxies/bh_accretion_rate_hh_history", macc_hh, comment);

			comment = "Black hole accretion rate due to starbursts [Msun/yr/h].";
			file_bh.write_dataset("galaxies/bh_accretion_rate_sb_history", macc_sb, comment);

			//Write masses and spin
			comment = "Black hole mass history (cumulative) [Msun/h].";
                        file_bh.write_dataset("galaxies/m_bh_history", bh_mass, comment);

			comment = "Black hole mass history coming from BH-BH mergers (cumulative) [Msun/h].";
                        file_bh.write_dataset("galaxies/m_bh_assembly_history", bh_assembly, comment);

			comment = "Black hole spin history [dimensionless].";
                        file_bh.write_dataset("galaxies/bh_spin", bh_spin, comment);

			comment = "Redshifts of the history outputs";
			file_bh.write_dataset("redshifts", redshifts, comment);

			comment = "Look back time to mean time between snapshots [Gyr]";
			file_bh.write_dataset("lbt_mean", age_mean, comment);

			comment = "Time interval covered between snapshots [Gyr]";
			file_bh.write_dataset("delta_t", delta_t, comment);

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

void ASCIIGalaxyWriter::write_galaxy(const Galaxy &galaxy, const SubhaloPtr &subhalo, int snapshot, std::ofstream &f, const molgas_per_galaxy &molgas_per_gal)
{
	auto mstars_disk = galaxy.disk_stars.mass;
	auto mstars_bulge = galaxy.bulge_stars.mass;
	auto mgas_disk = galaxy.disk_gas.mass;
	auto mgas_metals_disk = galaxy.disk_gas.mass_metals;
	auto mBH = galaxy.smbh.mass;
	auto rdisk = galaxy.disk_stars.rscale;
	auto rbulge = galaxy.bulge_stars.rscale;
	auto &molecular_gas = molgas_per_gal.at(galaxy.id);

	f << mstars_disk << " " << mstars_bulge << " " <<  molecular_gas.m_atom + molecular_gas.m_atom_b
	  << " " << mBH << " " << mgas_metals_disk / mgas_disk << " "
	  << mstars_disk + mstars_bulge << " " << rdisk << " " << rbulge << " "
	  << subhalo->id << " " << subhalo->host_halo->id << "\n";

}

}// namespace shark

