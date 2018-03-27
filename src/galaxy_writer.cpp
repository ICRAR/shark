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

void HDF5GalaxyWriter::write(int snapshot, const std::vector<HaloPtr> &halos, TotalBaryon &AllBaryons){

	using std::string;
	using std::vector;

	string comment;

	hdf5::Writer file(get_output_directory(snapshot) + "/galaxies.hdf5");

	//Write header
	write_header(file, snapshot);

	//Write galaxies
	write_galaxies(file, snapshot, halos);

	//Write total baryon components
	write_global_properties(file, snapshot, AllBaryons);

	// Write star formation histories.
	write_histories(snapshot, halos);

}

void HDF5GalaxyWriter::write_header(hdf5::Writer file, int snapshot){

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
std::size_t report_vsize(std::vector<T> v, std::ostringstream &os, const char *name) {
	const std::size_t amount = sizeof(T) * v.size();
	os << " " << name << ": " << memory_amount(amount);
	return amount;
};

void HDF5GalaxyWriter::write_galaxies(hdf5::Writer file, int snapshot, const std::vector<HaloPtr> &halos){



	using std::string;
	using std::vector;

	string comment;

	// Crate all subhalo properties to write.

	vector<long> descendant_id;
	vector<int> main;
	vector<long> id;
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

	long j = 1;
	long gal_id = 1;
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

				// Gas components separated into HI and H2.
				mmol_disk.push_back(m_mol);
				mmol_bulge.push_back(m_mol_b);
				matom_disk.push_back(m_atom);
				matom_bulge.push_back(m_atom_b);

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
				rdisk.push_back(galaxy->disk_stars.rscale);
				rbulge.push_back(galaxy->bulge_stars.rscale);
				sAM_disk.push_back(galaxy->disk_stars.sAM);
				sAM_bulge.push_back(galaxy->bulge_stars.sAM);

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

				// Galaxy position and velocity.
				position_x.push_back(pos.x);
				position_y.push_back(pos.y);
				position_z.push_back(pos.z);

				velocity_x.push_back(vel.x);
				velocity_y.push_back(vel.y);
				velocity_z.push_back(vel.z);

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
	REPORT(rdisk);
	REPORT(rbulge);
	REPORT(sAM_disk);
	REPORT(sAM_bulge);
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
	REPORT(position_x);
	REPORT(position_y);
	REPORT(position_z);
	REPORT(velocity_x);
	REPORT(velocity_y);
	REPORT(velocity_z);
	REPORT(type);
	REPORT(id_halo);
	REPORT(id_subhalo);

	LOG(info) << "Total amount of memory used by the writing process: " << memory_amount(total);
	LOG(debug) << "Detailed amounts follow: " << os.str();

	//Write subhalo properties.
	comment = "Subhalo id";
	file.write_dataset("Subhalo/id", id, comment);

	comment = "=1 if subhalo is the main progenitor' =0 otherwise.";
	file.write_dataset("Subhalo/main_progenitor", main, comment);

	comment = "id of the subhalo that is the descendant of this subhalo";
	file.write_dataset("Subhalo/descendant_id", descendant_id, comment);

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

	comment = "half-mass radius of the disk [cMpc/h]";
	file.write_dataset("Galaxies/rdisk", rdisk, comment);

	comment = "half-mass radius of the bulge [cMpc/h]";
	file.write_dataset("Galaxies/rbulge", rbulge, comment);

	comment = "specific angular momentum of the disk [km/s * cMpc/h]";
	file.write_dataset("Galaxies/specific_angular_momentum_disk", sAM_disk, comment);

	comment = "specific angular momentum of the bulge [km/s * cMpc/h]";
	file.write_dataset("Galaxies/specific_angular_momentum_bulge", sAM_bulge, comment);

	comment = "hot gas mass in the halo [Msun/h]";
	file.write_dataset("Galaxies/mhot", mhot, comment);

	comment = "mass of metals locked in the hot halo gas [Msun/h]";
	file.write_dataset("Galaxies/mhot_metals", mhot_metals, comment);

	comment = "gas mass in the ejected gas component [Msun/h]";
	file.write_dataset("Galaxies/mreheated", mreheated, comment);

	comment = "mass of metals locked in the ejected gas component [Msun/h]";
	file.write_dataset("Galaxies/mreheated_metals", mreheated_metals, comment);

	comment = "cooling rate of the hot halo compoent [Msun/Gyr/h].";
	file.write_dataset("Galaxies/cooling_rate", cooling_rate, comment);

	comment = "Dark matter mass of the host halo in which this galaxy resides [Msun/h]";
	file.write_dataset("Galaxies/mvir_hosthalo", mvir_hosthalo, comment);

	comment = "Dark matter mass of the subhalo in which this galaxy resides [Msun/h]";
	file.write_dataset("Galaxies/mvir_subhalo", mvir_subhalo, comment);

	comment = "Maximum circular velocity of the dark matter subhalo in which this galaxy resides [km/s]";
	file.write_dataset("Galaxies/vmax_subhalo", vmax_subhalo, comment);

	comment = "Virial velocity of the dark matter halo in which this galaxy resides [km/s]";
	file.write_dataset("Galaxies/vvir_hosthalo", vvir_hosthalo, comment);

	comment = "NFW concentration parameter of the dark matter halo in which this galaxy resides [dimensionless]";
	file.write_dataset("Galaxies/cnfw_subhalo", cnfw_subhalo, comment);

	//Galaxy position
	comment = "position component x of galaxy [cMpc/h]";
	file.write_dataset("Galaxies/position_x", position_x, comment);
	comment = "position component y of galaxy [cMpc/h]";
	file.write_dataset("Galaxies/position_y", position_y, comment);
	comment = "position component z of galaxy [cMpc/h]";
	file.write_dataset("Galaxies/position_z", position_z, comment);

	//Galaxy velocity
	comment = "peculiar velocity component x of galaxy [km/s]";
	file.write_dataset("Galaxies/velocity_x", velocity_x, comment);
	comment = "peculiar velocity component y of galaxy [km/s]";
	file.write_dataset("Galaxies/velocity_y", velocity_y, comment);
	comment = "peculiar velocity component z of galaxy [km/s]";
	file.write_dataset("Galaxies/velocity_z", velocity_z, comment);

	//Galaxy type.
	comment = "galaxy type; =0 for centrals; =1 for satellites that reside in well identified subhalos; =2 for orphan satellites";
	file.write_dataset("Galaxies/type", type, comment);

	//Galaxy IDs.
	comment = "subhalo ID. Unique to this snapshot.";
	file.write_dataset("Galaxies/id_subhalo", id_subhalo, comment);
	comment = "halo ID. Unique to this snapshot.";
	file.write_dataset("Galaxies/id_halo", id_halo, comment);
	comment = "galayx ID. Unique to this snapshot.";
	file.write_dataset("Galaxies/id_galaxy", id_galaxy, comment);
}

void HDF5GalaxyWriter::write_global_properties (hdf5::Writer file, int snapshot, TotalBaryon &AllBaryons){

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
	file.write_dataset("Global/SFR_burst",AllBaryons.SFR_bulge), comment;

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

			float defl_value = -1;

			long gal_id = 1;

			for (auto &halo: halos){
				for (auto &subhalo: halo->all_subhalos()){
					for (auto &galaxy: subhalo->galaxies){

						vector<float> sfh_gal_disk;
						vector<float> star_gal_disk;
						vector<float> star_metals_gal_disk;
						vector<float> gas_gal_disk;
						vector<float> gas_metals_gal_disk;
						vector<float> sfh_gal_bulge;
						vector<float> star_gal_bulge;
						vector<float> star_metals_gal_bulge;
						vector<float> gas_gal_bulge;
						vector<float> gas_metals_gal_bulge;

						bool star_gal_bulge_exists = false;
						for(int s=sim_params.min_snapshot; s <= snapshot; s++) {

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
								}
								sfh_gal_disk.push_back(defl_value);
								star_gal_disk.push_back(defl_value);
								star_metals_gal_disk.push_back(defl_value);
								gas_gal_disk.push_back(defl_value);
								gas_metals_gal_disk.push_back(defl_value);

								sfh_gal_bulge.push_back(defl_value);
								star_gal_bulge.push_back(defl_value);
								star_metals_gal_bulge.push_back(defl_value);
								gas_gal_bulge.push_back(defl_value);
								gas_metals_gal_bulge.push_back(defl_value);
							}
							else {
								star_gal_bulge_exists = true;
								auto item = *it;
								// assign disk properties
								sfh_gal_disk.push_back(item.sfr_disk);
								star_gal_disk.push_back(item.stellar_disk.mass);
								star_metals_gal_disk.push_back(item.stellar_disk.mass_metals);
								gas_gal_disk.push_back(item.gas_disk.mass);
								gas_metals_gal_disk.push_back(item.gas_disk.mass_metals);

								// assign bulge properties
								sfh_gal_bulge.push_back(item.sfr_bulge);
								star_gal_bulge.push_back(item.stellar_bulge.mass);
								star_metals_gal_bulge.push_back(item.stellar_bulge.mass_metals);
								gas_gal_bulge.push_back(item.gas_bulge.mass);
								gas_metals_gal_bulge.push_back(item.gas_bulge.mass_metals);
							}
						}

						sfhs_disk.emplace_back(std::move(sfh_gal_disk));
						stellar_mass_disk.emplace_back(std::move(star_gal_disk));
						stellar_metals_disk.emplace_back(std::move(star_metals_gal_disk));
						gas_hs_disk.emplace_back(std::move(gas_gal_disk));
						gas_metals_hs_disk.emplace_back(std::move(gas_metals_gal_disk));

						sfhs_bulge.emplace_back(std::move(sfh_gal_bulge));
						stellar_mass_bulge.emplace_back(std::move(star_gal_bulge));
						stellar_metals_bulge.emplace_back(std::move(star_metals_gal_bulge));
						gas_hs_bulge.emplace_back(std::move(gas_gal_bulge));
						gas_metals_hs_bulge.emplace_back(std::move(gas_metals_gal_bulge));

					}
				}
			}

			vector<float> redshifts;

			for (int i=sim_params.min_snapshot+1; i <= snapshot; i++){
				redshifts.push_back(sim_params.redshifts[i]);
			}

			//Write header
			write_header(file_sfh, snapshot);

			//Write disk component history.
			comment = "Star formation history of stars formed that by this output time end up in the disk [Msun/Gyr/h]";
			file_sfh.write_dataset("Disks/StarFormationHistories", sfhs_disk, comment);

			comment = "History of stellar mass that by this output time ends up in the disk [Msun/h]";
			file_sfh.write_dataset("Disks/StellarMassHistories", stellar_mass_disk, comment);

			comment = "History of stellar mass in metals that by this output time ends up in the disk [Msun/h]";
			file_sfh.write_dataset("Disks/StellarMassMetalsHistories", stellar_metals_disk, comment);

			//Write bulge component history.
			comment = "Star formation history of stars formed that by this output time end up in the bulge [Msun/Gyr/h]";
			file_sfh.write_dataset("Bulges/StarFormationHistories", sfhs_bulge, comment);

			comment = "History of stellar mass that by this output time ends up in the bulge [Msun/h]";
			file_sfh.write_dataset("Bulges/StellarMassHistories", stellar_mass_bulge, comment);

			comment = "History of stellar mass in metals that by this output time ends up in the bulge [Msun/h]";
			file_sfh.write_dataset("Bulges/StellarMassMetalsHistories", stellar_metals_bulge, comment);

			comment = "Redshifts of the history outputs.";
			file_sfh.write_dataset("Redshifts", redshifts, comment);

		}

	}
}

void ASCIIGalaxyWriter::write(int snapshot, const std::vector<HaloPtr> &halos, TotalBaryon &AllBaryons)
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

