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
 */

#include <cmath>
#include <memory>

#include "config.h"
#ifdef SHARK_OPENMP
#include <omp.h>
#endif // SHARK_OPENMP

#include "components.h"
#include "evolve_halos.h"
#include "logging.h"
#include "numerical_constants.h"

using namespace std;

namespace shark {

void transfer_galaxies_to_next_snapshot(const std::vector<HaloPtr> &halos, Cosmology cosmology, TotalBaryon &AllBaryons, int snapshot){

	/**
	 * This function transfer galaxies of the subhalos of this snapshot into the subhalos of the next snapshot, and baryon components from subhalo to subhalo.
	 */

	// First make sure central subhalos at this snapshot have only one central galaxy.
	for(auto &halo: halos){
		for(SubhaloPtr &subhalo: halo->all_subhalos()) {
			// Make sure all SFRs (in mass and metals) are set to 0 for the next snapshot
			for (GalaxyPtr & galaxy: subhalo->galaxies){
				galaxy->sfr_bulge  = 0;
				galaxy->sfr_z_disk = 0;
				galaxy->sfr_disk   = 0;
				galaxy->sfr_z_bulge= 0;
			}

			// Check if this is the last snapshot this subhalo is identified. If so, galaxies have already been transferred in galaxy_mergers.cpp
			if(subhalo->last_snapshot_identified == subhalo->snapshot){
				continue;
			}

			vector<float> mbaryon;
			if(subhalo->subhalo_type == Subhalo::SATELLITE){
				for (auto &galaxy: subhalo->galaxies){
					if(galaxy->galaxy_type == Galaxy::CENTRAL){
						std::ostringstream os;
						os << "Satellite subhalo " << subhalo << " has at least 1 central galaxy";
						throw invalid_argument(os.str());
					}
				}
			}
			else{
				int i = 0;
				for (auto &galaxy: subhalo->galaxies){
					if(galaxy->galaxy_type == Galaxy::CENTRAL){
						i++;
						mbaryon.push_back(galaxy->baryon_mass());
					}
				}
				if(i==0){
					if (subhalo->galaxy_count() > 0) {
						std::ostringstream os;
						os << "Central Subhalo " << subhalo << " has no central galaxy";
						throw invalid_data(os.str());
					}
				}
				if(i>1){
					std::ostringstream os;
					os << "Central Subhalo " << subhalo << " has " << i <<" central galaxies";
					os << "Baryon masses:" << mbaryon[0] << " " << mbaryon[1];
					throw invalid_argument(os.str());
				}
			}
		}
	}

	unsigned int subhalos_without_descendant = 0;
	double baryon_mass_loss = 0;

	for(auto &halo: halos){
		for(auto &subhalo: halo->all_subhalos()) {

			auto descendant_subhalo = subhalo->descendant;

			if (!descendant_subhalo) {
				subhalos_without_descendant++;
				baryon_mass_loss += subhalo->total_baryon_mass();
				continue;
			}

			// If the current subhalo is not a main progenitor, then make its central galaxy type 1. The other ones should already be type 1 or 2.
			if(!subhalo->main_progenitor or descendant_subhalo->subhalo_type == Subhalo::SATELLITE){
				auto central = subhalo->central_galaxy();
				if(central){
					central->galaxy_type = Galaxy::TYPE1;
				}
			}

			// Transfer galaxies to descendant subhalo.
			subhalo->transfer_galaxies_to(descendant_subhalo);

			if(subhalo->snapshot != descendant_subhalo->snapshot-1){
				std::ostringstream os;
				os << "Descendant subhalo is not in the subsequent snapshot";
				throw invalid_argument(os.str());
			}

			// Transfer subhalo baryon components.
			descendant_subhalo->cold_halo_gas += subhalo->cold_halo_gas;
			descendant_subhalo->hot_halo_gas += subhalo->hot_halo_gas;
			descendant_subhalo->ejected_galaxy_gas += subhalo->ejected_galaxy_gas;
			if (subhalo->main_progenitor) {
				descendant_subhalo->cooling_subhalo_tracking = subhalo->cooling_subhalo_tracking;
			}
		}
	}

	// Make sure there is only one central galaxy per halo/central in the descendant subhalos.
	for(auto &halo: halos){
		for(SubhaloPtr &subhalo: halo->all_subhalos()) {
			auto descendant_subhalo = subhalo->descendant;
			if(!descendant_subhalo){
				continue;
			}
			if(descendant_subhalo->subhalo_type == Subhalo::SATELLITE){
				for (auto &galaxy: descendant_subhalo->galaxies){
					if(galaxy->galaxy_type == Galaxy::CENTRAL){
						std::ostringstream os;
						os << "Satellite subhalo " << descendant_subhalo << " has at least 1 central galaxy";
						throw invalid_argument(os.str());
					}
				}
			}
			else{ // subhalo is central.
				int i = 0;
				for (auto &galaxy: descendant_subhalo->galaxies){
					if(galaxy->galaxy_type == Galaxy::CENTRAL){
						i++;
					}
				}
				if(i==0){
					if (descendant_subhalo->galaxy_count() > 0) {
						std::ostringstream os;
						os << "Central Subhalo " << descendant_subhalo << " has no central galaxy";
						throw invalid_data(os.str());
					}
				}
				if(i>1){
					std::ostringstream os;
					os << "Central Subhalo " << descendant_subhalo << " has " << i <<" central galaxies";
					throw invalid_argument(os.str());
				}
			}
		}
	}

	if (subhalos_without_descendant) {
		AllBaryons.baryon_total_lost[snapshot] = baryon_mass_loss;
		LOG(warning) << "Found " << subhalos_without_descendant << " subhalos without descendant while transferring galaxies.";
	}

}

void track_total_baryons(StarFormation &starformation, Cosmology &cosmology, ExecutionParameters execparams, SimulationParameters simulation_params, const std::vector<HaloPtr> &halos,
		TotalBaryon &AllBaryons, int snapshot, const molgas_per_galaxy &molgas, double deltat){


	BaryonBase mcold_total;
	BaryonBase mhothalo_total;
	BaryonBase mcoldhalo_total;
	BaryonBase mejectedhalo_total;
	BaryonBase mstars_total;
	BaryonBase mstars_bursts_galaxymergers;
	BaryonBase mstars_bursts_diskinstabilities;
	BaryonBase MBH_total;
	BaryonBase mHI_total;
	BaryonBase mH2_total;
	BaryonBase mDM_total;

	double SFR_total_disk = 0;
	double SFR_total_burst = 0;

	double total_baryons = 0;

	double z1 = simulation_params.redshifts[snapshot];
	double z2 = simulation_params.redshifts[snapshot+1];

	double mean_age = 0.5 * (cosmology.convert_redshift_to_age(z1) + cosmology.convert_redshift_to_age(z2));

	// Loop over all halos and subhalos to write galaxy properties
	for (auto &halo: halos){

		// accumulate dark matter mass
		mDM_total.mass += halo->Mvir;
        
		for (auto &subhalo: halo->all_subhalos()){
        
			// Accumulate subhalo baryons
			mhothalo_total.mass += subhalo->hot_halo_gas.mass;
			mhothalo_total.mass_metals += subhalo->hot_halo_gas.mass_metals;
        
			mcoldhalo_total.mass += subhalo->cold_halo_gas.mass;
			mcoldhalo_total.mass_metals += subhalo->cold_halo_gas.mass_metals;
        
			mejectedhalo_total.mass += subhalo->ejected_galaxy_gas.mass;
			mejectedhalo_total.mass_metals += subhalo->ejected_galaxy_gas.mass_metals;
        
			for (auto &galaxy: subhalo->galaxies){
        
				if(execparams.output_sf_histories){
        
					galaxy->mean_stellar_age += (galaxy->sfr_disk + galaxy->sfr_bulge) * deltat * mean_age;
					galaxy->total_stellar_mass_ever_formed += (galaxy->sfr_disk + galaxy->sfr_bulge) * deltat;
					HistoryItem hist_galaxy;
					hist_galaxy.sfr_disk    = galaxy->sfr_disk;
					hist_galaxy.sfr_bulge   = galaxy->sfr_bulge;
					hist_galaxy.sfr_z_disk  = galaxy->sfr_z_disk;
					hist_galaxy.sfr_z_bulge = galaxy->sfr_z_bulge;
					hist_galaxy.snapshot    = snapshot;
					galaxy->history.emplace_back(std::move(hist_galaxy));
				}
        
				//Accumulate galaxy baryons
				auto &molecular_gas = molgas.at(galaxy);
        
				mHI_total.mass += molecular_gas.m_atom + molecular_gas.m_atom_b;
				mH2_total.mass += molecular_gas.m_mol + molecular_gas.m_mol_b;
        
				mcold_total.mass += galaxy->disk_gas.mass + galaxy->bulge_gas.mass;
				mcold_total.mass_metals += galaxy->disk_gas.mass_metals + galaxy->bulge_gas.mass_metals;
        
				mstars_total.mass += galaxy->disk_stars.mass + galaxy->bulge_stars.mass;
				mstars_total.mass_metals += galaxy->disk_stars.mass_metals + galaxy->bulge_stars.mass_metals;
        
				mstars_bursts_galaxymergers.mass += galaxy->galaxymergers_burst_stars.mass;
				mstars_bursts_galaxymergers.mass_metals += galaxy->galaxymergers_burst_stars.mass_metals;
				mstars_bursts_diskinstabilities.mass += galaxy->diskinstabilities_burst_stars.mass;
				mstars_bursts_diskinstabilities.mass_metals += galaxy->diskinstabilities_burst_stars.mass_metals;

				SFR_total_disk  += galaxy->sfr_disk;
				SFR_total_burst += galaxy->sfr_bulge;
        
				MBH_total.mass += galaxy->smbh.mass;
        
			}
		}
	}

	total_baryons = mstars_total.mass + mcold_total.mass + MBH_total.mass + mhothalo_total.mass + mcoldhalo_total.mass + mejectedhalo_total.mass;

	AllBaryons.mstars.push_back(mstars_total);
	AllBaryons.mstars_burst_galaxymergers.push_back(mstars_bursts_galaxymergers);
	AllBaryons.mstars_burst_diskinstabilities.push_back(mstars_bursts_diskinstabilities);
	AllBaryons.mcold.push_back(mcold_total);
	AllBaryons.mHI.push_back(mHI_total);
	AllBaryons.mH2.push_back(mH2_total);
	AllBaryons.mBH.push_back(MBH_total);
	AllBaryons.SFR_disk.push_back(SFR_total_disk);
	AllBaryons.SFR_bulge.push_back(SFR_total_burst);

	AllBaryons.mhot_halo.push_back(mhothalo_total);
	AllBaryons.mcold_halo.push_back(mcoldhalo_total);
	AllBaryons.mejected_halo.push_back(mejectedhalo_total);

	AllBaryons.mDM.push_back(mDM_total);

	// Test for mass conservation.
	double all_bar = AllBaryons.baryon_total_created[snapshot] - AllBaryons.baryon_total_lost[snapshot];
	double frac = std::abs((total_baryons-all_bar)/all_bar);
			//std::abs(total_baryons / mDM_total.mass /  cosmology.universal_baryon_fraction() - 1.0);

	// If more than 1% of baryons have been lost, then stop the code.
	/*if(frac > 0.01){
		std::ostringstream os;
		os << "Accumulated baryon mass differs by " << frac << " with the total baryons that have been created by this snapshot.";
		throw invalid_data(os.str());
	}*/

}

}
