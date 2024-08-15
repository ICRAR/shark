//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2018
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

#include <algorithm>
#include <iomanip>
#include <iterator>
#include <memory>
#include <numeric>
#include <vector>

#include "components/algorithms.h"
#include "cosmology.h"
#include "dark_matter_halos.h"
#include "exceptions.h"
#include "halo.h"
#include "logging.h"
#include "merger_tree.h"
#include "omp_utils.h"
#include "ranges.h"
#include "subhalo.h"
#include "timer.h"
#include "total_baryon.h"
#include "tree_builder.h"

namespace shark {

TreeBuilder::TreeBuilder(ExecutionParameters exec_params, unsigned int threads) :
	exec_params(std::move(exec_params)), threads(threads)
{
	// no-op
}

TreeBuilder::~TreeBuilder() = default;

ExecutionParameters &TreeBuilder::get_exec_params()
{
	return exec_params;
}

void TreeBuilder::ensure_trees_are_self_contained(const std::vector<MergerTreePtr> &trees) const
{
	omp_static_for(trees, threads, [&](const MergerTreePtr &tree, unsigned int thread_idx) {
		for (auto &halo: tree->halos) {
			if (halo->merger_tree != tree) {
				std::ostringstream os;
				os << halo << " is not actually part of " << tree;
				throw invalid_data(os.str());
			}
		}
	});
}

void TreeBuilder::ignore_late_massive_halos(std::vector<MergerTreePtr> &trees, SimulationParameters sim_params, ExecutionParameters exec_params) 
{
	omp_static_for(trees, threads, [&](MergerTreePtr &tree, unsigned int thread_idx) {
		for (auto &root: tree->roots()) {
			if(root->Mvir > exec_params.ignore_npart_threshold * sim_params.particle_mass && sim_params.redshifts[root->snapshot] < exec_params.ignore_below_z){
				root->ignore_gal_formation = true;
			}
		}
	});

}


std::vector<MergerTreePtr> TreeBuilder::build_trees(std::vector<HaloPtr> &halos,
		SimulationParameters sim_params,
		GasCoolingParameters gas_cooling_params,
		DarkMatterHaloParameters dark_matter_params,
		const DarkMatterHalosPtr &darkmatterhalos,
		const CosmologyPtr &cosmology,
		TotalBaryon &AllBaryons)
{

	auto last_snapshot_to_consider = exec_params.last_output_snapshot();

	// Find roots and create Trees for each of them
	std::vector<MergerTreePtr> trees;
	MergerTree::id_t tree_counter = 0;
	for(const auto &halo: halos) {
		if (halo->snapshot == last_snapshot_to_consider) {
			auto tree = std::make_shared<MergerTree>(tree_counter++);
			if (LOG_ENABLED(debug)) {
				LOG(debug) << "Creating MergerTree at " << halo;
			}
			halo->merger_tree = tree;
			halo->merger_tree->add_halo(halo);
			trees.emplace_back(std::move(tree));
		}
	}

	// No halos found at desired snapshot, end now
	if (trees.empty()) {

		std::ostringstream os;
		os << "No Halo definitions found at snapshot " << last_snapshot_to_consider;
		os << ", cannot proceed any further with merger trees creation. " << std::endl;

		os << "Halos found at these snapshots: ";
		std::set<int> snapshots_found;
		for (const auto &halo: halos) {
			snapshots_found.insert(halo->snapshot);
		}
		std::copy(snapshots_found.begin(), snapshots_found.end(), std::ostream_iterator<int>(os, " "));
		os << std::endl;

		os << "Considering these snapshots during this run: ";
		auto &output_snaps = exec_params.output_snapshots;
		std::copy(output_snaps.begin(), output_snaps.end(), std::ostream_iterator<int>(os, " "));

		throw invalid_data(os.str());
	}


	loop_through_halos(halos, sim_params, exec_params);
	{
		Timer t;
		omp_static_for(trees, threads, [](MergerTreePtr &tree, unsigned int thread_idx) {
			tree->consolidate();
		});
		LOG(info) << "Took " << t << " to consolidate all trees";
	}

	// Ignore massive halos that pop in for the first time at low redshift. 
	// These generally are flaws of the merger tree builder.
	if(exec_params.ignore_late_massive_halos){
		ignore_late_massive_halos(trees, sim_params, exec_params);
	}

	// Make sure merger trees are fully self-contained
	ensure_trees_are_self_contained(trees);

	if(exec_params.ensure_mass_growth){
		// Ensure halos only grow in mass.
		LOG(info) << "Making sure halos only grow in mass";
		ensure_halo_mass_growth(trees, sim_params);
	}

	// Redefine angular momentum in the case of interpolated halos.
	// spin_interpolated_halos(trees, sim_params);

	// Define central subhalos
	LOG(info) << "Defining central subhalos";
	define_central_subhalos(trees, sim_params, dark_matter_params, darkmatterhalos);

	/*if(dark_matter_params.apply_fix_to_mass_swapping_events){
		// Function to redefine central subhalo properties from host halo
		LOG(info) << "Defining velocity, concentration and lambda of central subhalos";
		define_properties_central_subhalos(trees, sim_params, dark_matter_params, darkmatterhalos);
	}*/

 	// Define accretion rate from DM in case we want this.
	LOG(info) << "Defining accretion rate using cosmology";
	define_accretion_rate_from_dm(trees, sim_params, gas_cooling_params, *cosmology, AllBaryons);

	// Define halo and subhalos ages and other relevant properties
	LOG(info) << "Defining ages of halos and subhalos";
	define_ages_halos(trees, sim_params, darkmatterhalos);

	/*if(dark_matter_params.apply_fix_to_mass_swapping_events){
		// Function to redefine satellite subhalo properties with properties at infall
		LOG(info) << "Defining velocity, concentration and lambda of satellite subhalos";
		define_properties_satellite_subhalos(trees, sim_params, darkmatterhalos);
	}*/

	return trees;
}

void TreeBuilder::link(const SubhaloPtr &parent_shalo, const SubhaloPtr &desc_subhalo,
		       const HaloPtr &parent_halo, const HaloPtr &desc_halo) {

	// Establish ascendant and descendant links at subhalo level
	// Fail if subhalo has more than one descendant
	if (LOG_ENABLED(trace)) {
		LOG(trace) << "Connecting " << parent_shalo << " as a parent of " << desc_subhalo;
	}
	desc_subhalo->ascendants.push_back(parent_shalo);

	if (parent_shalo->descendant) {
		std::ostringstream os;
		os << parent_shalo << " already has a descendant " << parent_shalo->descendant;
		os << " but " << desc_subhalo << " is claiming to be its descendant as well";
		throw invalid_data(os.str());
	}
	parent_shalo->descendant = desc_subhalo;

	add_parent(desc_halo, parent_halo);
	
}

SubhaloPtr TreeBuilder::define_central_subhalo(HaloPtr &halo, SubhaloPtr &subhalo, SimulationParameters &sim_params,
		DarkMatterHaloParameters &dark_matter_params, const DarkMatterHalosPtr &darkmatterhalos)
{
	// point central subhalo to this subhalo.
	halo->central_subhalo = subhalo;
	halo->position = subhalo->position;
	halo->velocity = subhalo->velocity;

	double z = sim_params.redshifts[subhalo->snapshot];
	double npart = subhalo->Mvir/sim_params.particle_mass;

	// calculate accurate subhalo properties based on halo mass in the case of applying the fix to mass swapping events:
	if(dark_matter_params.apply_fix_to_mass_swapping_events){

	        auto mvir = halo->Mvir;

		subhalo->concentration = darkmatterhalos->nfw_concentration(mvir, z);
		
		if (subhalo->concentration < 1) {
		        throw invalid_argument("concentration is <1, cannot continue. Please check input catalogue");
		}
		subhalo->lambda = darkmatterhalos->halo_lambda(*subhalo, mvir, z, npart);
		subhalo->Vvir = darkmatterhalos->halo_virial_velocity(mvir, z);
	}

	halo->lambda = subhalo->lambda;
	// Calculate halos' vvir and concentration
	halo->Vvir = darkmatterhalos->halo_virial_velocity(halo->Mvir, z);
	halo->concentration = darkmatterhalos->nfw_concentration(halo->Mvir,z);

	/** If virial velocity of halo (which is calculated from the total mass
	and redshift) is smaller than the virial velocity of the central subhalo, which is
	directly calculated in VELOCIraptor, then adopt the VELOCIraptor one.**/
	if(halo->Vvir < subhalo->Vvir){
		halo->Vvir = subhalo->Vvir;
	}

	//remove subhalo from satellite list.
	remove_satellite(halo, subhalo);

	
	//define subhalo as central.
	subhalo->subhalo_type = Subhalo::CENTRAL;

	return subhalo;
}

void TreeBuilder::define_central_subhalos(const std::vector<MergerTreePtr> &trees, SimulationParameters &sim_params, DarkMatterHaloParameters &dark_matter_params,
					  const DarkMatterHalosPtr &darkmatterhalos){

	//This function loops over merger trees and halos to define central galaxies in a self-consistent way. The loop starts at z=0.

	//Loop over trees.
	omp_static_for(trees, threads, [&](const MergerTreePtr &tree, unsigned int thread_idx) {
		for (int snapshot=sim_params.max_snapshot; snapshot >= sim_params.min_snapshot; snapshot--) {

			for (auto &halo: tree->halos_at(snapshot)) {

				// First check if halo has a central subhalo, if yes, then continue with loop.
				if (halo->central_subhalo) {
					continue;
				}

				auto central_subhalo = halo->all_subhalos()[0];
				auto subhalo = define_central_subhalo(halo, central_subhalo, sim_params, dark_matter_params, darkmatterhalos);

				// save value of lambda to make sure that all main progenitors of this subhalo have the same lambda value. This is done for consistency 
				// throughout time.
				// NOTE: these are central subhalos, so use the host halo information
				auto z = sim_params.redshifts[subhalo->snapshot];
				auto npart = subhalo->Mvir/sim_params.particle_mass;
				auto lambda = darkmatterhalos->halo_lambda(*subhalo, halo->Mvir, z, npart);;

				// Now walk backwards through the main progenitor branch until subhalo has no more progenitors. This is done only in the case the ascendant
				// halo does not have a central already.

				// Loop going backwards through history:
				//  * Find the main progenitor of this subhalo and its host Halo.
				//  * Define that main progenitor as the central subhalo for the Halo (if none defined).
				//  * Define last_snapshot_identified for non-central ascendants.
				//  * Repeat
				auto ascendants = subhalo->ascendants;

				while (!ascendants.empty()) {

					// Check that there is a main progenitor first
					// If none is formally defined, we declare the most massive
					// ascendant to be the main progenitor
					auto main_prog = subhalo->main();
					if (!main_prog) {
						auto it = std::max_element(ascendants.begin(), ascendants.end(), [](const SubhaloPtr &s1, const SubhaloPtr &s2) {
							return s1->Mvir < s2->Mvir;
						});
						main_prog = *it;
						main_prog->main_progenitor = true;
						LOG(warning) << "No main progenitor defined for " << subhalo << ", defined "
									 << main_prog << " based on its Mvir";
					}

					auto ascendant_halo = main_prog->host_halo;

					// If a central subhalo has been defined, then its whole branch
					// has been processed, so there's no point on continuing.
					if (ascendant_halo->central_subhalo) {
						break;
					}

					// Redefine lambda of main progenitor to have the same one as its descendant, only if this halo is not reliable.
					if (!dark_matter_params.use_converged_lambda_catalog || (dark_matter_params.use_converged_lambda_catalog && main_prog->Mvir/sim_params.particle_mass < dark_matter_params.min_part_convergence)) {
						main_prog->lambda = lambda;

					}
					subhalo = define_central_subhalo(ascendant_halo, main_prog, sim_params, dark_matter_params, darkmatterhalos);

					// Define property last_identified_snapshot for all the ascendants that are not the main progenitor of the subhalo.
					for (auto &sub: ascendants) {
						if(!sub->main_progenitor){
							sub->last_snapshot_identified = sub->snapshot;
						}
					}
					// Now move to the ascendants of the main progenitor subhalo and repeat process.
					ascendants = subhalo->ascendants;

				}

			}

		}
	});

	// Make sure each halo has only one central subhalo and that the rest are satellites.
	omp_static_for(trees, threads, [&](const MergerTreePtr &tree, unsigned int thread_idx) {
		for (int snapshot=sim_params.min_snapshot; snapshot >= sim_params.max_snapshot; snapshot++) {

			for (auto &halo: tree->halos_at(snapshot)) {
				int i = 0;
				for (auto &subhalo: halo->all_subhalos()) {
					if(subhalo->subhalo_type == Subhalo::CENTRAL){
						i++;
						if (i > 1) {
							std::ostringstream os;
							os << "Halo " << halo << " has more than 1 central subhalo at snapshot " << snapshot;
							throw invalid_argument(os.str());
						}
					}
				}
				if (i == 0) {
					std::ostringstream os;
					os << "Halo " << halo << " has no central subhalo at snapshot " << snapshot;
					throw invalid_argument(os.str());
				}
			}
		}
	});
}

void TreeBuilder::ensure_halo_mass_growth(const std::vector<MergerTreePtr> &trees, SimulationParameters &sim_params){

	//This function loops over merger trees and halos to make sure that descendant halos are at least as massive as their progenitors.
	omp_static_for(trees, threads, [&](const MergerTreePtr &tree, unsigned int thread_idx) {
		/*for(int snapshot=sim_params.min_snapshot; snapshot < sim_params.max_snapshot; snapshot++) {

			for(auto &halo: tree->halos_at(snapshot)){
				// Check if current mass of halo is larger than descendant. If so, redefine descendant Mvir to that of the progenitor.
				if(halo->Mvir > halo->descendant->Mvir){
					halo->descendant->Mvir = halo->Mvir;
				}
			}
		}*/
		for(int snapshot=sim_params.max_snapshot; snapshot > sim_params.min_snapshot; snapshot--) {

			for(auto &halo: tree->halos_at(snapshot)){
				// Check if current mass of halo is smaller than the total of its ascendants. If so, redefine Mvir to that total.
				if(halo->Mvir < halo->total_mass_ascendants()){
					halo->Mvir = halo->total_mass_ascendants();
				}
			}
		}
	});
}

void TreeBuilder::spin_interpolated_halos(const std::vector<MergerTreePtr> &trees, SimulationParameters &sim_params){

	// This function loops over merger trees and halos and subhalos to make sure that interpolated halos have the right angular momentum and concentration information.
	// This has to be done starting from the first snapshot forward so that the angular momentum and concentration are propagated correctly if subhalo is interpolated over many snapshots.

	//Loop over trees.
	omp_static_for(trees, threads, [&](const MergerTreePtr &tree, unsigned int thread_idx) {
		for (int snapshot=sim_params.max_snapshot; snapshot >=sim_params.min_snapshot; snapshot--) {

			for (auto &halo: tree->halos_at(snapshot)) {

				for (auto &subhalo: halo->all_subhalos()) {
					//Check if subhalo is there because of interpolation. If so, redefine its angular momentum and concentration to that of its progenitor.
					if (subhalo->IsInterpolated) {
						auto main_progenitor = subhalo->main();
						subhalo->L = main_progenitor->L;
						subhalo->concentration = main_progenitor->concentration;
						subhalo->host_halo->concentration = main_progenitor->concentration;

						if (subhalo->concentration <= 0) {
							std::ostringstream os;
							os << "subhalo " << subhalo << " has concentration =0";
							throw invalid_argument(os.str());
						}
					}
				}
			}
		}
	});
}


void TreeBuilder::define_accretion_rate_from_dm(const std::vector<MergerTreePtr> &trees,
		SimulationParameters &sim_params,
		GasCoolingParameters &gas_cooling_params,
		Cosmology &cosmology,
		TotalBaryon &AllBaryons){


	//Loop over trees.
	auto universal_baryon_fraction = cosmology.universal_baryon_fraction();
	for(auto &tree: trees) {
		for(int snapshot=sim_params.max_snapshot; snapshot >= sim_params.min_snapshot; snapshot--) {
				for(auto &halo: tree->halos_at(snapshot)){

					auto Mvir_asc = halo->total_mass_ascendants();

					// Define accreted baryonic mass.
					halo->central_subhalo->accreted_mass = (halo->Mvir - Mvir_asc) * universal_baryon_fraction;

					// Avoid negative numbers
					if(halo->central_subhalo->accreted_mass < 0){
						halo->central_subhalo->accreted_mass = 0;
					}
				}
		}
	}

	// Now accummulate baryons staring from the highest redshift.
	double total_baryon_accreted = 0;

	for(int snapshot=sim_params.min_snapshot; snapshot <= sim_params.max_snapshot; snapshot++) {
		for(auto &tree: trees) {
				for(auto &halo: tree->halos_at(snapshot)){
					total_baryon_accreted += halo->central_subhalo->accreted_mass;
				}
		}
		// Keep track of the integral of the baryons mass accreted.
		AllBaryons.baryon_total_created[snapshot] = total_baryon_accreted;
	}

}

void TreeBuilder::define_ages_halos(const std::vector<MergerTreePtr> &trees,
		SimulationParameters &sim_params,
		const DarkMatterHalosPtr &darkmatterhalos){


	//Loop over trees.
	for(auto &tree: trees) {
		for(int snapshot=sim_params.max_snapshot; snapshot >= sim_params.min_snapshot; snapshot--) {
				for(auto &halo: tree->halos_at(snapshot)){

					auto prog = halo->main_progenitor();

					/*
					 * Define assembly ages of halos by going backwards in time and checking when the main progenitors had
					 * 50% and 80% of the mass of the current halo.
					 */
					auto snap = snapshot - 1;
					while(prog && (halo->age_50 == 0 || halo->age_80 == 0)){
						if(prog->Mvir <= 0.8* halo->Mvir && halo->age_80 == 0){
							halo->age_80 = sim_params.redshifts[snap];
						}
						if(prog->Mvir <= 0.5* halo->Mvir && halo->age_50 == 0){
							halo->age_50 = sim_params.redshifts[snap];
						}
						snap --;
						prog = prog->main_progenitor();
					}

					for (auto &subhalo: halo->satellite_subhalos) {

						auto main_prog = subhalo->main();
						auto snap = snapshot - 1;

						while(main_prog && subhalo->infall_t == 0){
							if(main_prog->subhalo_type == Subhalo::CENTRAL){

								subhalo->infall_t = sim_params.redshifts[snap];
								subhalo->Mvir_infall = main_prog->host_halo->Mvir;
								subhalo->rvir_infall = darkmatterhalos->halo_virial_radius(main_prog->host_halo->Mvir, sim_params.redshifts[snap]);

								subhalo->concentration_infall = darkmatterhalos->nfw_concentration(main_prog->host_halo->Mvir, sim_params.redshifts[snap]);
														
								if (subhalo->concentration_infall < 1) {
									throw invalid_argument("concentration is <1, cannot continue. Please check input catalogue");
								}

								subhalo->lambda_infall = darkmatterhalos->halo_lambda(*main_prog, main_prog->host_halo->Mvir,
															sim_params.redshifts[snap], main_prog->host_halo->Mvir/sim_params.particle_mass);
								subhalo->Vvir_infall = darkmatterhalos->halo_virial_velocity(main_prog->host_halo->Mvir, sim_params.redshifts[snap]);

								// properties directly taken from the catalogue: use central subhalo
								subhalo->Vcirc_infall = main_prog->host_halo->central_subhalo->Vcirc;
 								subhalo->L_infall.x = main_prog->host_halo->central_subhalo->L.x;
 								subhalo->L_infall.y = main_prog->host_halo->central_subhalo->L.y;
 								subhalo->L_infall.z = main_prog->host_halo->central_subhalo->L.z;

								//assume the stripping radius is equal to the virial radius at infall (which the largest it can be).
								subhalo->hot_halo_gas_r_rps = subhalo->rvir_infall;
							}
							snap --;
							main_prog = main_prog->main();
						}
					}
				}
		}

	}

}


void TreeBuilder::define_properties_central_subhalos(const std::vector<MergerTreePtr> &trees,
					SimulationParameters &sim_params,
					DarkMatterHaloParameters &dark_matter_params,
					const DarkMatterHalosPtr &darkmatterhalos){

		//Loop over trees
		for(auto &tree: trees) {
			for(int snapshot=sim_params.max_snapshot; snapshot >= sim_params.min_snapshot; snapshot--) {
				for(auto &halo: tree->halos_at(snapshot)){

					auto subhalo = halo->central_subhalo;
					// Calculate different properties: concentration, Vvir and lambda
					double mvir = halo->Mvir;
					double z= sim_params.redshifts[subhalo->snapshot];
					double npart = mvir/sim_params.particle_mass;

					subhalo->concentration = darkmatterhalos->nfw_concentration(mvir, z);

					if (subhalo->concentration < 1) {
						throw invalid_argument("concentration is <1, cannot continue. Please check input catalogue");
					}
					subhalo->lambda = darkmatterhalos->halo_lambda(*subhalo, mvir, z, npart);
					subhalo->Vvir = darkmatterhalos->halo_virial_velocity(mvir, z);

					//redefine halo properties based on new calculated central subhalo properties:
					halo->concentration = subhalo->concentration;
					halo->lambda = subhalo->lambda;
					halo->Vvir = subhalo->Vvir;
			}
		}
	}


	// now we will loop again to redefine lambda to the value of the main descendant in the cases of subhalos with unreliable values:
	for(auto &tree: trees) {
		for(int snapshot=sim_params.max_snapshot; snapshot >= sim_params.min_snapshot; snapshot--) {
			for(auto &halo: tree->halos_at(snapshot)){
				// First check if halo has a central subhalo, if yes, then continue with loop.
				if (!halo->central_subhalo) {
					continue;
				}
				auto main_prog = halo->central_subhalo->main();

				// Redefine lambda of main progenitor to have the same one as its descendant, only if this halo is not reliable.
				if (!dark_matter_params.use_converged_lambda_catalog || (dark_matter_params.use_converged_lambda_catalog && main_prog->Mvir/sim_params.particle_mass < dark_matter_params.min_part_convergence)) {
					main_prog->lambda = halo->central_subhalo->lambda;
					darkmatterhalos->redefine_angular_momentum(*main_prog, main_prog->lambda, sim_params.redshifts[main_prog->snapshot]);
				}
			}
		}
	}

}

void TreeBuilder::define_properties_satellite_subhalos(const std::vector<MergerTreePtr> &trees,
					  SimulationParameters &sim_params,
					  const DarkMatterHalosPtr &darkmatterhalos){

		//Loop over trees
		for(auto &tree: trees) {
			for(int snapshot=sim_params.max_snapshot; snapshot >= sim_params.min_snapshot; snapshot--) {
				for(auto &halo: tree->halos_at(snapshot)){
					for(auto &subhalo: halo->all_subhalos()){
						if(subhalo->subhalo_type == Subhalo::SATELLITE){

							double mvir = subhalo->Mvir;
							double z = sim_params.redshifts[subhalo->snapshot];

							// in the case of satellite subhalos with a well defined infall mass, we use the virial mass of its host at the time of infall, and the time of infall
							// to define other properties.
							if(subhalo->Mvir_infall > 0){
								mvir = subhalo->Mvir_infall;
								z = subhalo->infall_t;
							}

							double npart = mvir/sim_params.particle_mass;

							subhalo->concentration = darkmatterhalos->nfw_concentration(mvir, z);

							if (subhalo->concentration < 1) {
								throw invalid_argument("concentration is <1, cannot continue. Please check input catalogue");
							}

							subhalo->lambda = darkmatterhalos->halo_lambda(*subhalo, mvir, z, npart);
							subhalo->Vvir = darkmatterhalos->halo_virial_velocity(mvir, z);


						}
					}

				}
			}

		}
	
}

void TreeBuilder::remove_satellite(HaloPtr &halo, SubhaloPtr &subhalo){

	auto it = std::find(halo->satellite_subhalos.begin(), halo->satellite_subhalos.end(), subhalo);

	if (it == halo->satellite_subhalos.end()){
		std::ostringstream os;
		os << "Halo " << halo << " does not have satellite subhalos.";
		throw invalid_data(os.str());
	}

	halo->satellite_subhalos.erase(it);

}

HaloBasedTreeBuilder::HaloBasedTreeBuilder(ExecutionParameters exec_params, unsigned int threads) :
	TreeBuilder(std::move(exec_params), threads)
{
	// no-op
}

static std::vector<SubhaloPtr>::const_iterator find_by_id(const std::vector<SubhaloPtr> &subhalos, Subhalo::id_t id)
{
	return std::find_if(subhalos.begin(), subhalos.end(), [id](const SubhaloPtr &subhalo)
	{
		return subhalo->id == id;
	});
}

SubhaloPtr HaloBasedTreeBuilder::find_descendant_subhalo(
    const HaloPtr &halo, const SubhaloPtr &subhalo, const HaloPtr &descendant_halo)
{

	// if the descendant subhalo is not found in the descendant halos'
	// subhalos then we error
	auto descendant_subhalos = descendant_halo->all_subhalos();
	auto descendant_subhalo_found = find_by_id(descendant_subhalos, subhalo->descendant_id);

	if (descendant_subhalo_found == descendant_subhalos.end()) {

		std::ostringstream os;
		auto exec_params = get_exec_params();
		if (exec_params.skip_missing_descendants || exec_params.warn_on_missing_descendants) {
			os << "Descendant Subhalo id=" << subhalo->descendant_id;
			os << " for " << subhalo << " (mass: " << subhalo->Mvir << ") not found";
			os << " in the Subhalo's descendant Halo " << descendant_halo << std::endl;
			os << "Subhalos in " << descendant_halo << ": " << std::endl << "  ";
			auto all_subhalos = descendant_halo->all_subhalos();
			std::copy(all_subhalos.begin(), all_subhalos.end(),
					  std::ostream_iterator<SubhaloPtr>(os, "\n  "));
		}

		// Users can choose whether to continue in these situations
		// (with or without a warning) or if it should be considered an error
		if (!get_exec_params().skip_missing_descendants) {
			throw subhalo_not_found(os.str(), subhalo->descendant_id);
		}

		if (get_exec_params().warn_on_missing_descendants) {
			LOG(warning) << os.str();
		}
		halo->remove_subhalo(subhalo);
		return nullptr;
	}

	// We support only direct parentage; that is, descendants must be
	// in the snapshot directly after ours
	auto descendant_subhalo = *descendant_subhalo_found;
	if (subhalo->snapshot != descendant_subhalo->snapshot - 1) {
		std::ostringstream os;
		os << "Subhalo " << descendant_subhalo << " (snapshot " << subhalo->snapshot << ") ";
		os << "is not a direct descendant of " << subhalo << " (" << descendant_subhalo->snapshot << ").";
		throw invalid_data(os.str());
	}

	return descendant_subhalo;
}

static std::vector<HaloPtr>::iterator find_by_id(std::vector<HaloPtr> &halos, Halo::id_t id)
{

	auto lo = std::lower_bound(halos.begin(), halos.end(), id, [](const HaloPtr &x, Halo::id_t id)
	{
		return x->id < id;
	});
	auto up = std::upper_bound(halos.begin(), halos.end(), id, [](const Halo::id_t id, const HaloPtr &x)
	{
		return id < x->id;
	});
	if (lo == up) {
		return halos.end();
	}
	return lo;
}

void HaloBasedTreeBuilder::loop_through_halos(std::vector<HaloPtr> &halos, SimulationParameters sim_params, ExecutionParameters exec_params)
{
	sort_by_id(halos);

	// To find subhalos/halos that correspond to each other, we do the following
	//  1. Iterate over snapshots in descending order
	//  2. For each snapshot S we iterate over its halos
	//  3. For each halo H we iterate over its subhalos
	//  4. For each subhalo SH we find the halo with subhalo.descendant_halo_id
	//  5. When the descendant halo DH is found, we find the particular subhalo
	//     DSH inside DH that matches SH's descendant_id
	//  6. Now we have SH, H, DSH and DH. We link them all together,
	//     and to their tree

	// Get all snapshots in the Halos and sort them in decreasing order
	// (but skip the first one, those were already processed and MergerTrees
	// were built for them)
	std::vector<int> sorted_halo_snapshots;
	int min_particle_subhalo = 0;
	int particle_subhalo = 0;

	{
		std::set<int> halo_snapshots;
		for(const auto &halo: halos) {
			halo_snapshots.insert(halo->snapshot);
		        // compute the minimum particle number for structures
			for(const auto &subhalo: halo->all_subhalos()) {
			        if (exec_params.define_transient == ExecutionParameters::ZDEP_3SIGMA || exec_params.define_transient == ExecutionParameters::CONST_10MINPART){
				        // compute the particle number for each subhalo
			                particle_subhalo = round(subhalo->Mvir/sim_params.particle_mass);
					// find the minimum and save it
					if(min_particle_subhalo == 0 || particle_subhalo < min_particle_subhalo){
					        min_particle_subhalo = particle_subhalo;
					}
				}
			}
		}
		// minimum particle number for structures
		if (exec_params.define_transient == ExecutionParameters::ZDEP_3SIGMA || exec_params.define_transient == ExecutionParameters::CONST_10MINPART){
		        LOG(info) << "Minimum particle subhalo: " << min_particle_subhalo;
		}

		sorted_halo_snapshots = std::vector<int>(++(halo_snapshots.rbegin()), halo_snapshots.rend());
	}

	// Loop as per instructions above
	Timer t;
	// counters for massive transient statistics
	int count_transient_central = 0;
	int count_transient_sat = 0;
	int count_transient_mp = 0;
	// last snapshot
	int last_snapshot = sorted_halo_snapshots[0] + 1;
	for(int snapshot: sorted_halo_snapshots) {

	        
	        // flag massive transient subhalos if the transient fix is implemented
	        // these are subhalos that belong to main branches that are born with huge particle number
	        if(exec_params.apply_fix_to_massive_transient_events){
		        flag_massive_transients(halos, snapshot, last_snapshot, min_particle_subhalo, sim_params, exec_params);
		}
		
		LOG(info) << "Linking Halos/Subhalos at snapshot " << snapshot;

		int ignored = 0;
		auto halos_in_snapshot = make_range_filter(halos, in_snapshot(snapshot));
		for(auto &halo: halos_in_snapshot) {
			bool halo_linked = false;
			for(const auto &subhalo: halo->all_subhalos()) {


				// this subhalo has no descendants, let's not even try
				if (!subhalo->has_descendant) {
					if (LOG_ENABLED(debug)) {
						LOG(debug) << subhalo << " has no descendant, not following";
					}
					halo->remove_subhalo(subhalo);
					continue;
				}

				// if the descendant halo is not found, we don't consider this
				// halo anymore (and all its progenitors)
				auto descendant_halo_position = find_by_id(halos, subhalo->descendant_halo_id);
				if (descendant_halo_position == halos.end()) {
					if (LOG_ENABLED(debug)) {
						LOG(debug) << subhalo << " points to descendant halo/subhalo "
						           << subhalo->descendant_halo_id << " / " << subhalo->descendant_id
						           << ", which doesn't exist. Ignoring this halo and the rest of its progenitors";
					}
					ignored++;
					break;
				}
				auto descendant_halo = *descendant_halo_position;
				// hasn't been put in any merger tree, so it was ignored
				if (!descendant_halo->merger_tree) {
					continue;
				}

				auto descendant_subhalo = find_descendant_subhalo(halo, subhalo, descendant_halo);
				if (descendant_subhalo) {

				  // if fixes for the massive transients are applied, we look at the descendant halo
				  // if the descendant halo hosts any transient (previously flagged), we check whether this subhalo decreases in mass
				  // when there is a mass decrease, we evaluate if there is a more accurate connection with the transient
				  if(exec_params.apply_fix_to_massive_transient_events && descendant_halo->transient &&
				     subhalo->Mvir > descendant_subhalo->Mvir){
				          massive_transient_fix(subhalo, descendant_subhalo, halo, descendant_halo, exec_params,
								sorted_halo_snapshots[0]+1, count_transient_central, count_transient_sat, count_transient_mp);
				  }

				  // no fixes for massive transients, so normal subhalo-descendant subhalo connection is carried out
				  else{
				          link(subhalo, descendant_subhalo, halo, descendant_halo);
				  }
				  halo_linked = true;
			       	
				}
			}

			// If no subhalos were linked, this Halo will not have been linked,
			// meaning that it also needs to be ignored
			if (!halo_linked) {
				if (LOG_ENABLED(debug)) {
					LOG(debug) << halo << " doesn't contain any Subhalo pointing to"
					           << " descendants, ignoring it (and the rest of its progenitors)";
				}
				ignored++;
			}

		}

		if (LOG_ENABLED(debug)) {
			auto n_snapshot_halos = halos_in_snapshot.size();
			LOG(debug) << ignored << "/" << n_snapshot_halos << " ("
			           << std::setprecision(2) << std::setiosflags(std::ios::fixed)
			           << ignored * 100. / n_snapshot_halos << "%)"
			           << " Halos ignored at snapshot " << snapshot << " due to"
			           << " missing descendants (i.e., they were either the last Halo of"
			           << " their Halo family line, or they only hosted Subhalos"
			           << " that were the last Subhalo of their Subhalo families)";
		}

	}

	LOG(info) << "Linked all Halos/Subhalos in " << t;

	// statistics for massive transients
	if(exec_params.apply_fix_to_massive_transient_events){
	        LOG(info) << "Transient events: "
			  << count_transient_central << " fixed as central,  "
			  << count_transient_sat << " fixed as satellites,  "
			  << count_transient_mp << " removed main progenitor flag.";
	}

}


void HaloBasedTreeBuilder::flag_massive_transients(std::vector<HaloPtr> &halos, int snapshot, int last_snapshot, int min_particle_subhalo,
						   SimulationParameters sim_params, ExecutionParameters exec_params){

	// Flag those main branch subhalos born with an unexpected particle number:
        // 1. Flag all the subhalos that do not have main progenitor (subhalos that are born) as 'has_mainprogenitor = false'
        // 2. The ones that are part of a main branch are flagged as 'transient_candidate = true' 
	// 3. The candidates with Npart at birth above the defined threshold are flagged as 'transient = true'
	// (defined thresholds: constant Npart=200 for CONST_200, Npart=10*min(Npart,subhalo) for CONST_10MINPART
        // or Npart>3*sigma(Npart_subhalo at z) for ZDEP_3SIGMA)
  
	int particle_subhalo = 0;
	double three_sigma = 0;
	std::vector<int> Npart_subhalos;

	
	// To find subhalos/halos that do have main progenitor, we do the following
	//  1. For each snapshot S we iterate over its halos
	//  2. For each halo H we iterate over its subhalos
	//  3. For each subhalo SH we find the descendant: descendant is flagged as having main progenitor
	// if the subhalo is the main progenitor of that descendant
	//  4. Those that do not have main progenitor will be considered when flagging transient candidates 

	LOG(info) << "Finding if subhalos have main progenitor at snapshot " << snapshot+1;

	// loop over subhalos
	auto halos_in_snapshot = make_range_filter(halos, in_snapshot(snapshot));
	for(auto &halo: halos_in_snapshot) {
	        for(const auto &subhalo: halo->all_subhalos()) {

		        // this subhalo has no descendants, let's not even try
		        if (!subhalo->has_descendant) {
		                if (LOG_ENABLED(debug)) {
			                LOG(debug) << subhalo << " has no descendant, not following";
				}
				halo->remove_subhalo(subhalo);
				continue;
			}

			// if the descendant halo is not found, we don't consider this
			// halo anymore (and all its progenitors)
			auto descendant_halo_position = find_by_id(halos, subhalo->descendant_halo_id);
			if (descendant_halo_position == halos.end()) {
		                if (LOG_ENABLED(debug)) {
			                LOG(debug) << subhalo << " points to descendant halo/subhalo "
						   << subhalo->descendant_halo_id << " / " << subhalo->descendant_id
						   << ", which doesn't exist. Ignoring this halo and the rest of its progenitors";
				}
				break;
			}
			auto descendant_halo = *descendant_halo_position;
			// hasn't been put in any merger tree, so it was ignored
			if (!descendant_halo->merger_tree) {
		                continue;
			}
			
			auto descendant_subhalo = find_descendant_subhalo(halo, subhalo, descendant_halo);
			// if the descendant_subhalo is found, we check if the subhalo is its main progenitor
			if (descendant_subhalo) {
		                if(subhalo->main_progenitor){
			                // if that is the case, flag the descendant_subhalo as they have main progenitor
			                descendant_subhalo->has_mainprogenitor = true;
				}
			}
		}
	}


	// To find subhalos with no main progenitor that belong to main branches, we do the following
	//  1. For each snapshot S+1 we iterate over its halos
	//  2. For each halo H we iterate over its subhalos
	//  3. For each subhalo SH we select the ones with no main progenitor: we follow them forward
	// through main progenitors and if they reach the final snapshot, then they are main branches
	// (flag transient subhalos only if they belong to main branches)
	
	LOG(info) << "Defining Main Branch Subhalos born at snapshot " << snapshot+1;

	// loop over subhalos
	auto halos_in_snapshot0 = make_range_filter(halos, in_snapshot(snapshot+1));
	for(auto &halo: halos_in_snapshot0) {
	        for(const auto &subhalo: halo->all_subhalos()) {

			// if the subhalo does not have main progenitor (at birth)
			if(!subhalo->has_mainprogenitor){
			        // only define transients as subhalos that belong to main branches
			        // find the main branches following the subhalo forward in time through main progenitors
			        // define a main branch when:
			        // 1. the transient reaches the final snapshot through main progenitors
			        // 2. the descendant at the final snapshot is a central subhalo
	                        auto desc_mt = subhalo;
				while(desc_mt->descendant && desc_mt->main_progenitor){
				        desc_mt = desc_mt->descendant;
				}
			      
				// if tracking the progenitor branch we reach the a central subhalo at the final snapshot
				// then we have found out a main branch: we should flag it
				if(desc_mt->snapshot == last_snapshot && desc_mt->id == desc_mt->host_halo->id){

				        // flag it as main branch transient candidate
				        subhalo->transient_candidate = true;
		
					// calculate percentiles if we define transients with standard deviation
					// the particle number for each subhalo at birth is stored in an array for each snapshot
					// 3 times the standard deviation (99.85-th percentile) is computed for that array
					if(exec_params.define_transient == ExecutionParameters::ZDEP_3SIGMA){
					        particle_subhalo = round(subhalo->Mvir/sim_params.particle_mass);
						Npart_subhalos.push_back(particle_subhalo);
					}
				}
			}
		}
	}


	// calculate percentiles at each snapshot after the loop
	if(exec_params.define_transient == ExecutionParameters::ZDEP_3SIGMA){
	        // calculate percentiles if there are elements in the array (subhalos being born at this snap)
	        if(!Npart_subhalos.empty()){
		        three_sigma = percentiles(Npart_subhalos,0.9985); // value equivalent to 3 sigma
			LOG(info) << "Defining 3sigma ("<< three_sigma << ") for transient Subhalos at snapshot " << snapshot+1;
			// if the value is below 10 times the minimum particle number
			// 10 times the minimum particle number is considered
			if(three_sigma < 10*min_particle_subhalo){
			        three_sigma = 10*min_particle_subhalo;
			}
			Npart_subhalos.clear();
		}
		// if there are no subhalos born at this snapshots
		// 10 times the minimum particle number is considered	
		else{
		        three_sigma = 10*min_particle_subhalo;
		}
	}

		  
	// To define transient subhalos finally
	//  1. For each snapshot S+1 we iterate over its halos
	//  2. For each halo H we iterate over its subhalos
	//  3. For each subhalo SH we select the candidates: if its particle number is above
	// the imposed threshold, it is finally defined as transient

	LOG(info) << "Defining transient Subhalos at snapshot " << snapshot+1;
	
	// final loop over subhalos 
	auto halos_in_snapshotp1 = make_range_filter(halos, in_snapshot(snapshot+1));
	for(auto &halo: halos_in_snapshotp1) {
	        for(const auto &subhalo: halo->all_subhalos()) {
		        // to be transient the subhalo should be a candidate (main branch subhalo at birth)
		        if(subhalo->transient_candidate){
			        particle_subhalo = round(subhalo->Mvir/sim_params.particle_mass);
				// when that subhalo is born with a particle number that exceeds the imposed threshold
				// it is flagged as massive transient, and it should be studied in further detail later

				// constant threshold of 200 particles
				if((exec_params.define_transient == ExecutionParameters::CONST_200) && (particle_subhalo > 200)){
				        subhalo->transient = true;
					subhalo->host_halo->transient = true;
				}
				// constant threshold of 10 times the minimum particle number to form a subhalo
				else if((exec_params.define_transient == ExecutionParameters::CONST_10MINPART) && (particle_subhalo > 10*min_particle_subhalo)){
				        subhalo->transient = true;
					subhalo->host_halo->transient = true;
				}
				// redshift dependant threshold that is 3 times the standard deviation of the particle number
				// for main branches subhalo at birth at a particular snapshot
				else if((exec_params.define_transient == ExecutionParameters::ZDEP_3SIGMA) && (particle_subhalo > three_sigma)){
				        subhalo->transient = true;
					subhalo->host_halo->transient = true;
				}
			}
		}
		
	}
}





void TreeBuilder::massive_transient_fix(const SubhaloPtr &subhalo, const SubhaloPtr &descendant_subhalo,
					const HaloPtr &halo, const HaloPtr &descendant_halo, ExecutionParameters exec_params,
					int last_snapshot, int &count_transient_central, int &count_transient_sat, int &count_transient_mp){

        // This function gathers the fix for the massive transients 
        // in case the descendant halo for our subhalo hosts any transient, we check whether this subhalo decreases in mass
        // and in case there is a mass decrease, we evaluate if there is a more accurate connection involving the transient

        // flag to indicate if the transient has been linked
        // otherwise a normal link is carried out
        bool transient_linked = false;
					  
	// loop over all the descendant subhalos until finding the transient
	for(const auto &descendant_check_transient_subhalo: descendant_halo->all_subhalos()){
	        if(descendant_check_transient_subhalo->transient){
						    
		        // several conditions are imposed before breaking the subhalo-descendant subhalo connection
		        // and creating the subhalo-transient connection:
		        // 1. transient is more massive than current descendant
		        // 2. subhalo-descendant subhalo connection (which will be broken) has
		        // a decrease in mass greater than a threshold, controlled by "transient_lostmass_ratio" 
		        // 3. subhalo-transient connenction has a reasonable change in mass
		        // between a lower and an upper limit, controlled by "transient_gainedmass_ratio_low" and "transient_gainedmass_ratio_up"
		        if(descendant_check_transient_subhalo->Mvir > descendant_subhalo->Mvir &&
			   (descendant_subhalo->Mvir/subhalo->Mvir) < exec_params.transient_lostmass_ratio &&
			   (descendant_check_transient_subhalo->Mvir/subhalo->Mvir) < exec_params.transient_gainedmass_ratio_up &&
			   (descendant_check_transient_subhalo->Mvir/subhalo->Mvir) > exec_params.transient_gainedmass_ratio_low){
			  
//			        // looking at longer branches: compute how many snapshots forward in time the subhalos live
//			        int snap_sh = 0;
//				int snap_mt = 0;
//				auto desc_sh = descendant_subhalo;
//				while(desc_sh->descendant && desc_sh->main_progenitor){
//				        snap_sh++;
//					desc_sh = desc_sh->descendant;
//				}
//				auto desc_mt = descendant_check_transient_subhalo;
//				while(desc_mt->descendant && desc_mt->main_progenitor){
//				        snap_mt++;
//					desc_mt = desc_mt->descendant;
//				}
//				
//			        // if the transient length is longer then we should break the link
//				if(snap_mt > snap_sh){

			        // subhalo-transient connection
			        link(subhalo, descendant_check_transient_subhalo, halo, descendant_halo);
				subhalo->descendant_id = descendant_check_transient_subhalo->id;
				// not dealing with this particular transient anymore
				descendant_check_transient_subhalo->transient = false;
				// flag that the transient has been linked
				transient_linked = true;
					  
				if(descendant_check_transient_subhalo->id == descendant_halo->id){
				        // transient link counted as central for statistics
				        count_transient_central++;
				}
				else{
				        // transient link counted as satellite for statistics
				        count_transient_sat++;
				}
							    
				// if the mp flag had been removed and counted as that
				// we remove it from that category for statistics
				if(descendant_check_transient_subhalo->remove_mp_flag){
				        // redefine again the flag if it was removed
				        descendant_check_transient_subhalo->main_progenitor = true; 
					count_transient_mp--;
				}
				// transient fixed, loop can be skipped 
				break;
//			        }
			}

			// if the transient does not fulfill the previous criteria, at least
			// the mp flag is removed so that it does not deviate branches when
			// it is not the most massive halo
			else if(descendant_check_transient_subhalo->main_progenitor){
			        // remove main progenitor flag
			        descendant_check_transient_subhalo->main_progenitor = false; 
				descendant_check_transient_subhalo->remove_mp_flag = true;
				// count it as removing the mp flag
				count_transient_mp++;
			}
		}
	}
	
	// if transient is not linked, then apply normal subhalo-descendant subhalo connection
	if(!transient_linked){
	        link(subhalo, descendant_subhalo, halo, descendant_halo);
	}
}


double HaloBasedTreeBuilder::percentiles(std::vector<int> data, double val){

  // function that computes the val-th percentile of the array in data
  
  std::vector<double> ws(data.size(),1.); // array with weights=1
  std::vector<double> cumsum(data.size()); // array to get the x-axis

  // percentile should be between 0 and 1
  if((val<0)||(val>1)){
    std::ostringstream os;
    os << " Percentile value is not actually in the [0-1] range: val = " << val;
    throw invalid_data(os.str());    
  }  

  // sort the data in ascending order
  sort(data.begin(), data.end());

  double x0 = 0.;
  double x1 = 1.;
  double y0 = data[0];
  double y1 = data.back();

  // cumulative sum of the weights
  partial_sum(ws.begin(),ws.end(),cumsum.begin());

  for(std::size_t i = 0; i < data.size(); i++) {
    cumsum[i] = cumsum[i] - 1.*ws[1];
    if(data.size() != 0){
      cumsum[i] = cumsum[i]/(data.size()-1); // value in the x-axis
      // interpolation in the "val" value for the x-axis
      if(((val-cumsum[i])>0) && ((val-cumsum[i])<(val-x0))){
	x0 = cumsum[i];
	y0 = data[i];
      }
      if(((cumsum[i]-val)>0) && ((cumsum[i]-val)<(x1-val))){
	x1 = cumsum[i];
	y1 = data[i];
      }
    }
    else{
      std::ostringstream os;
      os << " STOP percentiles: problem with the function";
      throw invalid_data(os.str());        
    }	
  }
  // final interpolation
  return (y0*(x1-val)+y1*(val-x0))/(x1-x0);
}


}// namespace shark
