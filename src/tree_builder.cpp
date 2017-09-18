#include <iomanip>
#include <iterator>
#include <memory>
#include <numeric>
#include <vector>

#include "cosmology.h"
#include "exceptions.h"
#include "logging.h"
#include "tree_builder.h"


namespace shark {

TreeBuilder::TreeBuilder(ExecutionParameters exec_params) :
	exec_params(exec_params)
{
	// no-op
}

TreeBuilder::~TreeBuilder()
{
	// no-op
}

ExecutionParameters &TreeBuilder::get_exec_params()
{
	return exec_params;
}

std::vector<MergerTreePtr> TreeBuilder::build_trees(const std::vector<HaloPtr> &halos, SimulationParameters sim_params)
{

	const auto &output_snaps = exec_params.output_snapshots;
	auto last_snapshot_to_consider = *std::begin(output_snaps);

	// Find roots and create Trees for each of them
	std::vector<MergerTreePtr> trees;
	int tree_counter = 0;
	for(const auto &halo: halos) {
		if (halo->snapshot == last_snapshot_to_consider) {
			auto tree = std::make_shared<MergerTree>();
			tree->id = tree_counter++;
			LOG(debug) << "Creating MergerTree at " << halo;
			halo->merger_tree = tree;
			halo->merger_tree->add_halo(halo);
			trees.push_back(tree);
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
		std::copy(output_snaps.begin(), output_snaps.end(), std::ostream_iterator<int>(os, " "));

		throw invalid_data(os.str());
	}

	loop_through_halos(halos);

	// Define central galaxies
	define_central_subhalos(trees, sim_params);

	// Define accretion rate
	define_accretion_rate(trees, sim_params);

	return trees;
}

void TreeBuilder::link(const SubhaloPtr &subhalo, const SubhaloPtr &d_subhalo,
                       const HaloPtr &halo, const HaloPtr &d_halo) {

	// Establish ascendant and descendant links at subhalo level
	// Fail if subhalo has more than one descendant
	LOG(trace) << "Connecting " << subhalo << " as a parent of " << d_subhalo;
	d_subhalo->ascendants.push_back(subhalo);

	if (subhalo->descendant) {
		std::ostringstream os;
		os << subhalo << " already has a descendant " << subhalo->descendant;
		os << " but " << d_subhalo << " is claiming to be its descendant as well";
		throw invalid_data(os.str());
	}
	subhalo->descendant = d_subhalo;

	// Establish ascendant and descendant link at halo level
	// Fail if a halo has more than one descendant
	LOG(trace) << "Connecting " << halo << " as a parent of " << d_halo;
	d_halo->ascendants.push_back(halo);

	if (halo->descendant and halo->descendant->id != d_halo->id) {
		std::ostringstream os;
		os << halo << " already has a descendant " << halo->descendant;
		os << " but " << d_halo << " is claiming to be its descendant as well";
		throw invalid_data(os.str());
	}
	halo->descendant = d_halo;

	// Link this halo to merger tree and back
	if (!d_halo->merger_tree) {
		std::ostringstream os;
		os << "Descendant " << d_halo << " has no MergerTree associated to it";
		throw invalid_data(os.str());
	}
	halo->merger_tree = d_halo->merger_tree;
	halo->merger_tree->add_halo(halo);

}

void TreeBuilder::define_central_subhalos(std::vector<MergerTreePtr> trees, SimulationParameters sim_params){

	//This function loops over merger trees and halos to define central galaxies in a self-consistent way. The loop starts at z=0.

	//Loop over trees.
		for(auto &tree: trees) {

			for(int snapshot=sim_params.max_snapshot; snapshot >= sim_params.min_snapshot; snapshot--) {

				for(auto &halo: tree->halos[snapshot]){

					//First check in halo has a central subhalo, if yes, then continue with loop.
					if(halo->central_subhalo){
						continue;
					}

					//Only cases with no central subhalos.
					auto subhalos = halo->all_subhalos();

					// The first subhalo is always the most massive, so it will be the central.
					auto &subhalo = subhalos[0];

					// point central subhalo to this subhalo.
					halo->central_subhalo = subhalo;
					halo->position = subhalo->position;
					halo->velocity = subhalo->velocity;

					//remove subhalo from satellite list.
					remove_satellite(halo, subhalo);

					//define subhalo as central.
					subhalo->subhalo_type = Subhalo::CENTRAL;

					//Now look at most massive ascendants to define which one will be the central.
					auto ascendants = subhalo->ordered_ascendants();

					while(ascendants.size() > 0){
						//First ascendant is always the most massive.
						auto ascendant_central = ascendants[0];

						auto host_asc = ascendant_central->host_halo;

						host_asc->central_subhalo = ascendant_central;

						remove_satellite(host_asc, ascendant_central);

						ascendant_central->subhalo_type = Subhalo::CENTRAL;

						for (auto &sub: ascendants){
							// If subhalo is a satellite, then define this snapshot as the last identified one.
							if(sub->subhalo_type == Subhalo::SATELLITE){
								sub-> last_snapshot_identified = sub->snapshot;
							}
						}

						ascendants = ascendant_central->ordered_ascendants();

					}
				}
			}
		}


}

void TreeBuilder::define_accretion_rate(std::vector<MergerTreePtr> trees, SimulationParameters sim_params){

	//Loop over trees.
		for(auto &tree: trees) {

			for(int snapshot=sim_params.max_snapshot; snapshot >= sim_params.min_snapshot; snapshot--) {

				for(auto &halo: tree->halos[snapshot]){

					auto ascendants = halo->ascendants;

					halo->central_subhalo->accreted_mass = halo->Mvir;

					auto Mvir_asc = std::accumulate(ascendants.begin(), ascendants.end(), 0., [](double mass, const HaloPtr &halo) {
						return mass + halo->Mvir;
					});

					//Define accreted mass of dark matter.

					halo->central_subhalo->accreted_mass -= Mvir_asc;

					//Avoid negative numbers
					if(halo->central_subhalo->accreted_mass < 0){
						halo->central_subhalo->accreted_mass = 0;
					}

				}
			}
		}

}


void TreeBuilder::remove_satellite(HaloPtr halo, SubhaloPtr subhalo){

	auto it = std::find(halo->satellite_subhalos.begin(), halo->satellite_subhalos.end(), subhalo);

	if (it == halo->satellite_subhalos.end()){
		std::ostringstream os;
		os << "Halo " << halo << " does not have satellite subhalos.";
		throw invalid_data(os.str());
	}

	halo->satellite_subhalos.erase(it);

}

HaloBasedTreeBuilder::HaloBasedTreeBuilder(ExecutionParameters exec_params) :
	TreeBuilder(exec_params)
{
	// no-op
}

void HaloBasedTreeBuilder::loop_through_halos(const std::vector<HaloPtr> &halos)
{

	// Index all halos by snapshot and by ID, we'll need them later
	std::map<int, std::vector<HaloPtr>> halos_by_snapshot;
	std::map<Halo::id_t, HaloPtr> halos_by_id;
	for(const auto &halo: halos) {
		halos_by_snapshot[halo->snapshot].push_back(halo);
		halos_by_id[halo->id] = halo;
	}

	// To find subhalos/halos that correspond to each other, we do the following
	//  1. Iterate over snapshots in descending order
	//  2. For each snapshot S we iterate over its halos
	//  3. For each halo H we iterate over its subhalos
	//  4. For each subhalo SH we find the halo with subhalo.descendant_halo_id
	//     (which we globally keep at the halos_by_id map)
	//  5. When the descendant halo DH is found, we find the particular subhalo
	//     DSH inside DH that matches SH's descendant_id
	//  6. Now we have SH, H, DSH and DH. We link them all together,
	//     and to their tree

	// Get all snapshots in the Halos and sort them in decreasing order
	// (but skip the first one, those were already processed and MergerTrees
	// were built for them)
	std::set<int> halo_snapshots;
	for(const auto &halo: halos) {
		halo_snapshots.insert(halo->snapshot);
	}
	std::vector<int> sorted_halo_snapshots(++(halo_snapshots.rbegin()), halo_snapshots.rend());

	// Loop as per instructions above
	for(int snapshot: sorted_halo_snapshots) {

		LOG(info) << "Linking Halos/Subhalos at snapshot " << snapshot;

		int ignored = 0;
		for(const auto &halo: halos_by_snapshot[snapshot]) {

			bool halo_linked = false;
			for(const auto &subhalo: halo->all_subhalos()) {

				// this subhalo has no descendants, let's not even try
				if (!subhalo->has_descendant) {
					LOG(debug) << subhalo << " has no descendant, not following";
					continue;
				}

				// if the descendant halo is not found, we don't consider this
				// halo anymore (and all its progenitors)
				auto it = halos_by_id.find(subhalo->descendant_halo_id);
				if (it == halos_by_id.end()) {
					LOG(debug) << subhalo << " points to descendant halo/subhalo "
					           << subhalo->descendant_halo_id << " / " << subhalo->descendant_id
					           << ", which doesn't exist. Ignoring this halo and the rest of its progenitors";
					halos_by_id.erase(halo->id);
					ignored++;
					break;
				}

				// if the descendant subhalo is not found in the descendant halos'
				// subhalos then we error
				bool subhalo_descendant_found = false;
				const auto &d_halo = halos_by_id[subhalo->descendant_halo_id];
				for(auto &d_subhalo: d_halo->all_subhalos()) {
					if (d_subhalo->id == subhalo->descendant_id) {
						link(subhalo, d_subhalo, halo, d_halo);
						subhalo_descendant_found = true;
						halo_linked = true;
						break;
					}
				}
				if (!subhalo_descendant_found) {

					std::ostringstream os;
					os << "Descendant Subhalo id=" << subhalo->descendant_id;
					os << " for " << subhalo << " (mass: " << subhalo->Mvir << ") not found";
					os << " in the Subhalo's descendant Halo " << d_halo << std::endl;
					os << "Subhalos in " << d_halo << ": " << std::endl << "  ";
					auto all_subhalos = d_halo->all_subhalos();
					std::copy(all_subhalos.begin(), all_subhalos.end(),
					          std::ostream_iterator<SubhaloPtr>(os, "\n  "));

					// Users can choose whether to continue in these situations
					// (with a warning) or if it should be considered an error
					if (!get_exec_params().skip_missing_descendants) {
						throw subhalo_not_found(os.str(), subhalo->descendant_id);
					}

					LOG(warning) << os.str();
				}
			}

			// If no subhalos were linked, this Halo will not have been linked,
			// meaning that it also needs to be ignored
			if (!halo_linked) {
				LOG(debug) << halo << " doesn't contain any Subhalo pointing to"
				           << " descendants, ignoring it (and the rest of its progenitors)";
				halos_by_id.erase(halo->id);
				ignored++;
			}

		}

		auto n_snapshot_halos = halos_by_snapshot[snapshot].size();
		LOG(debug) << ignored << "/" << n_snapshot_halos << " ("
		          << std::setprecision(2) << std::setiosflags(std::ios::fixed)
		          << ignored * 100. / n_snapshot_halos << "%)"
		          << " Halos ignored at snapshot " << snapshot << " due to"
		          << " missing descendants (i.e., they were either the last Halo of"
		          << " their Halo family line, or they only hosted Subhalos"
		          << " that were the last Subhalo of their Subhalo families)";
	}

}


void HaloBasedTreeBuilder::create_galaxies(const std::vector<HaloPtr> &halos, Cosmology &cosmology, DarkMatterHalos &darkmatterhalos)
{

	//This function finds subhalos that have no progenitors (so first time they are identified) and are central, and creates a galaxy there.

	for (const auto &halo: halos){
		for(const auto &subhalo: halo->all_subhalos()) {
			if(subhalo->ascendants.empty() and subhalo->subhalo_type == Subhalo::CENTRAL and subhalo->galaxies.empty()){

				auto galaxy = std::make_shared<Galaxy>();

				galaxy->galaxy_type = Galaxy::CENTRAL;

				//assign an ad-hoc half-mass radius to start with.
				galaxy->disk_gas.rscale = darkmatterhalos.disk_size_theory(*subhalo);

				subhalo->galaxies.push_back(galaxy);

				subhalo->hot_halo_gas.mass = subhalo->host_halo->Mvir * cosmology.universal_baryon_fraction();


			}
		}
	}

}


}// namespace shark
