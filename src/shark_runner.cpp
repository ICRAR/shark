//
// Main shark runner class
//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2018
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

#include <memory>
#include <numeric>
#include <ostream>
#include <vector>

#include "config.h"
#ifdef SHARK_OPENMP
#include <omp.h>
#endif // SHARK_OPENMP

#include "components.h"
#include "evolve_halos.h"
#include "execution.h"
#include "disk_instability.h"
#include "galaxy_creator.h"
#include "galaxy_mergers.h"
#include "galaxy_writer.h"
#include "logging.h"
#include "merger_tree_reader.h"
#include "options.h"
#include "physical_model.h"
#include "shark_runner.h"
#include "timer.h"
#include "tree_builder.h"

namespace shark {

// A simple structure that holds all objects that need to exist separately for each thread.
// Only those objects that maintain a state need to be here; the rest can be safely
// used across many threads.
// Most notoriously, the GSL ODE system and the GSL integrator (and their shark wrappers)
// maintain a state, so any class
// (which keeps a state, and therefore needs to be separate on each thread).
// The
struct PerThreadObjects
{
	PerThreadObjects(std::shared_ptr<BasicPhysicalModel> &&physical_model, GalaxyMergers &&galaxy_megers, DiskInstability &&disk_instability):
		physical_model(std::move(physical_model)), galaxy_mergers(std::move(galaxy_megers)), disk_instability(std::move(disk_instability)) {}
	std::shared_ptr<BasicPhysicalModel> physical_model;
	GalaxyMergers galaxy_mergers;
	DiskInstability disk_instability;
};

/// impl class definition
class SharkRunner::impl {
public:

	/// @see SharkRunner::SharkRunner(const Options &, unsigned int)
	impl(const Options &options, unsigned int threads) :
	    options(options), threads(threads),
	    cosmo_params(options), dark_matter_halo_params(options),
	    exec_params(options), gas_cooling_params(options),
	    recycling_params(options), reincorporation_params(options),
	    simulation_params(options), star_formation_params(options),
	    cosmology(make_cosmology(cosmo_params)),
	    dark_matter_halos(make_dark_matter_halos(dark_matter_halo_params, cosmology, simulation_params)),
	    writer(make_galaxy_writer(exec_params, cosmo_params, cosmology, dark_matter_halos, simulation_params)),
	    simulation(simulation_params, cosmology),
	    star_formation(star_formation_params, recycling_params, cosmology)
	{
		create_per_thread_objects();
	}

	/// @see SharkRunner::run
	void run();

private:
	Options options;
	unsigned int threads;
	CosmologicalParameters cosmo_params;
	DarkMatterHaloParameters dark_matter_halo_params;
	ExecutionParameters exec_params;
	GasCoolingParameters gas_cooling_params;
	RecyclingParameters recycling_params;
	ReincorporationParameters reincorporation_params;
	SimulationParameters simulation_params;
	StarFormationParameters star_formation_params;
	CosmologyPtr cosmology;
	DarkMatterHalosPtr dark_matter_halos;
	GalaxyWriterPtr writer;
	Simulation simulation;
	StarFormation star_formation;
	std::vector<PerThreadObjects> thread_objects {};
	TotalBaryon all_baryons {};

	void create_per_thread_objects();
	std::vector<MergerTreePtr> import_trees();
	void evolve_merger_trees(const std::vector<MergerTreePtr> &merger_trees, int snapshot);
	void evolve_merger_tree(const MergerTreePtr &tree, int snapshot, double z, double delta_t);
	molgas_per_galaxy get_molecular_gas(const std::vector<HaloPtr> &halos, double x, bool calc_j);
};

// Wiring pimpl to the original class
SharkRunner::SharkRunner(const Options &options, unsigned int threads) :
    pimpl(std::unique_ptr<impl>(new impl(options, threads)))
{
}

SharkRunner::~SharkRunner() = default;

void SharkRunner::run()
{
	pimpl->run();
}

struct SnapshotStatistics {

	int snapshot;
	unsigned long starform_integration_intervals;
	unsigned long galaxy_ode_evaluations;
	unsigned long starburst_ode_evaluations;
	unsigned long n_halos;
	unsigned long n_subhalos;
	unsigned long n_galaxies;
	unsigned long duration_millis;

	double galaxy_ode_evaluations_per_galaxy() const {
		if (n_galaxies == 0) {
			return 0;
		}
		return static_cast<double>(galaxy_ode_evaluations) / n_galaxies;
	}

	double starburst_ode_evaluations_per_galaxy() const {
		if (n_galaxies == 0) {
			return 0;
		}
		return static_cast<double>(starburst_ode_evaluations) / n_galaxies;
	}

	double starform_integration_intervals_per_galaxy_ode_evaluations() const {
		if (galaxy_ode_evaluations == 0) {
			return 0;
		}
		return static_cast<double>(starform_integration_intervals) / galaxy_ode_evaluations;
	}
};

template <typename T>
std::basic_ostream<T> &operator<<(std::basic_ostream<T> &os, const SnapshotStatistics &stats)
{
	os << "Snapshot " << stats.snapshot << "\n"
	   << "  Number of halos:                      " << stats.n_halos << "\n"
	   << "  Number of subhalos:                   " << stats.n_subhalos << "\n"
	   << "  Number of galaxies:                   " << stats.n_galaxies << "\n"
	   << "  Galaxy evolution ODE evaluations:     " << stats.galaxy_ode_evaluations
	   << " (" << fixed<3>(stats.galaxy_ode_evaluations_per_galaxy()) << " [evals/gal])" << "\n"
	   << "  Starburst ODE evaluations:            " << stats.starburst_ode_evaluations
	   << " (" << fixed<3>(stats.starburst_ode_evaluations_per_galaxy()) << " [evals/gal])" << "\n"
	   << "  Star formation integration intervals: " << stats.starform_integration_intervals
	   << " (" << fixed<3>(stats.starform_integration_intervals_per_galaxy_ode_evaluations()) << " [ints/eval])\n"
	   << "  Time:                                 " << fixed<3>(stats.duration_millis / 1000.) << " [s]";
	return os;
}


void SharkRunner::impl::create_per_thread_objects()
{
	AGNFeedbackParameters agn_params(options);
	DiskInstabilityParameters disk_instability_params(options);
	GalaxyMergerParameters merger_parameters(options);
	ReionisationParameters reio_params(options);
	ReincorporationParameters reinc_params(options);
	StellarFeedbackParameters stellar_feedback_params(options);

	auto agnfeedback = make_agn_feedback(agn_params, cosmology);
	auto reionisation = make_reionisation(reio_params);
	auto reincorporation = make_reincorporation(reinc_params, dark_matter_halos);
	StellarFeedback stellar_feedback {stellar_feedback_params};
	GasCooling gas_cooling {gas_cooling_params, star_formation_params, reionisation, cosmology, agnfeedback, dark_matter_halos, reincorporation};

	for(unsigned int i = 0; i != threads; i++) {
		auto physical_model = std::make_shared<BasicPhysicalModel>(exec_params.ode_solver_precision, gas_cooling, stellar_feedback, star_formation, recycling_params, gas_cooling_params);
		GalaxyMergers galaxy_mergers(merger_parameters, cosmology, simulation_params, dark_matter_halos, physical_model, agnfeedback);
		DiskInstability disk_instability(disk_instability_params, merger_parameters, simulation_params, dark_matter_halos, physical_model, agnfeedback);
		thread_objects.emplace_back(std::move(physical_model), std::move(galaxy_mergers), std::move(disk_instability));
	}
}

std::vector<MergerTreePtr> SharkRunner::impl::import_trees()
{
	Timer t;
	SURFSReader reader(simulation_params.tree_files_prefix, dark_matter_halos, simulation_params, threads);
	HaloBasedTreeBuilder tree_builder(exec_params, threads);
	auto halos = reader.read_halos(exec_params.simulation_batches);
	auto trees = tree_builder.build_trees(halos, simulation_params, cosmology, all_baryons);
	LOG(info) << "Merger trees imported in " << t;
	return trees;
}

void _get_molecular_gas(const HaloPtr &halo, molgas_per_galaxy &molgas, StarFormation &star_formation, double z, bool calc_j)
{
	for (auto &subhalo: halo->all_subhalos()) {
		for (auto &galaxy: subhalo->galaxies) {
			molgas[galaxy] = star_formation.get_molecular_gas(galaxy, z, calc_j);
		}
	}
}

molgas_per_galaxy SharkRunner::impl::get_molecular_gas(const std::vector<HaloPtr> &halos, double z, bool calc_j)
{

	std::vector<StarFormation> star_formations(threads, star_formation);
	std::vector<molgas_per_galaxy> local_molgas(threads);

#ifdef SHARK_OPENMP
	#pragma omp parallel for num_threads(threads) schedule(static)
#endif // SHARK_OPENMP
	for (auto it = halos.begin(); it < halos.end(); it++) {

#ifdef SHARK_OPENMP
		auto idx = omp_get_thread_num();
#else
		auto idx = 0;
#endif // SHARK_OPENMP

		_get_molecular_gas(*it, local_molgas[idx], star_formations[idx], z, calc_j);
	}

	if (threads == 1) {
		return local_molgas[0];
	}

	return std::accumulate(local_molgas.begin(), local_molgas.end(), molgas_per_galaxy(), [](molgas_per_galaxy &mg1, molgas_per_galaxy &mg2) {
		mg1.insert(mg2.begin(), mg2.end());
		return mg1;
	});

}

void SharkRunner::impl::evolve_merger_tree(const MergerTreePtr &tree, int snapshot, double z, double delta_t)
{
	// Get the thread-specific objects needed to run the evolution
	// In the non-OpenMP case we simply have one
#ifdef SHARK_OPENMP
	auto which = omp_get_thread_num();
#else
	auto which = 0;
#endif // SHARK_OPENMP
	auto &objs = thread_objects[which];
	auto &physical_model = objs.physical_model;
	auto &galaxy_mergers = objs.galaxy_mergers;
	auto &disk_instability = objs.disk_instability;

	/*here loop over the halos this merger tree has at this time.*/
	for(auto &halo: tree->halos[snapshot]) {

		/*Evaluate which galaxies are merging in this halo.*/
		if (LOG_ENABLED(debug)) {
			LOG(debug) << "Merging galaxies in halo " << halo;
		}
		galaxy_mergers.merging_galaxies(halo, snapshot, delta_t);

		/*Evaluate disk instabilities.*/
		if (LOG_ENABLED(debug)) {
			LOG(debug) << "Evaluating disk instability in halo " << halo;
		}
		disk_instability.evaluate_disk_instability(halo, snapshot, delta_t);

		if (LOG_ENABLED(debug)) {
			LOG(debug) << "Evolving content in halo " << halo;
		}
		for(auto &subhalo: halo->all_subhalos()) {
			for(auto &galaxy: subhalo->galaxies) {
				physical_model->evolve_galaxy(*subhalo, *galaxy, z, delta_t);
			}
		}

		/*Determine which subhalos are disappearing in this snapshot and calculate dynamical friction timescale and change galaxy types accordingly.*/
		if (LOG_ENABLED(debug)) {
			LOG(debug) << "Merging subhalos in halo " << halo;
		}
		galaxy_mergers.merging_subhalos(halo, z);
	}

}

void SharkRunner::impl::evolve_merger_trees(const std::vector<MergerTreePtr> &merger_trees, int snapshot)
{
	auto z = simulation_params.redshifts[snapshot];
	Timer t;
	LOG(info) << "Will evolve galaxies in snapshot " << snapshot << " corresponding to redshift "<< z;

	for(auto &o: thread_objects) {
		o.physical_model->reset_ode_evaluations();
	}

	//Calculate the initial and final time of this snapshot.
	double ti = simulation.convert_snapshot_to_age(snapshot);
	double tf = simulation.convert_snapshot_to_age(snapshot + 1);
	auto delta_t = tf - ti;

	Timer evolution_t;
#ifdef SHARK_OPENMP
	#pragma omp parallel for num_threads(threads) schedule(static)
#endif // SHARK_OPENMP
	for (auto it = merger_trees.begin(); it < merger_trees.end(); it++) {
		evolve_merger_tree(*it, snapshot, simulation_params.redshifts[snapshot], delta_t);
	}
	LOG(info) << "Evolved galaxies in " << evolution_t;

	std::vector<HaloPtr> all_halos_this_snapshot;
	for (auto &tree: merger_trees) {
		all_halos_this_snapshot.insert(all_halos_this_snapshot.end(), tree->halos[snapshot].begin(), tree->halos[snapshot].end());
	}

	bool write_galaxies = std::find(exec_params.output_snapshots.begin(), exec_params.output_snapshots.end(), snapshot+1) != exec_params.output_snapshots.end();

	Timer molgas_t;
	auto molgas_per_gal = get_molecular_gas(all_halos_this_snapshot, z, write_galaxies);
	LOG(info) << "Calculated molecular gas in " << molgas_t;

	/*track all baryons of this snapshot*/
	Timer tracking_t;
	track_total_baryons(star_formation, *cosmology, exec_params, simulation_params, all_halos_this_snapshot, all_baryons, snapshot, molgas_per_gal, delta_t);
	LOG(info) << "Total baryon amounts tracked in " << tracking_t;

	/*Here you could include the physics that allow halos to speak to each other. This could be useful e.g. during reionisation.*/
	//do_stuff_at_halo_level(all_halos_this_snapshot);

	/*write snapshots only if the user wants outputs at this time (note that what matters here is snapshot+1).*/
	if (write_galaxies)
	{
		LOG(info) << "Will write output file for snapshot " << snapshot+1;
		writer->write(snapshot, all_halos_this_snapshot, all_baryons, molgas_per_gal);
	}

	auto duration_millis = t.get();

	// Some high-level ODE and integration iteration count statistics
	auto starform_integration_intervals = std::accumulate(thread_objects.begin(), thread_objects.end(), 0UL, [](unsigned long x, const PerThreadObjects &o) {
		return x + o.physical_model->get_star_formation_integration_intervals();
	});
	auto galaxy_ode_evaluations = std::accumulate(thread_objects.begin(), thread_objects.end(), 0UL, [](unsigned long x, const PerThreadObjects &o) {
		return x + o.physical_model->get_galaxy_ode_evaluations();
	});
	auto starburst_ode_evaluations = std::accumulate(thread_objects.begin(), thread_objects.end(), 0UL, [](unsigned long x, const PerThreadObjects &o) {
		return x + o.physical_model->get_galaxy_starburst_ode_evaluations();
	});
	auto n_halos = all_halos_this_snapshot.size();
	auto n_subhalos = std::accumulate(all_halos_this_snapshot.begin(), all_halos_this_snapshot.end(), 0UL, [](unsigned long n_subhalos, const HaloPtr &halo) {
		return n_subhalos + halo->subhalo_count();
	});
	auto n_galaxies = std::accumulate(all_halos_this_snapshot.begin(), all_halos_this_snapshot.end(), 0UL, [](unsigned long n_galaxies, const HaloPtr &halo) {
		return n_galaxies + halo->galaxy_count();
	});

	SnapshotStatistics stats {snapshot, starform_integration_intervals, galaxy_ode_evaluations, starburst_ode_evaluations,
							  n_halos, n_subhalos, n_galaxies, duration_millis};
	LOG(info) << "Statistics for snapshot " << snapshot << std::endl << stats;


	/*transfer galaxies from this halo->subhalos to the next snapshot's halo->subhalos*/
	LOG(debug) << "Transferring all galaxies for snapshot " << snapshot << " into next snapshot";
	transfer_galaxies_to_next_snapshot(all_halos_this_snapshot, *cosmology, all_baryons, snapshot);

}

void SharkRunner::impl::run() {

	std::vector<MergerTreePtr> merger_trees = import_trees();

	/* Create the first generation of galaxies if halo is first appearing.*/
	LOG(info) << "Creating initial galaxies in central subhalos across all merger trees";
	GalaxyCreator galaxy_creator(cosmology, gas_cooling_params, simulation_params);
	galaxy_creator.create_galaxies(merger_trees, all_baryons);

	// Go, go, go!
	for(int snapshot = simulation_params.min_snapshot; snapshot <= simulation_params.max_snapshot-1; snapshot++) {
		evolve_merger_trees(merger_trees, snapshot);
	}
}

} // namespace shark
