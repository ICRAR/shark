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
 *
 * Main shark runner class
 */

#include <algorithm>
#include <cassert>
#include <functional>
#include <memory>
#include <numeric>
#include <ostream>
#include <vector>

#include "components/algorithms.h"
#include "evolve_halos.h"
#include "execution.h"
#include "disk_instability.h"
#include "environment.h"
#include "galaxy_creator.h"
#include "galaxy_mergers.h"
#include "galaxy_writer.h"
#include "logging.h"
#include "merger_tree.h"
#include "merger_tree_reader.h"
#include "omp_utils.h"
#include "options.h"
#include "physical_model.h"
#include "shark_runner.h"
#include "subhalo.h"
#include "timer.h"
#include "total_baryon.h"
#include "tree_builder.h"
#include "utils.h"

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

/// Structure containing detailed runtimes for the impl::evolve_merger_tree routine
struct evolution_times {
	Timer::duration galaxy_mergers = 0;
	Timer::duration subhalos_mergers = 0;
	Timer::duration galaxy_evolution = 0;
	Timer::duration disk_instability_evaluation = 0;

	evolution_times &operator +=(const evolution_times &rhs)
	{
		galaxy_mergers += rhs.galaxy_mergers;
		subhalos_mergers += rhs.subhalos_mergers;
		galaxy_evolution += rhs.galaxy_evolution;
		disk_instability_evaluation += rhs.disk_instability_evaluation;
		return *this;
	}

	evolution_times operator +(const evolution_times &rhs) const
	{
		evolution_times sum = *this;
		return sum += rhs;
	}
};

/// impl class definition
class SharkRunner::impl {
public:

	/// @see SharkRunner::SharkRunner(const Options &, unsigned int)
	impl(const Options &options, unsigned int threads) :
	    options(options), threads(threads),
	    total_evolution_times(threads),
	    cosmo_params(options), dark_matter_halo_params(options),
	    environment_params(options), exec_params(options),
	    gas_cooling_params(options),recycling_params(options), reincorporation_params(options),
	    simulation_params(options), star_formation_params(options),
	    cosmology(make_cosmology(cosmo_params)),
	    dark_matter_halos(make_dark_matter_halos(dark_matter_halo_params, cosmology, simulation_params, exec_params)),
	    writer(make_galaxy_writer(exec_params, cosmo_params, cosmology, dark_matter_halos, simulation_params, AGNFeedbackParameters(options))),
	    simulation(simulation_params, cosmology),
	    star_formation(star_formation_params, recycling_params, cosmology)
	{
		create_per_thread_objects();
	}

	/// @see SharkRunner::run
	void run();

	/// @see SharkRunner::report_total_times
	void report_total_times();

private:
	Options options;
	unsigned int threads;
	std::vector<evolution_times> total_evolution_times;
	CosmologicalParameters cosmo_params;
	DarkMatterHaloParameters dark_matter_halo_params;
	EnvironmentParameters environment_params;
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
	std::vector<PerThreadObjects> thread_objects;
	TotalBaryon all_baryons;
	Timer::duration evolution_time_total = 0;

	void create_per_thread_objects();
	std::vector<MergerTreePtr> import_trees();
	void log_snapshot_statistics(int snapshot, const std::vector<HaloPtr> &halos, const Timer &t) const;
	void evolve_merger_trees(const std::vector<std::vector<MergerTreePtr>> &all_trees, int snapshot);
	evolution_times evolve_merger_tree(const MergerTreePtr &tree, unsigned int thread_idx, int snapshot, double z, double delta_t);
	molgas_per_galaxy get_molecular_gas(const std::vector<HaloPtr> &halos, double z, bool calc_j);
	void add_to_total(const std::vector<evolution_times> &snapshot_evolution_times);
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

void SharkRunner::report_total_times()
{
	pimpl->report_total_times();
}

struct SnapshotStatistics {

	int snapshot;
	unsigned int threads;
	std::size_t starform_integration_intervals;
	std::size_t galaxy_ode_evaluations;
	std::size_t starburst_ode_evaluations;
	std::size_t n_halos;
	std::size_t n_subhalos;
	std::size_t n_galaxies;
	Timer::duration duration_millis;

	double galaxy_ode_evaluations_per_galaxy() const {
		if (n_galaxies == 0) {
			return 0;
		}
		return static_cast<double>(galaxy_ode_evaluations) / n_galaxies;
	}

	double evolution_rate() const
	{
		if (n_galaxies == 0) {
			return 0;
		}
		return double(duration_millis) / n_galaxies;
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
	   << "  Evolution rate, global:               " << fixed<3>(stats.evolution_rate()) << " [ms/gal]\n"
	   << "  Evolution rate, per thread:           " << fixed<3>(stats.evolution_rate() * stats.threads) << " [ms/gal]\n"
	   << "  Galaxy evolution ODE evaluations:     " << stats.galaxy_ode_evaluations
	   << " (" << fixed<3>(stats.galaxy_ode_evaluations_per_galaxy()) << " [evals/gal])" << "\n"
	   << "  Starburst ODE evaluations:            " << stats.starburst_ode_evaluations
	   << " (" << fixed<3>(stats.starburst_ode_evaluations_per_galaxy()) << " [evals/gal])" << "\n"
	   << "  Star formation integration intervals: " << stats.starform_integration_intervals
	   << " (" << fixed<3>(stats.starform_integration_intervals_per_galaxy_ode_evaluations()) << " [ints/eval])\n"
	   << "  Time:                                 " << fixed<3>(stats.duration_millis / 1000.) << " [s]";
	return os;
}

template <typename T>
std::basic_ostream<T> &operator<<(std::basic_ostream<T> &os, const evolution_times &times)
{
	os << "galaxy mergers: " << ns_time(times.galaxy_mergers)
	   << ", disk instability: " << ns_time(times.disk_instability_evaluation)
	   << ", galaxy evolution: " << ns_time(times.galaxy_evolution)
	   << ", subhalos mergers: " << ns_time(times.subhalos_mergers);
	return os;
}

void SharkRunner::impl::create_per_thread_objects()
{
	AGNFeedbackParameters agn_params(options);
	DiskInstabilityParameters disk_instability_params(options);
	EnvironmentParameters environment_params(options);
	GalaxyMergerParameters merger_parameters(options);
	ReionisationParameters reio_params(options);
	ReincorporationParameters reinc_params(options);
	StellarFeedbackParameters stellar_feedback_params(options);

	auto agnfeedback = make_agn_feedback(agn_params, cosmology, recycling_params, exec_params);
	auto environment = make_environment(environment_params, dark_matter_halos, cosmology, cosmo_params, simulation_params);
	auto reionisation = make_reionisation(reio_params);
	auto reincorporation = make_reincorporation(reinc_params, dark_matter_halos);
	StellarFeedback stellar_feedback {stellar_feedback_params};
	GasCooling gas_cooling {gas_cooling_params, star_formation_params, exec_params, reionisation, cosmology, agnfeedback, dark_matter_halos, reincorporation, environment};

	for(unsigned int i = 0; i != threads; i++) {
		auto physical_model = std::make_shared<BasicPhysicalModel>(exec_params.ode_solver_precision, gas_cooling, stellar_feedback, star_formation, *agnfeedback,
				recycling_params, gas_cooling_params, agn_params);
		GalaxyMergers galaxy_mergers(merger_parameters, cosmology, cosmo_params, exec_params, agn_params, simulation_params, dark_matter_halos, physical_model, agnfeedback);
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
	auto trees = tree_builder.build_trees(halos, simulation_params, gas_cooling_params, dark_matter_halo_params, dark_matter_halos, cosmology, all_baryons);
	LOG(info) << trees.size() << " Merger trees imported in " << t;

	// Create the first generation of galaxies if halo is first appearing
	LOG(info) << "Creating initial galaxies in central subhalos across all merger trees";
	GalaxyCreator galaxy_creator(cosmology, gas_cooling_params, simulation_params);
	galaxy_creator.create_galaxies(trees, all_baryons);

	return trees;
}

void _get_molecular_gas(const HaloPtr &halo, molgas_per_galaxy &molgas, StarFormation &star_formation, double z, bool calc_j)
{
	for (auto &subhalo: halo->all_subhalos()) {
		for (auto &galaxy: subhalo->galaxies) {
			molgas[galaxy.id] = star_formation.get_molecular_gas(galaxy, z, calc_j);
		}
	}
}

molgas_per_galaxy SharkRunner::impl::get_molecular_gas(const std::vector<HaloPtr> &halos, double z, bool calc_j)
{

	std::vector<StarFormation> star_formations(threads, star_formation);
	std::vector<molgas_per_galaxy> local_molgas(threads);

	omp_static_for(halos, threads, [&](const HaloPtr &halo, unsigned int idx){
		_get_molecular_gas(halo, local_molgas[idx], star_formations[idx], z, calc_j);
	});

	if (threads == 1) {
		return local_molgas[0];
	}

	return std::accumulate(local_molgas.begin(), local_molgas.end(), molgas_per_galaxy(), [](molgas_per_galaxy &mg1, molgas_per_galaxy &mg2) {
		mg1.insert(mg2.begin(), mg2.end());
		return mg1;
	});

}

void SharkRunner::impl::add_to_total(const std::vector<evolution_times> &times)
{
	transform(times.begin(), times.end(),
	          total_evolution_times.begin(), total_evolution_times.begin(),
	          std::plus<evolution_times>{}
	);
}

evolution_times SharkRunner::impl::evolve_merger_tree(const MergerTreePtr &tree, unsigned int thread_idx, int snapshot, double z, double delta_t)
{
	// Get the thread-specific objects needed to run the evolution
	// In the non-OpenMP case we simply have one
	auto &objs = thread_objects[thread_idx];
	auto &physical_model = objs.physical_model;
	auto &galaxy_mergers = objs.galaxy_mergers;
	auto &disk_instability = objs.disk_instability;

	evolution_times times;

	/*here loop over the halos this merger tree has at this time.*/
	for(auto &halo: tree->halos_at(snapshot)) {


		/*Evaluate which galaxies are merging in this halo.*/
		if (LOG_ENABLED(debug)) {
			LOG(debug) << "Merging galaxies in halo " << halo;
		}
		Timer t1;
		galaxy_mergers.merging_galaxies(halo, snapshot, delta_t);
		times.galaxy_mergers += t1.get();

		/*Evaluate disk instabilities.*/
		if (LOG_ENABLED(debug)) {
			LOG(debug) << "Evaluating disk instability in halo " << halo;
		}
		Timer t2;
		disk_instability.evaluate_disk_instability(halo, snapshot, delta_t);
		times.disk_instability_evaluation += t2.get();

		if (LOG_ENABLED(debug)) {
			LOG(debug) << "Evolving content in halo " << halo;
		}
		Timer t3;
		for(auto &subhalo: halo->all_subhalos()) {
			for(auto &galaxy: subhalo->galaxies) {
				physical_model->evolve_galaxy(*subhalo, galaxy, z, delta_t);
			}
		}
		times.galaxy_evolution += t3.get();

		/*Determine which subhalos are disappearing in this snapshot and calculate dynamical friction timescale and change galaxy types accordingly.*/
		if (LOG_ENABLED(debug)) {
			LOG(debug) << "Merging subhalos in halo " << halo;
		}
		Timer t4;
		galaxy_mergers.merging_subhalos(halo, z, snapshot);
		times.subhalos_mergers += t4.get();
	}

	return times;
}

void SharkRunner::impl::log_snapshot_statistics(int snapshot, const std::vector<HaloPtr> &halos, const Timer &t) const
{
	auto duration_millis = t.get() / 1000 / 1000;

	// Some high-level ODE and integration iteration count statistics
	auto starform_integration_intervals = std::accumulate(thread_objects.begin(), thread_objects.end(), std::size_t(0), [](std::size_t x, const PerThreadObjects &o) {
		return x + o.physical_model->get_star_formation_integration_intervals();
	});
	auto galaxy_ode_evaluations = std::accumulate(thread_objects.begin(), thread_objects.end(), std::size_t(0), [](std::size_t x, const PerThreadObjects &o) {
		return x + o.physical_model->get_galaxy_ode_evaluations();
	});
	auto starburst_ode_evaluations = std::accumulate(thread_objects.begin(), thread_objects.end(), std::size_t(0), [](std::size_t x, const PerThreadObjects &o) {
		return x + o.physical_model->get_galaxy_starburst_ode_evaluations();
	});
	auto n_halos = halos.size();
	auto n_subhalos = std::accumulate(halos.begin(), halos.end(), std::size_t(0), [](std::size_t n_subhalos, const HaloPtr &halo) {
		return n_subhalos + halo->subhalo_count();
	});
	auto n_galaxies = std::accumulate(halos.begin(), halos.end(), std::size_t(0), [](std::size_t n_galaxies, const HaloPtr &halo) {
		return n_galaxies + halo->galaxy_count();
	});

	SnapshotStatistics stats {snapshot, threads, starform_integration_intervals, galaxy_ode_evaluations, starburst_ode_evaluations,
							  n_halos, n_subhalos, n_galaxies, duration_millis};
	LOG(info) << "Statistics for snapshot " << snapshot << "\n" << stats;
}

std::vector<HaloPtr> all_halos_at_snapshot(const std::vector<std::vector<MergerTreePtr>> &all_trees, int snapshot)
{
	std::vector<HaloPtr> halos_at_snapshot;
	for (auto &merger_trees: all_trees) {
		for (auto &tree: merger_trees) {
			const auto &halos = tree->halos_at(snapshot);
			halos_at_snapshot.insert(halos_at_snapshot.end(), halos.begin(), halos.end());
		}
	}
	return halos_at_snapshot;
}

void SharkRunner::impl::evolve_merger_trees(const std::vector<std::vector<MergerTreePtr>> &all_trees, int snapshot)
{
	Timer snapshot_evolution_t;

	for(auto &o: thread_objects) {
		o.physical_model->reset_ode_evaluations();
	}

	// Calculate the initial and final time for the evolution start at this snapshot.
	auto z = simulation_params.redshifts[snapshot];
	auto z_end = simulation_params.redshifts[snapshot + 1];
	double ti = simulation.convert_snapshot_to_age(snapshot);
	double tf = simulation.convert_snapshot_to_age(snapshot + 1);
	auto delta_t = tf - ti;

	std::ostringstream os;
	os << "Will evolve galaxies from snapshot " << snapshot << " to " << snapshot + 1;
	os << ". Redshift: " << z << " -> " << z_end << ", time: " << ti << " -> " << tf;
	LOG(info) << os.str();

	Timer evolution_t;
	std::vector<evolution_times> times(threads);
	omp_static_for(all_trees, threads, [&](const std::vector<MergerTreePtr> &merger_trees, unsigned int thread_idx) {
		for (auto &tree: merger_trees) {
			times[thread_idx] += evolve_merger_tree(tree, thread_idx, snapshot, simulation_params.redshifts[snapshot], delta_t);
		}
	});
	auto evolution_duration = evolution_t.get();
	evolution_time_total += evolution_duration;
	LOG(info) << "Evolved galaxies in " << ns_time(evolution_duration);
	LOG(info) << "Detailed times: " << sum(times);
	add_to_total(times);

	// Collect this snapshot's halos across all merger trees
	auto all_halos_this_snapshot = all_halos_at_snapshot(all_trees, snapshot);

	bool write_galaxies = exec_params.output_snapshot(snapshot + 1);

	Timer molgas_t;
	auto molgas_per_gal = get_molecular_gas(all_halos_this_snapshot, z, write_galaxies);
	LOG(info) << "Calculated molecular gas in " << molgas_t;

	/*track all baryons of this snapshot*/
	Timer tracking_t;
	track_total_baryons(*cosmology, exec_params, simulation_params, all_halos_this_snapshot, all_baryons, snapshot, molgas_per_gal, delta_t);
	LOG(info) << "Total baryon amounts tracked in " << tracking_t;

	log_snapshot_statistics(snapshot, all_halos_this_snapshot, snapshot_evolution_t);

	/*transfer galaxies from this halo->subhalos to the next snapshot's halo->subhalos*/
	LOG(debug) << "Transferring all galaxies for snapshot " << snapshot << " into next snapshot";
	transfer_galaxies_to_next_snapshot(all_halos_this_snapshot, snapshot, all_baryons);

	// Collect next snapshot's halos across all merger trees
	// We keep them sorted so when output files are created the order in which
	// information appears is the same regardless of how many threads were used
	auto all_halos_next_snapshot = all_halos_at_snapshot(all_trees, snapshot + 1);
	sort_by_id(all_halos_next_snapshot);

	if (write_galaxies)
	{
		// Note that the output is being done at "snapshot + 1". This is because
		// we don't evolve galaxies AT snapshot "i" but FROM snapshot "i" TO
		// snapshot "i+1", and therefore at this point in time (after the actual
		// evolution) we consider our galaxies to be at snapshot "i+1"
		LOG(info) << "Write output files for evolution from snapshot " << snapshot << " to " << snapshot + 1;
		writer->write(snapshot + 1, all_halos_next_snapshot, all_baryons, molgas_per_gal);
	}

	/*reset instantaneous galaxy properties to 0 to initiate calculation at subsequent snapshot*/
	LOG(debug) << "Reseting all instantaneous galaxy properties to 0 at snapshot " << snapshot;
	reset_instantaneous_galaxy_properties(all_halos_next_snapshot, snapshot);

}

void SharkRunner::impl::report_total_times()
{
	for (unsigned int thread = 0; thread != threads; thread++) {
		LOG(info) << "Total evolution times in thread " << thread << ": "
				  << total_evolution_times[thread];
	}
	LOG(info) << "Total evolution times, combined: " << sum(total_evolution_times);
	LOG(info) << "Total evolution walltime: " << ns_time(evolution_time_total);
}

/// Produce similarly-weighted partitions of merger trees based on galaxy count
static std::vector<std::vector<MergerTreePtr>> partition_trees(std::vector<MergerTreePtr> &trees, unsigned int n_partitions)
{
	Timer partitioning_t;

	// Small structure to cache the result from tree->galaxy_count()
	struct tree_and_count {
		tree_and_count(MergerTreePtr &tree)
		    : tree(tree), galaxy_count(tree->galaxy_count())
		{
		}
		MergerTreePtr tree;
		size_t galaxy_count;
	};

	// Sort merger trees by galaxy count
	std::vector<tree_and_count> trees_and_counts(trees.begin(), trees.end());
	sort(trees_and_counts.begin(), trees_and_counts.end(), [](const tree_and_count &lhs, const tree_and_count &rhs) {
		return lhs.galaxy_count > rhs.galaxy_count;
	});

	// simple greedy partitioning
	std::vector<std::vector<MergerTreePtr>> partitions(n_partitions);
	std::vector<size_t> galaxy_counts(n_partitions);
	for (auto &tree_and_count: trees_and_counts) {
		auto distance = std::distance(galaxy_counts.begin(), std::min_element(galaxy_counts.begin(), galaxy_counts.end()));
		assert(distance >= 0);
		auto target = size_t(distance);
		partitions[target].push_back(tree_and_count.tree);
		galaxy_counts[target] += tree_and_count.galaxy_count;
	}

	LOG(info) << "Created tree partitions in " << partitioning_t;
	return partitions;
}

void SharkRunner::impl::run() {

	std::vector<MergerTreePtr> merger_trees = import_trees();
	auto tree_partitions = partition_trees(merger_trees, threads);

	// Go, go, go!
	// Note that we evolve galaxies in merger tress in the snapshot range [min, max)
	// This is because at snapshot "i" we don't evolve galaxies AT snapshot "i",
	// but rather FROM snapshot "i" TO snapshot "i+1".
	for(int snapshot = simulation_params.min_snapshot; snapshot <= simulation_params.max_snapshot - 1; snapshot++) {
		evolve_merger_trees(tree_partitions, snapshot);
	}

	report_total_times();
}

} // namespace shark
