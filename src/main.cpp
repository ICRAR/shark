//
// The main function for the shark executable
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

#include <algorithm>
#include <iostream>
#include <memory>
#include <unistd.h>
#include <vector>

#include "config.h"
#ifdef SHARK_OPENMP
#include <omp.h>
#endif // SHARK_OPENMP

#include <boost/filesystem/convenience.hpp>
#include <boost/program_options.hpp>
#include <gsl/gsl_errno.h>

#include "agn_feedback.h"
#include "components.h"
#include "cosmology.h"
#include "dark_matter_halos.h"
#include "disk_instability.h"
#include "evolve_halos.h"
#include "exceptions.h"
#include "galaxy_creator.h"
#include "galaxy_mergers.h"
#include "galaxy_writer.h"
#include "logging.h"
#include "numerical_constants.h"
#include "physical_model.h"
#include "simulation.h"
#include "execution.h"
#include "recycling.h"
#include "reincorporation.h"
#include "reionisation.h"
#include "merger_tree_reader.h"
#include "tree_builder.h"
#include "timer.h"
#include "utils.h"

using namespace std;

namespace shark {

void show_help(const char *prog, const boost::program_options::options_description &desc) {
	cout << endl;
	cout << "SHArk: Semianalytic Halos Ark" << endl;
	cout << endl;
	cout << "Usage: " << prog << " [options] config-file [... config-file]" << endl;
	cout << endl;
	cout << "Options are loaded from the given configuration files in order. If an option" << endl;
	cout << "is present in more than one configuration file, the last one takes precedence." << endl;
	cout << "Options specified via -o take precedence, in order." << endl;
	cout << endl;
	cout << desc << endl;
	cout << "Example:" << endl;
	cout << endl;
	cout << endl;
	cout << " $> " << prog << " -o group1.option1=a group1.option2=b config_file1.txt config_file2.txt" << endl;
	cout << endl;
	cout << " It loads options from config_file1.txt first and then from config_file2.txt. On top of that" << endl;;
	cout << " it also loads options 'group1.option1' and 'group1.option2' from the command-line." << endl;;
	cout << endl;
}

void setup_logging(int lvl) {

	namespace log = ::boost::log;
	namespace trivial = ::boost::log::trivial;

	trivial::severity_level sev_lvl = logging_level = trivial::severity_level(lvl);
	log::core::get()->set_filter([sev_lvl](log::attribute_value_set const &s) {
		return s["Severity"].extract<trivial::severity_level>() >= sev_lvl;
	});
}

void throw_exception_gsl_handler(const char *reason, const char *file, int line, int gsl_errno)
{
	throw gsl_error(reason, file, line, gsl_errno, gsl_strerror(gsl_errno));
}

void install_gsl_error_handler() {
	gsl_set_error_handler(&throw_exception_gsl_handler);
}

void log_startup_information(int argc, char **argv)
{
	char the_hostname[100];
	gethostname(the_hostname, 100);
	LOG(info) << "shark is starting in " << the_hostname;
	LOG(info) << "shark was built on " << __DATE__ << " " __TIME__;
	LOG(info) << "shark running at: " << boost::filesystem::current_path().string();

	std::ostringstream os;
	std::copy(argv, argv + argc, std::ostream_iterator<char *>(os, " "));
	LOG(info) << "shark started with command line: " << os.str();
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

//
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

std::vector<PerThreadObjects> create_per_thread_objects(
	std::size_t count, const Options &options,
	ExecutionParameters &exec_params, GasCoolingParameters &gas_cooling_params,
	RecyclingParameters &recycling_params, SimulationParameters &sim_params,
	StarFormationParameters &star_formation_params,
	std::shared_ptr<Cosmology> &cosmology, std::shared_ptr<DarkMatterHalos> &dark_matter_halos,
	StarFormation &star_formation)
{
	AGNFeedbackParameters agn_params(options);
	DiskInstabilityParameters disk_instability_params(options);
	GalaxyMergerParameters merger_parameters(options);
	ReionisationParameters reio_params(options);
	ReincorporationParameters reinc_params(options);
	StellarFeedbackParameters stellar_feedback_params(options);

	std::shared_ptr<AGNFeedback> agnfeedback = std::make_shared<AGNFeedback>(agn_params, cosmology);
	std::shared_ptr<Reionisation> reionisation = std::make_shared<Reionisation>(reio_params);
	std::shared_ptr<Reincorporation> reincorporation = std::make_shared<Reincorporation>(reinc_params, dark_matter_halos);
	StellarFeedback stellar_feedback {stellar_feedback_params};
	GasCooling gas_cooling {gas_cooling_params, star_formation_params, reionisation, cosmology, agnfeedback, dark_matter_halos, reincorporation};

	std::vector<PerThreadObjects> objects;
	for(std::size_t i = 0; i != count; i++) {
		auto physical_model = std::make_shared<BasicPhysicalModel>(exec_params.ode_solver_precision, gas_cooling, stellar_feedback, star_formation, recycling_params, gas_cooling_params);
		GalaxyMergers galaxy_mergers(merger_parameters, sim_params, dark_matter_halos, physical_model, agnfeedback);
		DiskInstability disk_instability(disk_instability_params, merger_parameters, sim_params, dark_matter_halos, physical_model, agnfeedback);
		objects.emplace_back(std::move(physical_model), std::move(galaxy_mergers), std::move(disk_instability));
	}
	return objects;
}

/**
 * Main SHArk routine.
 *
 * Here we load the relevant information and do the basic loops to solve galaxy formation
 */
int run(int argc, char **argv) {

	using std::string;
	using std::vector;
	namespace po = boost::program_options;

	po::options_description visible_opts("SHArk options");
	visible_opts.add_options()
		("help,h",      "Show this help message")
		("version,V",   "Show version and exit")
		("verbose,v",   po::value<int>()->default_value(3), "Verbosity level. Higher is more verbose")
#ifdef SHARK_OPENMP
		("threads,t",   po::value<std::size_t>()->default_value(1), "OpenMP threads, defaults to 1. 0 means use OpenMP default number of threads")
#endif // SHARK_OPENMP
		("options,o",   po::value<vector<string>>()->multitoken()->default_value({}, ""),
		                "Space-separated additional options to override config file");

	po::positional_options_description pdesc;
	pdesc.add("config-file", -1);

	po::options_description all_opts;
	all_opts.add(visible_opts);
	all_opts.add_options()
		("config-file", po::value<vector<string>>()->multitoken(), "SHArk config file(s)");

	// Read command-line options
	po::variables_map vm;
	po::command_line_parser parser(argc, argv);
	parser.options(all_opts).positional(pdesc);
	po::store(parser.run(), vm);
	po::notify(vm);

	if (vm.count("help")) {
		show_help(argv[0], visible_opts);
		return 0;
	}
	if (vm.count("version")) {
		cout << "SHArk version " << SHARK_VERSION << endl;
		cout << "OpenMP support: ";
#ifdef SHARK_OPENMP
		cout << "Yes" << endl;
#else
		cout << "No" << endl;
#endif
		return 0;
	}
	if (vm.count("config-file") == 0 ) {
		cerr << "At least one <config-file> option must be given" << endl;
		return 1;
	}

	// Set up logging with indicated log level
	int verbosity = vm["verbose"].as<int>();
	verbosity = min(max(verbosity, 0), 5);
	verbosity = 5 - verbosity;
	setup_logging(verbosity);

	log_startup_information(argc, argv);
#ifdef SHARK_OPENMP
auto threads = vm["threads"].as<std::size_t>();
if (threads == 0) {
	threads = omp_get_max_threads();
}
#else
std::size_t threads = 1;
#endif // SHARK_OPENMP
LOG(info) << "shark using " << threads << " thread(s)";

	install_gsl_error_handler();

	// Read the configuration file, and override options with any given
	// on the command-line
	Options options;
	for (auto &config_file: vm["config-file"].as<vector<string>>()) {
		options.add_file(config_file);
	}
	for(auto &opt_spec: vm["options"].as<vector<string>>()) {
		options.add(opt_spec);
	}

	/**
	 * We load all relevant parameters and implement all relevant physical processes needed by the physical model.
	 */
	CosmologicalParameters cosmo_parameters(options);
	DarkMatterHaloParameters dark_matter_halo_parameters(options);
	ExecutionParameters exec_params(options);
	GasCoolingParameters gas_cooling_params(options);
	RecyclingParameters recycling_parameters(options);
	ReincorporationParameters reinc_params(options);
	SimulationParameters sim_params(options);
	StarFormationParameters star_formation_params(options);

	std::shared_ptr<Cosmology> cosmology = std::make_shared<Cosmology>(cosmo_parameters);

	// TODO: Move this logic away from the main
	std::shared_ptr<DarkMatterHalos> dark_matter_halos;
	if (dark_matter_halo_parameters.haloprofile == DarkMatterHaloParameters::NFW) {
		dark_matter_halos = std::make_shared<NFWDarkMatterHalos>(dark_matter_halo_parameters, cosmology, sim_params);
	}
	else if (dark_matter_halo_parameters.haloprofile == DarkMatterHaloParameters::EINASTO) {
		dark_matter_halos = std::make_shared<EinastoDarkMatterHalos>(dark_matter_halo_parameters, cosmology, sim_params);
	}
	else {
		std::ostringstream os;
		os << "Dark Matter halo profile " << dark_matter_halo_parameters.haloprofile
		   << " not currently supported";
		throw invalid_argument(os.str());
	}

	Simulation simulation{sim_params, cosmology};
	StarFormation star_formation{star_formation_params, recycling_parameters, cosmology};

	auto per_thread_objects = create_per_thread_objects(threads, options, exec_params,
	    gas_cooling_params, recycling_parameters, sim_params, star_formation_params,
	    cosmology, dark_matter_halos, star_formation);

	// Inform the size of things
	std::ostringstream os;
	os << "Main structure/class sizes follow. ";
	os << "Baryon: " << memory_amount(sizeof(Baryon)) << ", Subhalo: " << memory_amount(sizeof(Subhalo)) << ", Halo: " << memory_amount(sizeof(Halo));
	os << ", Galaxy: " << memory_amount(sizeof(Galaxy)) << ", MergerTree: " << memory_amount(sizeof(MergerTree));
	LOG(info) << os.str();

	// Create class to track all the baryons of the simulation in its different components.
	auto AllBaryons = std::make_shared<TotalBaryon>();

	// Read the merger tree files.
	// Each merger tree will be a construction of halos and subhalos
	// with their growth history.
	std::vector<MergerTreePtr> merger_trees;
	{
		HaloBasedTreeBuilder tree_builder(exec_params, threads);
		auto halos = SURFSReader(sim_params.tree_files_prefix, threads).read_halos(exec_params.simulation_batches, *dark_matter_halos, sim_params);
		merger_trees = tree_builder.build_trees(halos, sim_params, cosmology, *AllBaryons);
		merger_trees.shrink_to_fit();
	}

	/* Create the first generation of galaxies if halo is first appearing.*/
	LOG(info) << "Creating initial galaxies in central subhalos across all merger trees";
	GalaxyCreator galaxy_creator(cosmology, gas_cooling_params, sim_params);
	galaxy_creator.create_galaxies(merger_trees, *AllBaryons);

	// TODO: move this logic away from the main
	// Also provide a std::make_unique
	std::unique_ptr<GalaxyWriter> writer;
	if (exec_params.output_format == Options::HDF5) {
		writer.reset(new HDF5GalaxyWriter(exec_params, cosmo_parameters, cosmology, dark_matter_halos, sim_params, star_formation));
	}
	else {
		writer.reset(new ASCIIGalaxyWriter(exec_params, cosmo_parameters, cosmology, dark_matter_halos, sim_params, star_formation));
	}

	// The way we solve for galaxy formation is snapshot by snapshot. The loop is performed out to max snapshot-1, because we
	// calculate evolution in the time from the current to the next snapshot.
	// We first loop over snapshots, and for a fixed snapshot,
	// we loop over merger trees.
	// Each merger trees has a set of halos at a given snapshot,
	// which in turn contain galaxies.
	for(int snapshot=sim_params.min_snapshot; snapshot <= sim_params.max_snapshot-1; snapshot++) {

		Timer t;
		LOG(info) << "Will evolve galaxies in snapshot " << snapshot << " corresponding to redshift "<< sim_params.redshifts[snapshot];

		for(auto &o: per_thread_objects) {
			o.physical_model->reset_ode_evaluations();
		}

		//Calculate the initial and final time of this snapshot.
		double ti = simulation.convert_snapshot_to_age(snapshot);
		double tf = simulation.convert_snapshot_to_age(snapshot+1);

		Timer evolution_t;
#ifdef SHARK_OPENMP
		#pragma omp parallel for num_threads(threads) schedule(static)
#endif // SHARK_OPENMP
		for (auto it = merger_trees.begin(); it < merger_trees.end(); it++) {

			const auto &tree = *it;

			// Get the thread-specific objects needed to run the evolution
			// In the non-OpenMP case we simply have one
#ifdef SHARK_OPENMP
			auto which = omp_get_thread_num();
#else
			auto which = 0;
#endif // SHARK_OPENMP
			auto &objs = per_thread_objects[which];
			auto &physical_model = objs.physical_model;
			auto &galaxy_mergers = objs.galaxy_mergers;
			auto &disk_instability = objs.disk_instability;

			/*here loop over the halos this merger tree has at this time.*/
			for(auto &halo: tree->halos[snapshot]) {

				/*Evaluate which galaxies are merging in this halo.*/
				LOG(debug) << "Merging galaxies in halo " << halo;
				galaxy_mergers.merging_galaxies(halo, snapshot, tf-ti);

				/*Evaluate disk instabilities.*/
				LOG(debug) << "Evaluating disk instability in halo " << halo;
				disk_instability.evaluate_disk_instability(halo, snapshot, tf-ti);

				/*populate halos. This function should evolve the subhalos inside the halo.*/
				LOG(debug) << "Evolving content in halo " << halo;
				populate_halos(physical_model, halo, snapshot,  sim_params.redshifts[snapshot], tf-ti);

				/*Determine which subhalos are disappearing in this snapshot and calculate dynamical friction timescale and change galaxy types accordingly.*/
				LOG(debug) << "Merging subhalos in halo " << halo;
				galaxy_mergers.merging_subhalos(halo, sim_params.redshifts[snapshot]);
			}
		}
		LOG(info) << "Evolved galaxies in " << evolution_t;

		vector<HaloPtr> all_halos_this_snapshot;
		for (auto &tree: merger_trees) {
			all_halos_this_snapshot.insert(all_halos_this_snapshot.end(), tree->halos[snapshot].begin(), tree->halos[snapshot].end());
		}

		/*track all baryons of this snapshot*/
		Timer tracking_t;
		track_total_baryons(star_formation, *cosmology, exec_params, all_halos_this_snapshot, *AllBaryons, sim_params.redshifts[snapshot], snapshot);
		LOG(info) << "Total baryon amounts tracked in " << tracking_t;

		/*Here you could include the physics that allow halos to speak to each other. This could be useful e.g. during reionisation.*/
		//do_stuff_at_halo_level(all_halos_this_snapshot);

		/*write snapshots only if the user wants outputs at this time (note that what matters here is snapshot+1).*/
		if(std::find(exec_params.output_snapshots.begin(), exec_params.output_snapshots.end(), snapshot+1) != exec_params.output_snapshots.end() )
		{
			LOG(info) << "Will write output file for snapshot " << snapshot+1;
			writer->write(snapshot, all_halos_this_snapshot, *AllBaryons);
		}

		auto duration_millis = t.get();

		// Some high-level ODE and integration iteration count statistics
		auto starform_integration_intervals = std::accumulate(per_thread_objects.begin(), per_thread_objects.end(), 0UL, [](unsigned long x, const PerThreadObjects &o) {
			return x + o.physical_model->get_star_formation_integration_intervals();
		});
		auto galaxy_ode_evaluations = std::accumulate(per_thread_objects.begin(), per_thread_objects.end(), 0UL, [](unsigned long x, const PerThreadObjects &o) {
			return x + o.physical_model->get_galaxy_ode_evaluations();
		});
		auto starburst_ode_evaluations = std::accumulate(per_thread_objects.begin(), per_thread_objects.end(), 0UL, [](unsigned long x, const PerThreadObjects &o) {
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
		transfer_galaxies_to_next_snapshot(all_halos_this_snapshot, *cosmology, *AllBaryons, snapshot);

	}

	return 0;
}

} // namespace shark

int main(int argc, char **argv) {
	try {
		shark::Timer timer;
		auto code = shark::run(argc, argv);
		if (code == 0) {
			LOG(info) << "Sucessfully finished in " << timer;
		}
		return code;
	} catch (const shark::missing_option &e) {
		std::cerr << "Missing option: " << e.what() << std::endl;
		return 1;
	} catch (const shark::exception &e) {
		std::cerr << "Unexpected shark exception found while running:" << std::endl << std::endl;
		std::cerr << e.what() << std::endl;
		return 1;
	} catch (const boost::program_options::error &e) {
		std::cerr << "Error while parsing command-line: " << e.what() << std::endl;
		return 1;
	} catch (const std::exception &e) {
		std::cerr << "Unexpected exception while running" << std::endl << std::endl;
		std::cerr << e.what() << std::endl;
		return 1;
	}
}
