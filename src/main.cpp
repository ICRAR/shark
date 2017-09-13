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
#include <vector>

#include <boost/program_options.hpp>
#include <gsl/gsl_errno.h>

#include "recycling.h"
#include "config.h"
#include "components.h"
#include "cosmology.h"
#include "dark_matter_halos.h"
#include "evolve_halos.h"
#include "exceptions.h"
#include "galaxy_mergers.h"
#include "logging.h"
#include "numerical_constants.h"
#include "physical_model.h"
#include "simulation.h"
#include "execution.h"
#include "reionisation.h"
#include "agn_feedback.h"
#include "merger_tree_reader.h"
#include "tree_builder.h"

using namespace std;

namespace shark {

void show_help(const char *prog, const boost::program_options::options_description &desc) {
	cout << endl;
	cout << "SHArk: Semianalytic Halos Ark" << endl;
	cout << endl;
	cout << "Usage: " << prog << " [options] config-file" << endl;
	cout << endl;
	cout << desc << endl;
}

void setup_logging(int lvl) {

	namespace log = ::boost::log;
	namespace trivial = ::boost::log::trivial;

	trivial::severity_level sev_lvl = trivial::severity_level(lvl);
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

/**
 * Main SHArk routine.
 *
 * Here we load the relevant information and do the basic loops to solve galaxy formation
 */
int run(int argc, char **argv) {

	namespace po = boost::program_options;

	po::options_description visible_opts("SHArk options");
	visible_opts.add_options()
		("help,h",      "Show this help message")
		("version,V",   "Show version and exit")
		("verbose,v",   po::value<int>()->default_value(3), "Verbosity level. Higher is more verbose");

	po::positional_options_description pdesc;
	pdesc.add("config-file", 1);

	po::options_description all_opts;
	all_opts.add(visible_opts);
	all_opts.add_options()
		("config-file", po::value<string>(), "SHArk config file");

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
		return 0;
	}
	if (vm.count("config-file") == 0 ) {
		cerr << "Missing mandatory <config-file> option" << endl;
		return 1;
	}

	// Set up logging with indicated log level
	int verbosity = vm["verbose"].as<int>();
	verbosity = min(max(verbosity, 0), 5);
	verbosity = 5 - verbosity;
	setup_logging(verbosity);

	install_gsl_error_handler();

	/* We read the parameters that have been given as input by the user.*/
	string config_file = vm["config-file"].as<string>();

	/**
	 * We load all relevant parameters and implement all relevant physical processes needed by the physical model.
	 */
	AGNFeedbackParameters agn_params(config_file);
	CosmologicalParameters cosmo_parameters(config_file);
	DarkMatterHaloParameters dark_matter_halo_parameters(config_file);
	ExecutionParameters exec_params(config_file);
	GalaxyMergerParameters merger_parameters(config_file);
	GasCoolingParameters gas_cooling_params(config_file);
	RecyclingParameters recycling_parameters(config_file);
	ReionisationParameters reio_params(config_file);
	SimulationParameters sim_params(config_file);
	StellarFeedbackParameters stellar_feedback_params(config_file);
	StarFormationParameters star_formation_params(config_file);

	std::shared_ptr<Cosmology> cosmology = std::make_shared<Cosmology>(cosmo_parameters);
	std::shared_ptr<DarkMatterHalos> dark_matter_halos = std::make_shared<DarkMatterHalos>(dark_matter_halo_parameters, cosmology, sim_params);
	std::shared_ptr<AGNFeedback> agnfeedback = std::make_shared<AGNFeedback>(agn_params, cosmology);

	Simulation simulation{sim_params, cosmology};
	GasCooling gas_cooling{gas_cooling_params, reio_params, cosmology, agnfeedback, dark_matter_halos};
	StellarFeedback stellar_feedback{stellar_feedback_params};
	StarFormation star_formation{star_formation_params, cosmology};

	std::shared_ptr<BasicPhysicalModel> basic_physicalmodel = std::make_shared<BasicPhysicalModel>(exec_params.ode_solver_precision, gas_cooling, stellar_feedback, star_formation, recycling_parameters);

	GalaxyMergers galaxy_mergers{merger_parameters, dark_matter_halos,basic_physicalmodel};

	HaloBasedTreeBuilder tree_builder(exec_params);

	// Read the merger tree files.
	// Each merger tree will be a construction of halos and subhalos
	// with their growth history.
	auto halos = SURFSReader(sim_params.tree_files_prefix).read_halos(exec_params.simulation_batches, *dark_matter_halos, sim_params);
	auto merger_trees = tree_builder.build_trees(halos, sim_params);

	// Create the first generation of galaxies in the first halos apprearing.
	tree_builder.create_galaxies(halos, *cosmology, *dark_matter_halos);

	// The way we solve for galaxy formation is snapshot by snapshot. The loop is performed out to max snapshot-1, because we
	// calculate evolution in the time from the current to the next snapshot.
	// We first loop over snapshots, and for a fixed snapshot,
	// we loop over merger trees.
	// Each merger trees has a set of halos at a given snapshot,
	// which in turn contain galaxies.
	for(int snapshot=sim_params.min_snapshot; snapshot <= sim_params.max_snapshot-1; snapshot++) {

		LOG(info) << "Will evolve galaxies in snapshot " << snapshot;

		//Calculate the initial and final time of this snapshot.
		double ti = simulation.convert_snapshot_to_age(snapshot);
		double tf = simulation.convert_snapshot_to_age(snapshot+1);

		vector<HaloPtr> all_halos_this_snapshot;

		for(auto &tree: merger_trees) {
			/*here loop over the halos this merger tree has at this time.*/
			for(auto &halo: tree->halos[snapshot]) {

				//Append this halo to the list of halos of this snapshot

				all_halos_this_snapshot.insert(all_halos_this_snapshot.end(), halo);

				/*Check if there are any mergers in this snapshot*/

				/*First, determine which subhalos are disappearing in this snapshot and calculate dynamical friction timescale.*/
				galaxy_mergers.merging_subhalos(halo);

				/*Second, evaluate which galaxies are merging in this halo.*/
				galaxy_mergers.merging_galaxies(halo, sim_params.redshifts[snapshot], tf-ti);

				/*populate halos. This function should evolve the subhalos inside the halo.*/
				populate_halos(basic_physicalmodel, halo, snapshot,  sim_params.redshifts[snapshot], tf-ti);

				/*transfer galaxies from this halo->subhalos to the next snapshot's halo->subhalos*/
				transfer_galaxies_to_next_snapshot(halo);
			}
		}

		/*Here you could include the physics that allow halos to speak to each other. This could be useful e.g. during reionisation.*/
		//do_stuff_at_halo_level(all_halos_this_snapshot);

//		/*write snapshots only if the user wants outputs at this time.*/
		if(std::find(exec_params.output_snapshots.begin(), exec_params.output_snapshots.end(), snapshot) != exec_params.output_snapshots.end() )
		{
			//write_output(snapshot, all_halos_this_snapshot);
		}

		destroy_galaxies_this_snapshot(all_halos_this_snapshot);

	}

	return 0;
}

} // namespace shark

int main(int argc, char **argv) {
	try {
		return shark::run(argc, argv);
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
