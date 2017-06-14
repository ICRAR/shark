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
#include "recycling.h"
#include "config.h"
#include "components.h"
#include "cosmology.h"
#include "evolve_halos.h"
#include "logging.h"
#include "numerical_constants.h"
#include "physical_model.h"
#include "simulation.h"
#include "execution.h"
#include "reionisation.h"

using namespace shark;
using namespace std;

// TODO: All these are stub, dummy methods that will go away one by one
//       into their specific parts of the code
vector<MergerTree> load_merger_trees() {
	return vector<MergerTree>();
}

void physical_processes() {
	return;
}

void do_stuff_at_halo_level(vector<shared_ptr<Halo>> halos) {
	return;
}

void write_output(int snapshot, const vector<MergerTree> &merger_trees) {
	return;
}

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

/**
 * Main SHArk routine.
 *
 * Here we load the relevant information and do the basic loops to solve galaxy formation
 */
int main(int argc, char **argv) {

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
	try {
		po::command_line_parser parser(argc, argv);
		parser.options(all_opts).positional(pdesc);
		po::store(parser.run(), vm);
	} catch (const boost::program_options::error &e) {
		cerr << "Error while parsing command-line: " << e.what() << endl;
		return 1;
	}
	notify(vm);

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

	/* We read the parameters that have been given as input by the user.*/
	string config_file = vm["config-file"].as<string>();

	ExecutionParameters exec_params(config_file);
	SimulationParameters sim_params(config_file);
	std::shared_ptr<Cosmology> cosmology = std::make_shared<Cosmology>(CosmologicalParameters(config_file));
	ReionisationParameters reio_params(config_file);
	Simulation simulation{sim_params, cosmology};
	GasCooling gas_cooling{GasCoolingParameters(config_file), reio_params};
	StellarFeedback stellar_feedback{StellarFeedbackParameters(config_file)};
	StarFormation star_formation{StarFormationParameters(config_file), cosmology};
	RecyclingParameters recycling_parameters;
	BasicPhysicalModel basic_physicalmodel(1e-6, gas_cooling, stellar_feedback, star_formation, recycling_parameters);

	// Read the merger tree files.
	// Each merger tree will be a construction of halos and subhalos
	// with their growth history.
	vector<MergerTree> merger_trees = load_merger_trees();

	// This function should return the system of differential equations
	// to be solved at each snapshot and in each galaxy.
	// This set of ODEs apply on the ideal case of a central galaxy, with no
	// AGN feedback. These equations are modified later if galaxies are
	// satellites or have an AGN.
	physical_processes();

	// The way we solve for galaxy formation is snapshot by snapshot. The loop is performed out to max snapshot-1, because we
	// calculate evolution in the time from the current to the next snapshot.
	// We first loop over snapshots, and for a fixed snapshot,
	// we loop over merger trees.
	// Each merger trees has a set of halos at a given snapshot,
	// which in turn contain galaxies.
	for(int snapshot=sim_params.min_snapshot; snapshot <= sim_params.max_snapshot-1; snapshot++) {
		//Calculate the initial and final time of this snapshot.
		double ti = simulation.convert_snapshot_to_age(snapshot);
		double tf = simulation.convert_snapshot_to_age(snapshot+1);

		for(MergerTree &tree: merger_trees) {
			/*here loop over the halos this merger tree has at this time.*/
			for(shared_ptr<Halo> halo: tree.halos[snapshot]) {
				/*populate halos. This function should evolve the subhalos inside the halo.*/
				populate_halos(basic_physicalmodel, halo, snapshot,  sim_params.redshifts[snapshot], tf-ti);
			}
		}

		/*Here you shoud include the physics that allow halos to speak to each other. This could be useful e.g. during reionisation.*/
		vector<shared_ptr<Halo>> all_halos_for_this_snapshot;
		do_stuff_at_halo_level(all_halos_for_this_snapshot);

//		/*write snapshots only if the user wants outputs at this time.*/
		if(std::find(exec_params.output_snapshots.begin(), exec_params.output_snapshots.end(), snapshot) != exec_params.output_snapshots.end() )
		{
			write_output(snapshot, merger_trees);
		}

	}

	return 0;
}
