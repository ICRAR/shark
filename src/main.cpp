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
 * The main function for the shark executable
 */

#include <algorithm>
#include <ios>
#include <iostream>
#include <ostream>
#include <vector>

#include "config.h"
#ifdef SHARK_OPENMP
#include <omp.h>
#endif // SHARK_OPENMP

#include <boost/filesystem/convenience.hpp>
#include <boost/program_options.hpp>
#include <gsl/gsl_errno.h>

#include "galaxy.h"
#include "merger_tree.h"
#include "logging.h"
#include "options.h"
#include "shark_runner.h"
#include "subhalo.h"
#include "git_revision.h"
#include "timer.h"

namespace shark {

void show_help(const char *prog, const boost::program_options::options_description &desc, std::ostream &out) {
	using std::endl;
	out << endl;
	out << "SHArk: Semianalytic Halos Ark" << endl;
	out << endl;
	out << "Usage: " << prog << " [options] config-file [... config-file]" << endl;
	out << endl;
	out << "Options are loaded from the given configuration files in order. If an option" << endl;
	out << "is present in more than one configuration file, the last one takes precedence." << endl;
	out << "Options specified via -o take precedence, in order." << endl;
	out << endl;
	out << desc << endl;
	out << "Example:" << endl;
	out << endl;
	out << endl;
	out << " $> " << prog << " -o group1.option1=a group1.option2=b config_file1.txt config_file2.txt" << endl;
	out << endl;
	out << " It loads options from config_file1.txt first and then from config_file2.txt. On top of that" << endl;;
	out << " it also loads options 'group1.option1' and 'group1.option2' from the command-line." << endl;;
	out << endl;
}

template <typename T>
void report_size_of_things(std::basic_ostream<T> &os)
{
	os << "Main structure/class sizes follow. ";
	os << "Baryon: " << memory_amount(sizeof(Baryon)) << ", Subhalo: " << memory_amount(sizeof(Subhalo)) << ", Halo: " << memory_amount(sizeof(Halo));
	os << ", Galaxy: " << memory_amount(sizeof(Galaxy)) << ", MergerTree: " << memory_amount(sizeof(MergerTree));
}

static
void setup_logging(const boost::program_options::variables_map &vm) {

	namespace log = ::boost::log;
	namespace trivial = ::boost::log::trivial;

	// Set up logging with indicated log level
	int verbosity = vm["verbose"].as<int>();
	verbosity = std::min(std::max(verbosity, 0), 5);
	verbosity = 5 - verbosity;

	trivial::severity_level sev_lvl = logging_level = trivial::severity_level(verbosity);
	log::core::get()->set_filter([sev_lvl](log::attribute_value_set const &s) {
		return s["Severity"].extract<trivial::severity_level>() >= sev_lvl;
	});
}

void install_gsl_error_handler() {
	gsl_set_error_handler_off();
}

void log_startup_information(int argc, char **argv)
{
	LOG(info) << "shark is starting in " << gethostname();
	LOG(info) << "shark version: " << SHARK_VERSION;
	LOG(info) << "shark git version: " << git_sha1();
	LOG(info) << "shark has local changes: " << std::boolalpha << git_has_local_changes() << std::noboolalpha;
	LOG(info) << "shark was built on " << __DATE__ << " " __TIME__;
	LOG(info) << "shark running at: " << boost::filesystem::current_path().string();

	std::ostringstream os;
	std::copy(argv, argv + argc, std::ostream_iterator<char *>(os, " "));
	LOG(info) << "shark started with command line: " << os.str();

	os = std::ostringstream();
	report_size_of_things(os);
	LOG(info) << os.str();
}


boost::program_options::variables_map parse_cmdline(int argc, char **argv) {

	using std::string;
	using std::vector;
	namespace po = boost::program_options;

	po::options_description visible_opts("SHArk options");
	visible_opts.add_options()
		("help,h",      "Show this help message")
		("version,V",   "Show version and exit")
		("verbose,v",   po::value<int>()->default_value(3), "Verbosity level. Higher is more verbose")
#ifdef SHARK_OPENMP
		("threads,t",   po::value<unsigned int>()->default_value(1), "OpenMP threads, defaults to 1. 0 means use OpenMP default number of threads")
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
	boost::program_options::variables_map vm;
	po::command_line_parser parser(argc, argv);
	parser.options(all_opts).positional(pdesc);
	po::store(parser.run(), vm);
	po::notify(vm);

	if (vm.count("help") != 0) {
		vm.clear();
		show_help(argv[0], visible_opts, std::cout);
	}
	else if (vm.count("version") != 0) {
		vm.clear();
		std::cout << "SHArk version " << SHARK_VERSION << std::endl;
		std::cout << "SHArk git revision " << git_sha1() << std::endl;
		std::cout << "SHArk has local changes: " << std::boolalpha << git_has_local_changes() << std::noboolalpha << std::endl;
		std::cout << "OpenMP support: ";
#ifdef SHARK_OPENMP
		std::cout << "Yes" << std::endl;
#else
		std::cout << "No" << std::endl;
#endif
		report_size_of_things(std::cout);
		std::cout << '\n';
	}
	else if (vm.count("config-file") == 0 ) {
		throw boost::program_options::error("At least one <config-file> option must be given. Use -h to see the help");
	}

	return vm;
}

Options read_options(const boost::program_options::variables_map &vm, unsigned int &threads) {

#ifdef SHARK_OPENMP
	threads = vm["threads"].as<unsigned int>();
	if (threads == 0) {
		threads = static_cast<unsigned int>(omp_get_max_threads());
	}
#else
	threads = 1;
#endif // SHARK_OPENMP
	LOG(info) << "shark using " << threads << " thread(s)";

	// Read the configuration file, and override options with any given
	// on the command-line
	Options options;
	for (auto &config_file: vm["config-file"].as<std::vector<std::string>>()) {
		options.add_file(config_file);
	}
	for(auto &opt_spec: vm["options"].as<std::vector<std::string>>()) {
		options.add(opt_spec);
	}

	return options;
}

int main(int argc, char **argv) {

	try {
		boost::program_options::variables_map vm = parse_cmdline(argc, argv);
		if (vm.empty()) {
			return 0;
		}

		setup_logging(vm);
		log_startup_information(argc, argv);
		install_gsl_error_handler();

		Timer timer;
		unsigned int threads;
		auto options = read_options(vm, threads);
		SharkRunner(options, threads).run();
		LOG(info) << "Successfully finished in " << timer;
		LOG(info) << "Maximum memory usage: " << memory_amount(peak_rss());

		return 0;
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

} // namespace shark

int main(int argc, char **argv) {
	return shark::main(argc, argv);
}