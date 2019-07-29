//
// Main routine for the shark-importer program
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

#include <iostream>
#include <vector>
#include <memory>

#include "timer.h"
#include "options.h"
#include "subhalo.h"
#include "importer/descendants.h"
#include "importer/velociraptor.h"

namespace shark {
namespace importer {

class ImporterParameters {

public:

	/**
	 * Ctor providing default values
	 */
	explicit ImporterParameters(const Options &options) :
		descendants_format(Options::ASCII),
		descendants_file("descendants.txt"),
		tree_format(TREES_VELOCIRAPTOR),
		tree_dir("."),
		first_snapshot(0),
		last_snapshot(0)
	{
		options.load("input.tree_dir", tree_dir);
		options.load("input.tree_format", tree_format);
		options.load("input.descendants_format", descendants_format);
		options.load("input.descendants", descendants_file);
		options.load("input.first_snapshot", first_snapshot);
		options.load("input.last_snapshot", last_snapshot);
	}

	enum tree_format_t {
		TREES_VELOCIRAPTOR,
		TREES_NIFTY
	};

	Options::file_format_t descendants_format;
	std::string descendants_file;
	tree_format_t tree_format;
	std::string tree_dir;
	int first_snapshot;
	int last_snapshot;
};

} // namespace importer

template <>
importer::ImporterParameters::tree_format_t
Options::get<importer::ImporterParameters::tree_format_t>(const std::string &name, const std::string &value) const {
	if ( value == "velociraptor" ) {
		return importer::ImporterParameters::TREES_VELOCIRAPTOR;
	}
	else if ( value == "nifty" ) {
		return importer::ImporterParameters::TREES_NIFTY;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value;
	throw invalid_option(os.str());
}

namespace importer {

int run(int argc, char **argv)
{
	if ( argc < 2 ) {
		std::cerr << "Usage: " << argv[0] << " <options-file>" << std::endl;
		return 1;
	}

	Options options(argv[1]);
	ImporterParameters importer_params(options);

	//
	// The reader for the descendants file
	//
	std::shared_ptr<DescendantReader> descendants_reader;
	if ( importer_params.descendants_format == shark::Options::HDF5 ) {
		descendants_reader = std::make_shared<HDF5DescendantReader>(importer_params.descendants_file);
	}
	else if ( importer_params.descendants_format == shark::Options::ASCII ) {
		descendants_reader = std::make_shared<AsciiDescendantReader>(importer_params.descendants_file);
	}

	//
	// The tree reader
	//
	if ( importer_params.tree_format != ImporterParameters::TREES_VELOCIRAPTOR ) {
		throw invalid_option("Only tree format currently supported is VELOCIraptor");
	}
	std::unique_ptr<Reader> reader(new VELOCIraptorReader(descendants_reader, importer_params.tree_dir));

	//
	// Go ahead and read all required snapshots
	//
	for(int snapshot=importer_params.last_snapshot; snapshot >= importer_params.first_snapshot; snapshot--) {
		Timer timer;
		auto subhalos = reader->read_subhalos(snapshot);
		std::cout << "Snapshot " << snapshot << " read and processed in " << timer << "\n";
	}

	return 0;
}

} // namespace importer
} // namespace shark

int main(int argc, char *argv[]) {
	try {
		shark::importer::run(argc, argv);
		return 0;
	} catch (shark::invalid_option &) {
		return 1;
	} catch (...) {
		return 2;
	}
}