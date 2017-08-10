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
#include "importer/options.h"
#include "importer/descendants.h"
#include "importer/velociraptor.h"

using namespace std;
using namespace shark;
using namespace shark::importer;

int main(int argc, char **argv)
{
	using importer::Options;

	if ( argc < 2 ) {
		cerr << "Usage: " << argv[0] << " <options-file>" << endl;
		return 1;
	}

	shark::Options options(argv[1]);
	Options importer_opts(options);

	//
	// The reader for the descendants file
	//
	shared_ptr<DescendantReader> descendants_reader;
	if ( importer_opts.descendants_format == shark::Options::HDF5 ) {
		descendants_reader = make_shared<HDF5DescendantReader>(importer_opts.descendants_file);
	}
	else if ( importer_opts.descendants_format == shark::Options::ASCII ) {
		descendants_reader = make_shared<AsciiDescendantReader>(importer_opts.descendants_file);
	}

	//
	// The tree reader
	//
	if ( importer_opts.tree_format != Options::TREES_VELOCIRAPTOR ) {
		throw invalid_option("Only tree format currently supported is VELOCIraptor");
	}
	unique_ptr<Reader> reader(new VELOCIraptorReader(descendants_reader, importer_opts.tree_dir));

	//
	// Go ahead and read all required snapshots
	//
	for(int snapshot=importer_opts.last_snapshot; snapshot >= importer_opts.first_snapshot; snapshot--) {
		Timer timer;
		auto subhalos = reader->read_subhalos(snapshot);
		cout << "Snapshot " << snapshot << " read and processed in " << timer.get() << " [ms]" << endl;
	}

}
