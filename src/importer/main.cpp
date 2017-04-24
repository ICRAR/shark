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

#include <chrono>
#include <iostream>
#include <vector>
#include <memory>

#include "importer/options.h"
#include "importer/descendants.h"
#include "importer/velociraptor.h"

using namespace std;
using namespace shark;
using namespace shark::importer;

int main(int argc, char **argv)
{
	using chrono::steady_clock;
	using importer::Options;

	if ( argc < 2 ) {
		cerr << "Usage: " << argv[0] << " <options-file>" << endl;
		return 1;
	}

	Options opts(argv[1]);

	unique_ptr<Reader> reader;
	shared_ptr<DescendantReader> descendants_reader;

	if ( opts.descendants_format == Options::HDF5 ) {
		descendants_reader = make_shared<HDF5DescendantReader>(opts.descendants_file);
	}
	else if ( opts.descendants_format == Options::ASCII ) {
		descendants_reader = make_shared<AsciiDescendantReader>(opts.descendants_file);
	}

	if ( opts.tree_format == Options::TREES_VELOCIRAPTOR ) {
		reader.reset(new VELOCIraptorReader(descendants_reader, opts.tree_dir));
	}
	else if ( opts.tree_format == Options::TREES_NIFTY ) {
		// something else
	}
	else {
		// something else
	}

	for(int snapshot=opts.last_snapshot; snapshot >= opts.first_snapshot; snapshot--) {
		steady_clock::time_point t0 = steady_clock::now();
		auto subhalos = reader->read_subhalos(snapshot);
		long duration = chrono::duration_cast<chrono::milliseconds>(steady_clock::now() - t0).count();
		cout << "Snapshot " << snapshot << " read and processed in " << duration << " [ms]" << endl;
	}

}
