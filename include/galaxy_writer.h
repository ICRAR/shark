//
// Definition of classes that output galaxies into disk
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

#ifndef SHARK_GALAXY_WRITER_H_
#define SHARK_GALAXY_WRITER_H_

#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "components.h"
#include "cosmology.h"
#include "execution.h"
#include "hdf5/writer.h"
#include "simulation.h"
#include "star_formation.h"

namespace shark {

class GalaxyWriter {

public:

	GalaxyWriter(ExecutionParameters exec_params, CosmologicalParameters cosmo_params,  std::shared_ptr<Cosmology> cosmology, SimulationParameters sim_params, StarFormation starformation);
	virtual ~GalaxyWriter() {};

	virtual void write(int snapshot, const std::vector<HaloPtr> &halos, TotalBaryon &AllBaryons) = 0;

	void track_total_baryons(int snapshot, const std::vector<HaloPtr> &halos);

protected:

	ExecutionParameters exec_params;
	CosmologicalParameters cosmo_params;
	std::shared_ptr<Cosmology> cosmology;
	SimulationParameters sim_params;
	StarFormation starformation;

	std::string get_output_directory(int snapshot);
};

class HDF5GalaxyWriter : public GalaxyWriter {

public:
	using GalaxyWriter::GalaxyWriter;
	void write(int snapshot, const std::vector<HaloPtr> &halos, TotalBaryon &AllBaryons) override;
	void write_header (hdf5::Writer file, int snapshot);
	void write_galaxies (hdf5::Writer file, int snapshot, const std::vector<HaloPtr> &halos);
	void write_global_properties (hdf5::Writer file, int snapshot, TotalBaryon &AllBaryons);
	void write_histories (int snapshot, const std::vector<HaloPtr> &halos);
};

class ASCIIGalaxyWriter : public GalaxyWriter {

public:
	using GalaxyWriter::GalaxyWriter;
	void write(int snapshot, const std::vector<HaloPtr> &halos, TotalBaryon &AllBaryons) override;

private:
	void write_galaxy(const GalaxyPtr &galaxy, const SubhaloPtr &subhalo, int snapshot, std::ofstream &f);

};

}


#endif /* SHARK_GALAXY_WRITER_H_ */
