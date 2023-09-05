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
 * Definition of classes that output galaxies into disk
 */

#ifndef SHARK_GALAXY_WRITER_H_
#define SHARK_GALAXY_WRITER_H_

#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "agn_feedback.h"
#include "components.h"
#include "cosmology.h"
#include "dark_matter_halos.h"
#include "execution.h"
#include "hdf5/io/writer.h"
#include "simulation.h"
#include "star_formation.h"

namespace shark {

class GalaxyWriter {

public:

	GalaxyWriter(ExecutionParameters exec_params,
			CosmologicalParameters cosmo_params,
			CosmologyPtr cosmology,
			DarkMatterHalosPtr darkmatterhalo,
			SimulationParameters sim_params,
			AGNFeedbackParameters agn_params);
	virtual ~GalaxyWriter() = default;

	virtual void write(int snapshot, const std::vector<HaloPtr> &halos, TotalBaryon &AllBaryons, const molgas_per_galaxy &molgas_per_gal) = 0;

	void track_total_baryons(int snapshot, const std::vector<HaloPtr> &halos);

protected:

	ExecutionParameters exec_params;
	CosmologicalParameters cosmo_params;
	CosmologyPtr cosmology;
	DarkMatterHalosPtr darkmatterhalo;
	SimulationParameters sim_params;
	AGNFeedbackParameters agn_params;

	std::string get_output_directory(int snapshot);
};

class HDF5GalaxyWriter : public GalaxyWriter {

public:
	using GalaxyWriter::GalaxyWriter;
	void write(int snapshot, const std::vector<HaloPtr> &halos, TotalBaryon &AllBaryons, const molgas_per_galaxy &molgas_per_gal) override;

private:
	void write_header (hdf5::Writer &file, int snapshot);
	void write_galaxies (hdf5::Writer &file, int snapshot, const std::vector<HaloPtr> &halos, const molgas_per_galaxy &molgas_per_gal);
	void write_global_properties (hdf5::Writer &file, int snapshot, TotalBaryon &AllBaryons);
	void write_sf_histories (int snapshot, const std::vector<HaloPtr> &halos);
	void write_bh_histories (int snapshot, const std::vector<HaloPtr> &halos);
};

class ASCIIGalaxyWriter : public GalaxyWriter {

public:
	using GalaxyWriter::GalaxyWriter;
	void write(int snapshot, const std::vector<HaloPtr> &halos, TotalBaryon &AllBaryons, const molgas_per_galaxy &molgas_per_gal) override;

private:
	void write_galaxy(const Galaxy &galaxy, const SubhaloPtr &subhalo, int snapshot, std::ofstream &f, const molgas_per_galaxy &molgas_per_gal);

};

using GalaxyWriterPtr = std::unique_ptr<GalaxyWriter>;

template <typename ...Ts>
GalaxyWriterPtr make_galaxy_writer(const ExecutionParameters &exec_params, Ts&&...ts)
{
	if (exec_params.output_format == Options::HDF5) {
		return GalaxyWriterPtr(new HDF5GalaxyWriter(exec_params, std::forward<Ts>(ts)...));
	}
	else if (exec_params.output_format == Options::ASCII) {
		return GalaxyWriterPtr(new ASCIIGalaxyWriter(exec_params, std::forward<Ts>(ts)...));
	}

	std::ostringstream os;
	os << "Output format " << exec_params.output_format << " not currently supported";
	throw invalid_argument(os.str());
}

} // namespace shark

#endif /* SHARK_GALAXY_WRITER_H_ */
