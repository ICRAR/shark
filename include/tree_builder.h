//
// Merger tree builder classes
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

#ifndef SHARK_TREE_BUILDER_H_
#define SHARK_TREE_BUILDER_H_

#include <vector>

#include "components.h"
#include "cosmology.h"
#include "dark_matter_halos.h"
#include "execution.h"
#include "gas_cooling.h"
#include "simulation.h"

namespace shark {

class TreeBuilder {

public:
	TreeBuilder(ExecutionParameters exec_params);
	virtual ~TreeBuilder();
	std::vector<MergerTreePtr> build_trees(const std::vector<HaloPtr> &halos, SimulationParameters sim_params, std::shared_ptr<Cosmology> cosmology);

protected:

	ExecutionParameters &get_exec_params();

	virtual void loop_through_halos(const std::vector<HaloPtr> &halos) = 0;

	void link(const SubhaloPtr &subhalo, const SubhaloPtr &d_subhalo,
	          const HaloPtr &halo, const HaloPtr &d_halo);

	void ensure_halo_mass_growth(std::vector<MergerTreePtr> trees, SimulationParameters sim_params);

	void spin_interpolated_halos(std::vector<MergerTreePtr> trees, SimulationParameters sim_params);

	void define_central_subhalos(std::vector<MergerTreePtr> trees, SimulationParameters sim_params);

	SubhaloPtr define_central_subhalo(HaloPtr &halo, SubhaloPtr &subhalo);

	void define_accretion_rate_from_dm(std::vector<MergerTreePtr> trees, SimulationParameters sim_params, Cosmology &cosmology);

	void remove_satellite(HaloPtr halo, SubhaloPtr subhalo);

private:
	ExecutionParameters exec_params;

};


class HaloBasedTreeBuilder : public TreeBuilder {

public:
	HaloBasedTreeBuilder(ExecutionParameters exec_params);

	void create_galaxies(HaloPtr halo,
			Cosmology &cosmology,
			DarkMatterHalos &darkmatterhalos,
			GasCoolingParameters &cool_params,
			SimulationParameters sim_params);

protected:
	virtual void loop_through_halos(const std::vector<HaloPtr> &halos) override;

};

}  // namespace shark

#endif // SHARK_TREE_BUILDER_H_
