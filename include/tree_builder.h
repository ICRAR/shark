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
 * Merger tree builder classes
 */

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
	TreeBuilder(ExecutionParameters exec_params, unsigned int threads);
	virtual ~TreeBuilder();
	std::vector<MergerTreePtr> build_trees(std::vector<HaloPtr> &halos,
			SimulationParameters sim_params,
			GasCoolingParameters gas_cooling_params,
			DarkMatterHaloParameters dark_matter_params,
			const DarkMatterHalosPtr &darkmatterhalos,
			const CosmologyPtr &cosmology,
			TotalBaryon &AllBaryons);

protected:

	ExecutionParameters &get_exec_params();

	virtual void loop_through_halos(std::vector<HaloPtr> &halos) = 0;

	void link(const SubhaloPtr &parent_shalo, const SubhaloPtr &desc_subhalo,
	          const HaloPtr &parent_halo, const HaloPtr &desc_halo);

private:
	void ensure_trees_are_self_contained(const std::vector<MergerTreePtr> &trees) const;
	void ensure_halo_mass_growth(const std::vector<MergerTreePtr> &trees, SimulationParameters &sim_params);
	void spin_interpolated_halos(const std::vector<MergerTreePtr> &trees, SimulationParameters &sim_params);
	void define_central_subhalos(const std::vector<MergerTreePtr> &trees, SimulationParameters &sim_params, DarkMatterHaloParameters &dark_matter_params);
	SubhaloPtr define_central_subhalo(HaloPtr &halo, SubhaloPtr &subhalo);
	void define_accretion_rate_from_dm(const std::vector<MergerTreePtr> &trees, SimulationParameters &sim_params, GasCoolingParameters &gas_cooling_params, Cosmology &cosmology, TotalBaryon &AllBaryons);
	void remove_satellite(HaloPtr &halo, SubhaloPtr &subhalo);
 	void define_ages_halos(const std::vector<MergerTreePtr> &trees, SimulationParameters &sim_params, const DarkMatterHalosPtr &darkmatterhalos);
	void ignore_late_massive_halos(std::vector<MergerTreePtr> &trees,  SimulationParameters sim_params, ExecutionParameters exec_params);

private:
	ExecutionParameters exec_params;
	unsigned int threads = 1;
};


class HaloBasedTreeBuilder : public TreeBuilder {

public:
	HaloBasedTreeBuilder(ExecutionParameters exec_params, unsigned int threads);

protected:
	void loop_through_halos(std::vector<HaloPtr> &halos) override;
	SubhaloPtr find_descendant_subhalo(const HaloPtr &halo, const SubhaloPtr &subhalo, const HaloPtr &descendant_halo);

};

}  // namespace shark

#endif // SHARK_TREE_BUILDER_H_
