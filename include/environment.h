//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2018
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
 */

#ifndef INCLUDE_ENVIRONMENT_H_
#define INCLUDE_ENVIRONMENT_H_



#include <memory>
#include <utility>

#include "baryon.h"
#include "components.h"
#include "cosmology.h"
#include "dark_matter_halos.h"
#include "options.h"
#include "root_solver.h"
#include "simulation.h"

namespace shark {

class EnvironmentParameters{

public:
	explicit EnvironmentParameters(const Options &options);

	bool gradual_stripping_halo = false;
	bool gradual_stripping_ism = false;
	bool stripping = true;
	bool tidal_stripping = false;
	float minimum_halo_mass_fraction = 0.01;
	float alpha_rps_halo = 1;
	float Accuracy_RPS = 0.05;

};

class Environment{

public:
	explicit Environment(const EnvironmentParameters &parameters,
			DarkMatterHalosPtr darkmatterhalos,
			CosmologyPtr cosmology,
			CosmologicalParameters cosmo_params,
			SimulationParameters simparams);

	void process_satellite_subhalo_environment (Subhalo &satellite_subhalo, SubhaloPtr &central_subhalo, double z);
	BaryonBase remove_tidal_stripped_stars(SubhaloPtr &subhalo, Galaxy &galaxy, BaryonBase lost_stellar);

	double process_ram_pressure_stripping_gas(const SubhaloPtr &primary,
			Subhalo &secondary,
			double z,
			double ram_press,
			bool halo_strip,
			bool ism_strip);

	double ram_pressure_stripping_hot_gas(const SubhaloPtr &primary,
			const Subhalo &secondary,
			double r,
			double z,
			double ram_press);

	double ram_pressure_stripping_galaxy_gas(const GalaxyPtr &galaxy,
			double r,
			double z,
			double ram_press);

	double ram_pressure(const SubhaloPtr &primary,
			const Subhalo &secondary,
			double z);

	using func_x = double (*)(double x, void *);

private:

	EnvironmentParameters parameters;
	DarkMatterHalosPtr darkmatterhalos;
	CosmologyPtr cosmology;
	CosmologicalParameters cosmo_params;
	SimulationParameters simparams;
	Root_Solver root_solver;

	float remove_gas(BaryonBase &component, double m_removed, float f_gas);

};

using EnvironmentPtr = std::shared_ptr<Environment>;

template <typename ...Ts>
EnvironmentPtr make_environment(Ts&&...ts)
{
	return std::make_shared<Environment>(std::forward<Ts>(ts)...);
}


}//end namespace shark


#endif /* INCLUDE_ENVIRONMENT_H_ */
