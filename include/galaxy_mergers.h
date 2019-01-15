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

#ifndef INCLUDE_GALAXY_MERGERS_H_
#define INCLUDE_GALAXY_MERGERS_H_

#include <memory>
#include <random>
#include <vector>

#include "agn_feedback.h"
#include "components.h"
#include "cosmology.h"
#include "dark_matter_halos.h"
#include "options.h"
#include "physical_model.h"
#include "simulation.h"
#include "execution.h"

namespace shark {

class GalaxyMergerParameters {

	public:
		GalaxyMergerParameters(const Options &options);

		/**
		 * Merger parameters:
		 * - major_merger_ratio: threshold M2/M1 to consider major mergers. In this case we convert disks to spheroids.
		 * - minor_merger_burst_ratio: threshold M2/M1 for triggering bursts in minor mergers.
		 * - tau_delay: controls delays from the standard merging timescale for testing purposes.
		 * - min_mass: minimum mass allowed in bulges. This is to avoid long tails in the star formation histories of bulges.
		 */
		float major_merger_ratio = 0;
		float minor_merger_burst_ratio = 0;
		float gas_fraction_burst_ratio = 0;
		float tau_delay = 0.05;
		float mass_min = 1e5;

		/**
		 * Sizes parameters:
		 * - f_orbit: orbital factor defining orbital energy. It should be =1 for two point masses in a circular orbit separation rgal,1+rgal,2. Lacey et al. (2016) Eq. 18.
		 * - cgal: parameter that defines internal energy of galaxy. It depends weekly on density profile. =0.49 for a pure exponential disk; =0.45 for a De Vacouleurs profile.
		 * - fgas_dissipation: parameter that defines how much dissipation there is when calculating the galaxy sizes in mergers. A value of 0 is adopted if no dissipation takes place.
		 * - merger_ratio_dissipation: parameter that defines the merger mass ratio above which dissipation is triggered.
		 */

		float f_orbit = 1;
		float cgal = 0.5;
		float merger_ratio_dissipation = 0;
		double fgas_dissipation = 0;

};



class GalaxyMergers{

public:
	GalaxyMergers(GalaxyMergerParameters parameters,
			const CosmologyPtr &cosmology,
			const ExecutionParameters &execparams,
			SimulationParameters simparams,
			const DarkMatterHalosPtr &darkmatterhalo,
			std::shared_ptr<BasicPhysicalModel> physicalmodel,
			const AGNFeedbackPtr &agnfeedback);

	void orbital_parameters(double &vr, double &vt, double f);

	double mass_ratio_function(double mp, double ms);

	double merging_timescale_mass(double mp, double ms);

	double merging_timescale_orbital();

	/**
	 * Calculates the dynamical friction timescale for the subhalo secondary to merge into the subhalo primary,
	 * or of the satellite galaxies type 2 of the satellite subhalo.
	 * This should be calculated only in the snapshot before the secondary mergers onto the primary (i.e. disappears from merger tree).
	 *
	 * @param primary the primary subhalo, which must contain the central galaxy
	 * @param secondary the secondary subhalo
	 * @param z redshift
	 * @param transfer_types2 whether we are merging a satellite subhalo or transfering type 2 galaxies.
	 */
	void merging_timescale(SubhaloPtr &primary, SubhaloPtr &secondary, double z, bool transfer_types2);

	/**
	 * Evaluates whether subhalos in each timestep are disappearing from the merger tree, and if they are
	 * it passes that satellite galaxy onto the central subhalo and calculates a dynamical friction timescale that the
	 * galaxy will save in case it's a satellite.
	 *
	 * In the case of satellite subhalos that will continue to exist, we transfer
	 * any existing type2 galaxies to the central subhalo, after recalculating
	 * their dynamical friction timescale with respect to the central galaxy
	 * of the central subhalo.
	 *
	 * @param halo the halo where subhalos are going to be possibly merged
	 * @param z the redshift
	 */
	void merging_subhalos(HaloPtr &halo, double z);

	void merging_galaxies(HaloPtr &halo, int snapshot, double delta_t);

	void create_merger(GalaxyPtr &central, GalaxyPtr &satellite, HaloPtr &halo, int snapshot);

	void create_starbursts(HaloPtr &halo, double z, double delta_t);

	double bulge_size_merger(double mass_ratio, double mgas_ratio, GalaxyPtr &central, GalaxyPtr &satellite, HaloPtr &halo);

	double r_remnant(double mc, double ms, double rc, double rs);

	void transfer_baryon_mass(SubhaloPtr central, SubhaloPtr satellite);

	void transfer_bulge_gas(SubhaloPtr &subhalo, GalaxyPtr &galaxy, double z);

	void transfer_history_satellite_to_bulge(GalaxyPtr &central, GalaxyPtr &satellite, int snapshot);

	void transfer_history_disk_to_bulge(GalaxyPtr &central, int snapshot);


private:
	GalaxyMergerParameters parameters;
	std::shared_ptr<Cosmology> cosmology;
	SimulationParameters simparams;
	DarkMatterHalosPtr darkmatterhalo;
	std::shared_ptr<BasicPhysicalModel> physicalmodel;
	AGNFeedbackPtr agnfeedback;

	std::default_random_engine generator;
	std::lognormal_distribution<double> distribution;

};

} // namespace shark

#endif /* INCLUDE_GALAXY_MERGERS_H_ */
