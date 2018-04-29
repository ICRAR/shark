/*
 * galaxy_mergers.h
 *
 *  Created on: 4Aug.,2017
 *      Author: clagos
 */

#ifndef INCLUDE_GALAXY_MERGERS_H_
#define INCLUDE_GALAXY_MERGERS_H_

#include <memory>
#include <vector>

#include "agn_feedback.h"
#include "components.h"
#include "dark_matter_halos.h"
#include "options.h"
#include "physical_model.h"
#include "simulation.h"

namespace shark {

class GalaxyMergerParameters {

	public:
		GalaxyMergerParameters(const Options &options);

		/**
		 * Merger parameters:
		 * - major_merger_ratio: threshold M2/M1 to consider major mergers. In this case we convert disks to spheroids.
		 * - minor_merger_burst_ratio: threshold M2/M1 for triggering bursts in minor mergers.
		 * - merger_random_seed: merger random seed to draw orbital parameters from Benson+05.
		 * - tau_delay: controls delays from the standard merging timescale for testing purposes.
		 * - jiang08: parameters of the best fit dynamical timescale of Jiang et al. (2008).
		 * - min_mass: minimum mass allowed in bulges. This is to avoid long tails in the star formation histories of bulges.
		 */

		float major_merger_ratio;
		float minor_merger_burst_ratio;
		float gas_fraction_burst_ratio;
		int merger_random_seed;
		float tau_delay;
		std::vector<double> jiang08;
		float mass_min;

		/**
		 * Sizes parameters:
		 * - f_orbit: orbital factor defining orbital energy. It should be =1 for two point masses in a circular orbit separation rgal,1+rgal,2. Lacey et al. (2016) Eq. 18.
		 * - cgal: parameter that defines internal energy of galaxy. It depends weekly on density profile. =0.49 for a pure exponential disk; =0.45 for a De Vacouleurs profile.
		 * - fgas_dissipation: parameter that defines how much dissipation there is when calculating the galaxy sizes in mergers. A value of 0 is adopted if no dissipation takes place.
		 * - merger_ratio_dissipation: parameter that defines the merger mass ratio above which dissipation is triggered.
		 */

		float f_orbit;
		float cgal;
		float merger_ratio_dissipation;
		double fgas_dissipation;

};



class GalaxyMergers{

public:
	GalaxyMergers(GalaxyMergerParameters parameters,
			SimulationParameters simparams,
			std::shared_ptr<DarkMatterHalos> darkmatterhalo,
			std::shared_ptr<BasicPhysicalModel> physicalmodel,
			std::shared_ptr<AGNFeedback> agnfeedback);

	void orbital_parameters(double &vr, double &vt, double f);

	double mass_ratio_function(double mp, double ms);

	double merging_timescale_mass(double mp, double ms);

	double merging_timescale_orbital(double vr, double vt, double f, double c);

	void merging_timescale(SubhaloPtr &primary, SubhaloPtr &secondary, double z);

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
	SimulationParameters simparams;
	std::shared_ptr<DarkMatterHalos> darkmatterhalo;
	std::shared_ptr<BasicPhysicalModel> physicalmodel;
	std::shared_ptr<AGNFeedback> agnfeedback;

};

} // namespace shark

#endif /* INCLUDE_GALAXY_MERGERS_H_ */
