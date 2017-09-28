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

namespace shark {

class GalaxyMergerParameters {

	public:
		GalaxyMergerParameters(const Options &options);

		/**
		 * Merger parameters:
		 * - major_merger_ratio: threshold M2/M1 to consider major mergers. In this case we convert disks to spheroids.
		 * - minor_merger_burst_ratio: threshold M2/M1 for triggering bursts in minor mergers.
		 * - merger_random_seed: merger random seed to draw orbital parameters from Benson+05.
		 */

		float major_merger_ratio;
		float minor_merger_burst_ratio;
		int merger_random_seed;
		std::vector<double> jiang08;

		/**
		 * Sizes parameters:
		 * - f_orbit: orbital factor defining orbital energy. It should be =1 for two point masses in a circular orbit separation rgal,1+rgal,2. Lacey et al. (2016) Eq. 18.
		 * - cgal: parameter that defines internal energy of galaxy. It depends weekly on density profile. =0.49 for a pure exponential disk; =0.45 for a De Vacouleurs profile.
		 */

		float f_orbit;
		float cgal;

};



class GalaxyMergers{

public:
	GalaxyMergers(GalaxyMergerParameters parameters, std::shared_ptr<DarkMatterHalos> darkmatterhalo, std::shared_ptr<BasicPhysicalModel> physicalmodel, std::shared_ptr<AGNFeedback> agnfeedback);

	void orbital_parameters(double &vr, double &vt, double f);

	double mass_ratio_function(double mp, double ms);

	double merging_timescale_mass(double mp, double ms);

	double merging_timescale_orbital(double vr, double vt, double f, double c);

	void merging_timescale(SubhaloPtr &primary, SubhaloPtr &secondary);

	void merging_subhalos(HaloPtr &halo);

	void merging_galaxies(HaloPtr &halo, double z, double delta_t);

	void create_merger(GalaxyPtr &central, GalaxyPtr &satellite, HaloPtr &halo);

	void create_starbursts(HaloPtr &halo, double z, double delta_t);

	double bulge_size_merger(double mass_ratio, GalaxyPtr &central, GalaxyPtr &satellite, HaloPtr &halo);

	double r_remnant(double mc, double ms, double rc, double rs);

	void transfer_baryon_mass(SubhaloPtr satellite, SubhaloPtr central);


private:
	GalaxyMergerParameters parameters;
	std::shared_ptr<DarkMatterHalos> darkmatterhalo;
	std::shared_ptr<BasicPhysicalModel> physicalmodel;
	std::shared_ptr<AGNFeedback> agnfeedback;

};

} // namespace shark

#endif /* INCLUDE_GALAXY_MERGERS_H_ */
