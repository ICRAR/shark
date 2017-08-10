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

#include "components.h"
#include "dark_matter_halos.h"
#include "options.h"

namespace shark {

class GalaxyMergerParameters {

	public:
		GalaxyMergerParameters(const Options &options);

		/**
		 * Merger parameters:
		 * - major_merger_ratio: threshold M2/M1 to consider major mergers. In this case we convert disks to spheroids.
		 * - minor_merger_burst_ratio: threshold M2/M1 for triggering bursts in minor mergers.
		 * - gas_fraction_minor_merger: suppress minor-merger bursts if fgas<gas_fraction_minor_merger in primary disk.
		 * - merger_random_seed: merger random seed to draw orbital parameters from Benson+05.
		 */

		float major_merger_ratio;
		float minor_merger_burst_ratio;
		float gas_fraction_minor_merger;
		int merger_random_seed;
		std::vector<double> jiang08;

};



class GalaxyMergers{

public:
	GalaxyMergers(GalaxyMergerParameters parameters, std::shared_ptr<DarkMatterHalos> darkmatterhalo);

	void orbital_parameters(double &vr, double &vt, double f);

	double mass_ratio_function(double mp, double ms);

	double merging_timescale_mass(double mp, double ms);

	double merging_timescale_orbital(double vr, double vt, double f, double c);

	double merging_timescale(std::shared_ptr<Subhalo> &primary, std::shared_ptr<Subhalo> &secondary);

	void merging_subhalos(std::shared_ptr<Halo> &halo);

	void merging_galaxies(std::shared_ptr<Halo> &halo);

private:
	GalaxyMergerParameters parameters;
	std::shared_ptr<DarkMatterHalos> darkmatterhalo;

};

} // namespace shark

#endif /* INCLUDE_GALAXY_MERGERS_H_ */
