/*
 * disk_instability.h
 *
 *  Created on: 29Sep.,2017
 *      Author: clagos
 */

#ifndef INCLUDE_DISK_INSTABILITY_H_
#define INCLUDE_DISK_INSTABILITY_H_

#include "agn_feedback.h"
#include "components.h"
#include "galaxy_mergers.h"
#include "physical_model.h"
#include "simulation.h"

namespace shark {

class DiskInstabilityParameters {

public:
	DiskInstabilityParameters(const Options &options);

	float stable = 0;
	float fint = 0;
};

class DiskInstability{

public:
	DiskInstability (DiskInstabilityParameters parameters,
			GalaxyMergerParameters merger_params,
			SimulationParameters simparams,
			const DarkMatterHalosPtr &darkmatterhalo,
			std::shared_ptr<BasicPhysicalModel> physicalmodel,
			const AGNFeedbackPtr &agnfeedback);

	double bulge_size(GalaxyPtr &galaxy);

	double toomre_parameter(GalaxyPtr &galaxy);

	void evaluate_disk_instability (HaloPtr &halo, int snapshot, double delta_t);

	void create_starburst(SubhaloPtr &subhalo, GalaxyPtr &galaxy, double z, double delta_t);

	void transfer_history_disk_to_bulge(GalaxyPtr &central, int snapshot);

	void effective_angular_momentum(GalaxyPtr &galaxy);

private:
	DiskInstabilityParameters parameters;
	GalaxyMergerParameters merger_params;
	SimulationParameters simparams;
	DarkMatterHalosPtr darkmatterhalo;
	std::shared_ptr<BasicPhysicalModel> physicalmodel;
	AGNFeedbackPtr agnfeedback;

};

} // end namespace


#endif /* INCLUDE_DISK_INSTABILITY_H_ */
