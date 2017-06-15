/*
 * agn_feedback.h
 *
 *  Created on: 14Jun.,2017
 *      Author: clagos
 */

#ifndef INCLUDE_AGN_FEEDBACK_H_
#define INCLUDE_AGN_FEEDBACK_H_


#include <memory>
#include <string>
#include <vector>

#include "cosmology.h"
#include "options.h"
#include "components.h"

namespace shark {

class AGNFeedbackParameters : public Options {

public:
	AGNFeedbackParameters(const std::string &filename);

	double mseed;
	double alpha_cool;
	double epsilon_smbh;
	double accretion_eff_cooling;
	double accretion_eff_bursts;
	double mhalo_seed;
};


class AGNFeedback {

public:
	AGNFeedback(AGNFeedbackParameters parameters, std::shared_ptr<Cosmology> cosmology);

	/**
	 * All input quantities should be in comoving units.
	 */

	void plant_seed_smbh(std::shared_ptr<Subhalo> &subhalo);
	double eddington_luminosity(double mbh);

private:
	AGNFeedbackParameters parameters;
	std::shared_ptr<Cosmology> cosmology;
};

}


#endif /* INCLUDE_AGN_FEEDBACK_H_ */
