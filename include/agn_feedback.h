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

class AGNFeedbackParameters {

public:
	AGNFeedbackParameters(const Options &options);

	double mseed;
	double mhalo_seed;
	double alpha_cool;
	double f_edd;
	double accretion_eff_cooling;
	double accretion_eff_bursts;
	double kappa_agn;
	double nu_smbh;

	enum AGNFeedbackModel {
		LGALAXIES = 0,
		GALFORM
	};

	AGNFeedbackModel model;
};


class AGNFeedback {

public:
	AGNFeedback(AGNFeedbackParameters parameters, std::shared_ptr<Cosmology> cosmology);

	/**
	 * All input quantities should be in comoving units.
	 */

	void plant_seed_smbh(Subhalo &subhalo);
	double eddington_luminosity(double mbh);
	double accretion_rate_hothalo_smbh(double Lcool, double mbh);
	double agn_bolometric_luminosity(double macc);

	// TODO: move this to private when possible
	AGNFeedbackParameters parameters;

private:
	std::shared_ptr<Cosmology> cosmology;
};

}


#endif /* INCLUDE_AGN_FEEDBACK_H_ */
