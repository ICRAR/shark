/*
 * agn_feedback.h
 *
 *  Created on: 14Jun.,2017
 *      Author: clagos
 */

#ifndef INCLUDE_AGN_FEEDBACK_H_
#define INCLUDE_AGN_FEEDBACK_H_

#include <memory>
#include <utility>

#include "cosmology.h"
#include "options.h"
#include "components.h"

namespace shark {

class AGNFeedbackParameters {

public:
	AGNFeedbackParameters(const Options &options);

	double mseed = 0;
	double mhalo_seed = 0;
	double alpha_cool = 0;
	double f_edd = 0;
	double f_smbh = 0;
	double v_smbh = 0;
	double tau_fold = 0;
	double accretion_eff_cooling = 0;
	double kappa_agn = 0;
	double nu_smbh = 0;
	double mass_thresh = 0;

	enum AGNFeedbackModel {
		CROTON16 = 0,
		BOWER06
	};

	AGNFeedbackModel model = BOWER06;
};


class AGNFeedback {

public:
	AGNFeedback(const AGNFeedbackParameters &parameters, const CosmologyPtr &cosmology);

	/**
	 * All input quantities should be in comoving units.
	 */

	void plant_seed_smbh(Halo &halo);
	double eddington_luminosity(double mbh);
	double accretion_rate_hothalo_smbh(double Lcool, double mbh);
	double agn_bolometric_luminosity(double macc);
	double smbh_growth_starburst(double mgas, double vvir);
	double smbh_accretion_timescale(Galaxy &galaxy, double z);
	double accretion_rate_hothalo_smbh_limit(double mheatrate, double vvir);

	// TODO: move this to private when possible
	AGNFeedbackParameters parameters;

private:
	CosmologyPtr cosmology;
};

/// Type used by users to handle an instance of AGNFeedback
typedef std::shared_ptr<AGNFeedback> AGNFeedbackPtr;

template <typename ...Ts>
AGNFeedbackPtr make_agn_feedback(Ts&&...ts)
{
	return std::make_shared<AGNFeedback>(std::forward<Ts>(ts)...);
}

} //end namespace shark

#endif /* INCLUDE_AGN_FEEDBACK_H_ */