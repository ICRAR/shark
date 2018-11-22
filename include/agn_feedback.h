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

#ifndef INCLUDE_AGN_FEEDBACK_H_
#define INCLUDE_AGN_FEEDBACK_H_

#include <memory>
#include <utility>

#include "cosmology.h"
#include "components.h"
#include "options.h"
#include "recycling.h"

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
	double kappa_qso = 0;
	double epsilon_qso = 0;

	enum AGNFeedbackModel {
		CROTON16 = 0,
		BOWER06
	};

	AGNFeedbackModel model = BOWER06;
};


class AGNFeedback {

public:
	AGNFeedback(const AGNFeedbackParameters &parameters, const CosmologyPtr &cosmology, const RecyclingParameters &recycle_params);

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
	double qso_critical_luminosity(Galaxy &galaxy);
	double salpeter_timescale(double Lbol, double mbh);
	double qso_outflow_velocity(double Lbol, double zgas, double mgas);
	void qso_outflow_rate(double mgas, double tsalp, double vout, double vcirc, double sfr, double beta_halo, double beta_ejec);


	// TODO: move this to private when possible
	AGNFeedbackParameters parameters;

private:
	CosmologyPtr cosmology;
	RecyclingParameters recycle_params;

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
