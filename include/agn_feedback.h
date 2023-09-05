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
#include <random>
#include <utility>

#include "components.h"
#include "cosmology.h"
#include "execution.h"
#include "galaxy.h"
#include "options.h"
#include "recycling.h"
#include "subhalo.h"

namespace shark {


class AGNFeedbackParameters {

public:
	explicit AGNFeedbackParameters(const Options &options);

	double mseed = 0;
	double mhalo_seed = 0;

	double alpha_cool = 0;
	double f_edd = 0;
	double f_smbh = 0;
	double v_smbh = 0;
	double tau_fold = 0;
	double accretion_eff_cooling = 0;
	double kappa_agn = 0;
	double nu_smbh = 0.1;

	bool qso_feedback = false;
	double epsilon_qso = 0;

	// Parameters relevant to Lagos22 model
	float kappa_radio = 1;
	float eta_superedd = 4;
	float hot_halo_threshold = 1;
	float alpha_adaf = 0.1;
	float alpha_td = 0.1;
	float delta_adaf = 0.2;
	float mdotcrit_adaf = 0.01;
	float low_accretion_adaf = 1e-4;
	float constant_lowlum_adaf = 0;
	float constant_highlum_adaf = 0;
	float nu2_nu1 = 0;
	int loop_limit_accretion = 10;

	enum AGNFeedbackModel {
		CROTON16 = 0,
		BOWER06,
		LAGOS23
	};

	enum SpinModel{
		VOLONTERI07 = 0,
		GRIFFIN19,
		CONSTANTSPIN
	};

	enum AccretionDiskModel{
		WARPEDDISK = 0,
		SELFGRAVITYDISK,
		PROLONGED
	};

	AGNFeedbackModel model = CROTON16;
	SpinModel spin_model = VOLONTERI07;
	AccretionDiskModel accretion_disk_model;
};


class AGNFeedback {

public:
	AGNFeedback(const AGNFeedbackParameters &parameters, CosmologyPtr cosmology, RecyclingParameters recycle_params, ExecutionParameters exec_params);

	/**
	 * All input quantities should be in comoving units.
	 */

	void plant_seed_smbh(Subhalo &subhalo);
	double eddington_luminosity(double mbh);
	double accretion_rate_hothalo_smbh(double Lcool, double tacc, double fhot, double vvir, Galaxy &galaxy);
	double accretion_rate_hothalo_smbh_limit(double mheatrate, double vvir, const BlackHole &smbh);
	double accretion_rate_ratio(double macc, double mBH);
	double agn_bolometric_luminosity(const BlackHole &smbh,  bool starburst);
	double agn_mechanical_luminosity(const BlackHole &smbh);
	double smbh_growth_starburst(double mgas, double vvir, double tacc, Galaxy &galaxy);
	double smbh_accretion_timescale(Galaxy &galaxy, double z);
	double qso_critical_luminosity(double mgas, double m, double r);
	double salpeter_timescale(double Lbol, double mbh);
	double qso_outflow_velocity(double Lbol, double mbh, double zgas, double mgas, double mbulge, double rbulge);
	void qso_outflow_rate(double mgas, const BlackHole &smbh, double zgas, double vcirc,
			double sfr, double mbulge, double rbulge, double &beta_halo, double &beta_ejec);
	void griffin19_spinup_accretion(double delta_mbh, double tau_acc, Galaxy &galaxy);
	void griffin19_spinup_mergers(BlackHole &smbh_primary, const BlackHole &smbh_secondary, const Galaxy &galaxy);
	void volonteri07_spin(BlackHole &smbh);
	std::vector<float> efficiency_luminosity_agn(float spin);

	// TODO: move this to private when possible
	AGNFeedbackParameters parameters;

private:
	CosmologyPtr cosmology;
	RecyclingParameters recycle_params;
	ExecutionParameters exec_params;

	std::uniform_real_distribution<double> distribution;

	double angle_acc_disk(const Galaxy &galaxy);
	double angle_acc_disk(std::default_random_engine &generator);
	double phi_acc_disk(std::default_random_engine &generator);
	double final_spin(const double mbh, const double mfin, const double r_lso);

};

/// Type used by users to handle an instance of AGNFeedback
using AGNFeedbackPtr = std::shared_ptr<AGNFeedback>;

template <typename ...Ts>
AGNFeedbackPtr make_agn_feedback(Ts&&...ts)
{
	return std::make_shared<AGNFeedback>(std::forward<Ts>(ts)...);
}

} //end namespace shark

#endif /* INCLUDE_AGN_FEEDBACK_H_ */
