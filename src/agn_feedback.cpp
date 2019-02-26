//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2017
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

#include <cmath>
#include <memory>

#include "agn_feedback.h"
#include "numerical_constants.h"
#include "utils.h"

namespace shark {

AGNFeedbackParameters::AGNFeedbackParameters(const Options &options)
{
	options.load("agn_feedback.mseed",mseed);
	options.load("agn_feedback.mhalo_seed",mhalo_seed);

	options.load("agn_feedback.model", model);

	//relevant for Bower06 model
	options.load("agn_feedback.alpha_cool",alpha_cool);
	options.load("agn_feedback.f_edd",f_edd);
	options.load("agn_feedback.accretion_eff_cooling",accretion_eff_cooling);

	//control accretion rate onto BHs during starbursts
	options.load("agn_feedback.f_smbh", f_smbh);
	options.load("agn_feedback.v_smbh", v_smbh);
	options.load("agn_feedback.tau_fold", tau_fold);

	// relevant for Croton16 model.
	options.load("agn_feedback.kappa_agn", kappa_agn);
	options.load("agn_feedback.accretion_eff_cooling", nu_smbh);

	// control QSO feedback.
	options.load("agn_feedback.qso_feedback", qso_feedback, false);
	options.load("agn_feedback.kappa_qso", kappa_qso);
	options.load("agn_feedback.epsilon_qso", epsilon_qso);

}

template <>
AGNFeedbackParameters::AGNFeedbackModel
Options::get<AGNFeedbackParameters::AGNFeedbackModel>(const std::string &name, const std::string &value) const {
	auto lvalue = lower(value);
	if (lvalue == "bower06") {
		return AGNFeedbackParameters::BOWER06;
	}
	else if (lvalue == "croton16") {
		return AGNFeedbackParameters::CROTON16;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value << ". Supported values are Bower06 and Croton16";
	throw invalid_option(os.str());
}

AGNFeedback::AGNFeedback(const AGNFeedbackParameters &parameters, const CosmologyPtr &cosmology, const RecyclingParameters &recycle_params) :
	parameters(parameters),
	cosmology(cosmology),
	recycle_params(recycle_params)
{
	// no-op
}

void AGNFeedback::plant_seed_smbh(Halo &halo){

	if (halo.Mvir > parameters.mhalo_seed) {
		auto central = halo.central_subhalo->central_galaxy();
		if (central && central->smbh.mass == 0) {
			central->smbh.mass = parameters.mseed;
			central->smbh.mass_metals = 0;
		}
	}

}

double AGNFeedback::eddington_luminosity(double mbh){

	// Numerical constants appearing in the expression for the Eddington luminosity (4 pi c G M_sun M_H/sigma_T)

	if(mbh >0){
		//Eddington luminosity in units of $10^{40} erg/s$
		return constants::Eddngtn_Lmnsty_Scale_Factor * mbh / constants::ERG2J;
	}
	else {
		return 0;
	}
}

double AGNFeedback::accretion_rate_hothalo_smbh(double Lcool, double mbh) {

	/**
	 * Function calculates the accretion rate onto the central black hole based on a cooling luminosity.
	 * Inputs:
	 * Lcool: cooling luminosity in units of 10^40 erg/s.
	 * mbh: mass of central supermassive black hole.
	 */

	using namespace constants;

	if (Lcool > 0 && Lcool < MAXLUM) {
		double macc = 0;
		if (parameters.model == AGNFeedbackParameters::BOWER06) {
			macc = Lcool * 1e40 / std::pow(c_light_cm,2.0) / parameters.accretion_eff_cooling;
		}
		else if (parameters.model == AGNFeedbackParameters::CROTON16) {
			macc = parameters.kappa_agn * 0.9375 * PI * G_cgs * M_Atomic_g * mu_Primordial * Lcool * 1e40 * (mbh * MSOLAR_g);
		}
		return macc * MACCRETION_cgs_simu; //accretion rate in units of Msun/Gyr.
	}
	else{
		return 0;
	}

}

double AGNFeedback::agn_bolometric_luminosity(double macc) {

	//return bolometric luminosity in units of 10^40 erg/s.
	using namespace constants;

	double Lbol = parameters.nu_smbh * (macc/MACCRETION_cgs_simu) * std::pow(c_light_cm,2.0) / 1e40;

	return Lbol;
}

double AGNFeedback::accretion_rate_hothalo_smbh_limit(double mheatrate, double vvir){

	using namespace constants;

	double Lbol = mheatrate * (0.5*std::pow(vvir*KM2CM,2.0)) / MACCRETION_cgs_simu / 1e40;

	double macc = Lbol / (parameters.nu_smbh) / std::pow(c_light_cm,2.0) * 1e40 * MACCRETION_cgs_simu;

	return macc;

}

double AGNFeedback::smbh_growth_starburst(double mgas, double vvir){

	double m = 0;

	if(mgas > 0){
		m =  parameters.f_smbh * mgas / (1 + std::pow(parameters.v_smbh/vvir, 2.0));
	}

	return m;
}

double AGNFeedback::smbh_accretion_timescale(Galaxy &galaxy, double z){

	double vbulge = std::sqrt(constants::G * galaxy.bulge_mass() / galaxy.bulge_gas.rscale);

	double tdyn = constants::MPCKM2GYR * cosmology->comoving_to_physical_size(galaxy.bulge_gas.rscale, z) / vbulge;

	return tdyn * parameters.tau_fold;

}

double AGNFeedback::qso_critical_luminosity(double mgas, double m, double r){

	double fgas =mgas/m;

	double sigma_bulge = std::sqrt(constants::G * m / r);

	//expression gives luminosity in 10^40 ergs/s.

	double Lm = 3e6 * (fgas/0.1) * std::pow((sigma_bulge/200.0), 4.0);

	return Lm;

}

double AGNFeedback::salpeter_timescale(double Lbol, double mbh){

	double ledd = eddington_luminosity(mbh);

	double edd_ratio = Lbol/ledd;

	//returns Salpeter timescale in Gyr.

	return 43.0 / edd_ratio / constants::KILO;
}

double AGNFeedback::qso_outflow_velocity(double Lbol, double zgas, double mgas){

	double vout  = 320.0 * std::pow(Lbol/(1e7 * constants::LSOLAR), 0.5) * std::pow(zgas/recycle_params.zsun, 0.25) * std::pow(mgas, -0.25);

	return vout;

}

void AGNFeedback::qso_outflow_rate(double mgas, double macc, double mBH, double zgas, double vcirc,
		double sfr, double mbulge, double rbulge, double &beta_halo, double &beta_ejec){

	// QSO feedback only acts if the accretion rate is >0, BH mass is > 0 and QSO feedback is activated by the user.
	if(macc > 0 and mBH > 0 and sfr > 0 and mgas > 0 and parameters.qso_feedback){
		double Lbol = agn_bolometric_luminosity(macc);
		double Lcrit = qso_critical_luminosity(mgas, mbulge, rbulge);

		// check if bolometric luminosity is larger than the critical luminosity and the gas mass in the bulge is positive. The latter is not always the case becaus equations are solved
		// numerically and hence negative solutions are in principle possible.
		if(Lbol > parameters.kappa_qso * Lcrit){

			double tsalp = salpeter_timescale(Lbol, mBH);
			double vout = qso_outflow_velocity(Lbol, zgas, mgas);

			double mout_rate= mgas/tsalp;

			double mejec_rate = (parameters.epsilon_qso * std::pow(vout/vcirc, 2.0) - 1) * mout_rate;

			// Apply boundary conditions to outflow and ejection rates
			if(mout_rate <  0 or std::isnan(mout_rate)){
				mout_rate = 0;
			}
			if(mejec_rate <  0 or std::isnan(mejec_rate)){
				mejec_rate = 0;
			}
			/*if(mejec_rate > mout_rate){
				mejec_rate = mout_rate;
			}*/

			beta_halo = mout_rate/sfr;
			beta_ejec = mejec_rate/sfr;
		}
		else{
			beta_halo = 0;
			beta_ejec = 0;
		}
	}
	else{
		beta_halo = 0;
		beta_ejec = 0;
	}

}

} // namespace shark
