/*
 * agn_feedback.cpp
 *
 *  Created on: 15Jun.,2017
 *      Author: clagos
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
	options.load("agn_feedback.alpha_cool",alpha_cool);
	options.load("agn_feedback.f_edd",f_edd);

	options.load("agn_feedback.f_smbh", f_smbh);
	options.load("agn_feedback.v_smbh", v_smbh);
	options.load("agn_feedback.tau_fold", tau_fold);

	options.load("agn_feedback.accretion_eff_cooling",accretion_eff_cooling);

	options.load("agn_feedback.kappa_agn", kappa_agn);
	options.load("agn_feedback.accretion_eff_cooling", nu_smbh);

	options.load("agn_feedback.mass_thresh", mass_thresh);

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

AGNFeedback::AGNFeedback(const AGNFeedbackParameters &parameters, const CosmologyPtr &cosmology) :
	parameters(parameters),
	cosmology(cosmology)
{
	// no-op
}

void AGNFeedback::plant_seed_smbh(Halo &halo){

	if (halo.Mvir > parameters.mhalo_seed) {
		auto central = halo.central_subhalo->central_galaxy();
		if (central and central->smbh.mass == 0) {
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

} // namespace shark
