/*
 * agn_feedback.cpp
 *
 *  Created on: 15Jun.,2017
 *      Author: clagos
 */

#include <cmath>
#include <memory>

#include "numerical_constants.h"
#include "agn_feedback.h"

namespace shark {

AGNFeedbackParameters::AGNFeedbackParameters(const std::string &filename) :
	Options(filename),
	mseed(0),
	mhalo_seed(0),
	alpha_cool(0),
	f_edd(0),
	accretion_eff_cooling(0),
	accretion_eff_bursts(0),
	kappa_agn(0),
	nu_smbh(0),
	model(GALFORM)
{
	load("agn_feedback.mseed",mseed);
	load("agn_feedback.mhalo_seed",mhalo_seed);

	load("agn_feedback.model", model);
	load("agn_feedback.alpha_cool",alpha_cool);
	load("agn_feedback.epsilon_smbh",f_edd);
	load("agn_feedback.accretion_eff_cooling",accretion_eff_cooling);
	load("agn_feedback.accretion_eff_bursts",accretion_eff_bursts);

	load("agn_feedback.kappa_agn", kappa_agn);
	load("agn_feedback.nu_smbh", nu_smbh);

}
namespace detail {

template <>
AGNFeedbackParameters::AGNFeedbackModel Helper<AGNFeedbackParameters::AGNFeedbackModel>::get(const std::string &name, const std::string &value) {
	if ( value == "galform" ) {
		return AGNFeedbackParameters::GALFORM;
	}
	else if ( value == "l-galaxies" ) {
		return AGNFeedbackParameters::LGALAXIES;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value << ". Supported values are galform and l-galaxies";
	throw invalid_option(os.str());
}
}

AGNFeedback::AGNFeedback(AGNFeedbackParameters parameters, std::shared_ptr<Cosmology> cosmology) :
	parameters(parameters),
	cosmology(cosmology)
{
	// no-op
}

void AGNFeedback::plant_seed_smbh(Subhalo &subhalo){

	if(subhalo.Mvir > parameters.mhalo_seed){
		for(auto &galaxies: subhalo.galaxies) {
			if(galaxies->galaxy_type == Galaxy::CENTRAL && galaxies->smbh.mass == 0){
				galaxies->smbh.mass = 0;
				galaxies->smbh.mass_metals = 0;
			}
		}
	}
}

double AGNFeedback::eddington_luminosity(double mbh){

	// Numerical constants appearing in the expression for the Eddington luminosity (4 pi c G M_sun M_H/sigma_T)

	if(mbh >0){
		//Eddington luminosity in units of $10^{40} erg/s$
		return constants::Eddngtn_Lmnsty_Scale_Factor*mbh/constants::ERG2J;
	}

	return 0;
}

double AGNFeedback::accretion_rate_hothalo_smbh(double Lcool, double mbh) {

	/**
	 * Input Lcool: cooling luminosity in units of 10^40 erg/s.
	 */
	using namespace constants;

	if (Lcool > 0) {
		double macc;
		if (parameters.model == AGNFeedbackParameters::GALFORM) {
			macc = Lcool * std::pow(10,40) / std::pow(c_light_cm,2) / parameters.accretion_eff_cooling;
		}
		else if (parameters.model == AGNFeedbackParameters::LGALAXIES) {
			macc = parameters.kappa_agn * 0.9375 * PI * G_cgs * M_Atomic_g * mu_Primordial * Lcool * (mbh/MSOLAR_g);
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

	return parameters.nu_smbh * (macc/MACCRETION_cgs_simu) * std::pow(c_light_cm,2) / std::pow(10,40);
}

} // namespace shark
