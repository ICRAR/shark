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
#include "galaxy.h"
#include "halo.h"
#include "subhalo.h"
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

	// relevant for Croton16 and Lagos 22 models.
	options.load("agn_feedback.kappa_agn", kappa_agn);
	options.load("agn_feedback.accretion_eff_cooling", nu_smbh);

	// relevant for Lagos 22 model.
	options.load("agn_feedback.kappa_radio", kappa_radio);
	options.load("agn_feedback.hot_halo_threshold", hot_halo_threshold);
	options.load("agn_feedback.spin_model", spin_model);
	options.load("agn_feedback.eta_superedd", eta_superedd);

	// parameters below are relevant for the Griffin 2020's spin model.
	options.load("agn_feedback.alpha_adaf", alpha_adaf);
	options.load("agn_feedback.alpha_td", alpha_td);
	options.load("agn_feedback.delta_adaf", delta_adaf);
	options.load("agn_feedback.mdotcrit_adaf", mdotcrit_adaf);
	options.load("agn_feedback.accretion_disk_model", accretion_disk_model);
	options.load("agn_feedback.loop_limit_accretion", loop_limit_accretion);

	auto beta = 1 - alpha_adaf / 0.55;
	low_accretion_adaf = 0.001 * (delta_adaf / 0.0005) * (1 - beta) / beta * std::pow(alpha_adaf, 2.0);
	constant_lowlum_adaf = (delta_adaf / 0.0005) * (1 - beta) / 0.5 * 6;
	constant_highlum_adaf = beta / 0.5 / std::pow(alpha_adaf, 2) * 6;
	nu2_nu1 = std::pow(alpha_td, 2.0);

	// control QSO feedback.
	options.load("agn_feedback.qso_feedback", qso_feedback);
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
	else if (lvalue == "lagos23") {
		return AGNFeedbackParameters::LAGOS23;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value << ". Supported values are bower06, croton16 and lagos23";
	throw invalid_option(os.str());
}

template <>
AGNFeedbackParameters::SpinModel
Options::get<AGNFeedbackParameters::SpinModel>(const std::string &name, const std::string &value) const {
	auto lvalue = lower(value);
	if (lvalue == "volonteri07") {
		return AGNFeedbackParameters::VOLONTERI07;
	}
	else if (lvalue == "griffin19") {
		return AGNFeedbackParameters::GRIFFIN19;
	}
	else if (lvalue == "constant") {
		return AGNFeedbackParameters::CONSTANTSPIN;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value << ". Supported values are constant, volonteri07, and griffin19";
	throw invalid_option(os.str());
}

template <>
AGNFeedbackParameters::AccretionDiskModel
Options::get<AGNFeedbackParameters::AccretionDiskModel>(const std::string &name, const std::string &value) const {
	auto lvalue = lower(value);
	if (lvalue == "warpeddisk") {
		return AGNFeedbackParameters::WARPEDDISK;
	}
	else if (lvalue == "selfgravitydisk") {
		return AGNFeedbackParameters::SELFGRAVITYDISK;
	}
	else if (lvalue == "prolonged") {
		return AGNFeedbackParameters::PROLONGED;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value << ". Supported values are warpeddisk, selfgravitydisk or prolonged";
	throw invalid_option(os.str());
}

AGNFeedback::AGNFeedback(const AGNFeedbackParameters &parameters, CosmologyPtr cosmology, RecyclingParameters recycle_params, ExecutionParameters exec_params) :
	parameters(parameters),
	cosmology(std::move(cosmology)),
	recycle_params(std::move(recycle_params)),
	exec_params(std::move(exec_params)),
	distribution(-1, 1)
{
	// no-op
}

void AGNFeedback::plant_seed_smbh(Subhalo &subhalo){

	auto central = subhalo.central_galaxy();

	if(central){
		float mvir = 0;
		if(subhalo.subhalo_type == Subhalo::CENTRAL){
			mvir = subhalo.host_halo->Mvir;
		}
		else{
			mvir = subhalo.Mvir_infall;
		}
		if (mvir > parameters.mhalo_seed && central->smbh.mass ==0){
				central->smbh.mass = parameters.mseed;
				central->smbh.spin = 0;
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

double AGNFeedback::accretion_rate_hothalo_smbh(double Lcool, double tacc, double fhot, double vvir, Galaxy &galaxy) {

	/**
	 * Function calculates the accretion rate onto the central black hole based on a cooling luminosity.
	 * Inputs:
	 * Lcool: cooling luminosity in units of 10^40 erg/s.
	 * mbh: mass of central supermassive black hole.
	 */

	using namespace constants;

	auto &smbh = galaxy.smbh;

	if (Lcool > 0 && Lcool < MAXLUM) {

		double macc = 0;
		if (parameters.model == AGNFeedbackParameters::BOWER06) {
			macc = Lcool * 1e40 / std::pow(c_light_cm , 2.0) / parameters.accretion_eff_cooling;
		}
		else if (parameters.model == AGNFeedbackParameters::CROTON16) {
			macc = parameters.kappa_agn * 0.9375 * PI * G_cgs * M_Atomic_g * mu_Primordial * Lcool * 1e40 * (smbh.mass * MSOLAR_g);
		}
		else if (parameters.model == AGNFeedbackParameters::LAGOS23) {
			// here we adopt Croton et al. (2006)
			macc = parameters.kappa_agn * (smbh.mass / 1e8) * (fhot / 0.1) * std::pow( vvir / 200.0, 3.0);
		}

		// calculate new spin if necessary
		if(parameters.model == AGNFeedbackParameters::LAGOS23){
			// in this case compute spin
			if(parameters.spin_model == AGNFeedbackParameters::VOLONTERI07){
				volonteri07_spin(smbh);
			}
			else if (parameters.spin_model == AGNFeedbackParameters::GRIFFIN19){
				auto deltam = macc * MACCRETION_cgs_simu * tacc;
				griffin19_spinup_accretion(deltam, tacc, galaxy);
			}
		}

		return macc * MACCRETION_cgs_simu; //accretion rate in units of Msun/Gyr.
	}
	else{
		return 0;
	}

}

double AGNFeedback::accretion_rate_ratio(double macc, double mBH){

	//return the ratio between the physical and Eddington mass accretion rates.
	using namespace constants;

	if(macc > 0){
		double LEdd = eddington_luminosity(mBH);
		double M_dot_Edd = 1e40 * LEdd / (parameters.nu_smbh * std::pow(c_light_cm,2.0));
		double m_dot = (macc/MACCRETION_cgs_simu) / M_dot_Edd;
		return m_dot;
	}
	else{
		return 0;
	}
}

double AGNFeedback::agn_bolometric_luminosity(const BlackHole &smbh, bool starburst){

	//return bolometric luminosity in units of 10^40 erg/s.
	//Follows equations from Griffin et al. (2020).

	using namespace constants;

	auto mBH = cosmology->comoving_to_physical_mass(smbh.mass);

	auto macc = cosmology->comoving_to_physical_mass(smbh.macc_hh);

	if(parameters.model == AGNFeedbackParameters::LAGOS23 || starburst){
		//In this case also sum the starburst accretion rate
		macc += cosmology->comoving_to_physical_mass(smbh.macc_sb);
	}

	// assume constant radiation efficiency unless this model is Lagos22
	double Lbol = 0;
	if (parameters.model == AGNFeedbackParameters::LAGOS23) {
		double LEdd = eddington_luminosity(mBH);
		double m_dot_norm = accretion_rate_ratio(macc,mBH);
		
		auto eff = efficiency_luminosity_agn(smbh.spin);

		if(m_dot_norm >= parameters.mdotcrit_adaf){
			Lbol = eff[0] * (macc / MACCRETION_cgs_simu) * std::pow(c_light_cm,2.0) / 1e40;
			if(m_dot_norm > parameters.eta_superedd * (0.1 / eff[0])){
				Lbol = parameters.eta_superedd * (1.0 + std::log(m_dot_norm / parameters.eta_superedd * eff[0] / 0.1)) * LEdd;
			}
		}
		else{
			if(m_dot_norm > parameters.low_accretion_adaf){
				Lbol = 0.2 * eff[0] * (macc / MACCRETION_cgs_simu)  * std::pow(c_light_cm, 2.0) * m_dot_norm * parameters.constant_highlum_adaf / eff[1] / 1e40;
			}
			else{
				Lbol = 0.0002 * eff[0] * (macc / MACCRETION_cgs_simu) * std::pow(c_light_cm, 2.0) * parameters.constant_lowlum_adaf / eff[1] / 1e40;
			}
		}
	}
	else if (parameters.model == AGNFeedbackParameters::CROTON16) {
		Lbol = parameters.nu_smbh * (macc / MACCRETION_cgs_simu) * std::pow(c_light_cm, 2.0) / 1e40;
	}

	return Lbol;
}

double AGNFeedback::agn_mechanical_luminosity(const BlackHole &smbh){
	
	//return mechanical luminosity in units of 10^40 erg/s.
	//Note that here we use the radio luminosities as described in Griffin et al. (2019; arxiv: 1912.09490) Eqs. 11 and 12.
	using namespace constants;
	
	auto macc = cosmology->comoving_to_physical_mass(smbh.macc_hh + smbh.macc_sb);
	auto mBH = cosmology->comoving_to_physical_mass(smbh.mass);

	double m_dotdiv0p01 = accretion_rate_ratio(macc, mBH) / 0.01;
	double Lmech = 0;

	if(m_dotdiv0p01 >= 1.0){
		Lmech = 2.5e3 * std::pow(mBH/1e9, 1.1) * std::pow(m_dotdiv0p01, 1.2) * std::pow(smbh.spin, 2);
	}
	else{
		Lmech = 2e5 * (mBH/1e9) * (m_dotdiv0p01)  * std::pow(smbh.spin, 2);
	}

	return Lmech;
}

double AGNFeedback::accretion_rate_hothalo_smbh_limit(double mheatrate, double vvir, const BlackHole &smbh){

	using namespace constants;

	double Lbol = mheatrate * (0.5*std::pow(vvir*KM2CM,2.0)) / MACCRETION_cgs_simu / 1e40;

	auto eff = efficiency_luminosity_agn(smbh.spin);
	double macc = Lbol / eff[0] / std::pow(c_light_cm,2.0) * 1e40 * MACCRETION_cgs_simu;

	return macc;

}

double AGNFeedback::smbh_growth_starburst(double mgas, double vvir, double tacc, Galaxy &galaxy){

	double m = 0;

	auto &smbh = galaxy.smbh;

	float TOL = 0.99;

	// Grow supermassive BH only if above the black hole seed (avoid formation of low BH masses without seeds).
        if(smbh.mass >= TOL * parameters.mseed){
		if(mgas > 0){
			m =  parameters.f_smbh * mgas / (1 + std::pow(parameters.v_smbh/vvir, 2.0));
		}
        
		// calculate new spin if necessary
		if(parameters.model == AGNFeedbackParameters::LAGOS23){
			// in this case compute spin
			if(parameters.spin_model == AGNFeedbackParameters::VOLONTERI07){
				volonteri07_spin(smbh);
			}
			else if (parameters.spin_model == AGNFeedbackParameters::GRIFFIN19){
				griffin19_spinup_accretion(m, tacc, galaxy);
			}
		}
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

double AGNFeedback::qso_outflow_velocity(double Lbol, double mbh, double zgas, double mgas, double mbulge, double rbulge){

	double vout  = 320.0 * std::pow(Lbol/(1e7 * constants::LSOLAR), 0.5) * std::pow(zgas/recycle_params.zsun, 0.25) * std::pow(mgas, -0.25);

	return vout;

}

void AGNFeedback::qso_outflow_rate(double mgas, const BlackHole &smbh, double zgas, double vcirc,
		double sfr, double mbulge, double rbulge, double &beta_halo, double &beta_ejec){

	double macc = cosmology->comoving_to_physical_mass(smbh.macc_sb + smbh.macc_hh);
	double mBH = cosmology->comoving_to_physical_mass(smbh.mass);

	// QSO feedback only acts if the accretion rate is >0, BH mass is > 0 and QSO feedback is activated by the user.
	if(macc > 0 && mBH > 0 && sfr > 0 && mgas > 0 && parameters.qso_feedback){
		double Lbol = agn_bolometric_luminosity(smbh, true);
		double Lcrit = qso_critical_luminosity(mgas, mbulge, rbulge);

		// check if bolometric luminosity is larger than the critical luminosity and the gas mass in the bulge is positive. The latter is not always the case becaus equations are solved
		// numerically and hence negative solutions are in principle possible.
		if(Lbol > Lcrit){

			double tsalp = salpeter_timescale(Lbol, mBH);
			double vout = qso_outflow_velocity(Lbol, mBH, zgas, mgas, mbulge, rbulge);

			double mout_rate= mgas/tsalp;

			double mejec_rate = (parameters.epsilon_qso * std::pow(vout/vcirc, 2.0) - 1) * mout_rate;

			// Apply boundary conditions to outflow and ejection rates
			if(mout_rate <  0 || std::isnan(mout_rate)){
				mout_rate = 0;
			}
			if(mejec_rate <  0 || std::isnan(mejec_rate)){
				mejec_rate = 0;
			}

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

void AGNFeedback::volonteri07_spin(BlackHole &smbh){

	if(smbh.mass > 0){
		auto logmbh = std::log10(smbh.mass);
		smbh.spin = 0.305 * logmbh - 1.7475;
	}
	else{
		smbh.spin = 0;
	}

	// putting physical limits to spin
	if(smbh.spin > 1){
		smbh.spin = 1;
	}
	else if (smbh.spin < -1){
		smbh.spin = -1;
	}

}

void AGNFeedback::griffin19_spinup_accretion(double delta_mbh, double tau_acc, Galaxy &galaxy){

	/*Function computes the black hole spin resulting from the black hole accretion during starbursts.
	Model follows that published by Griffin et al. (2020).*/

	//if no accretion takes place then return.
	if(delta_mbh <= 0 || galaxy.smbh.mass == 0){
		return;
	}

	double m_in = 0;
	int n_accretion_chunks = 10;
	int loop_inner = 0;
	int loop_outer = 0;
	double eff = 0;
	double mdot_norm = 0;
	double mfin = 0;
	double delta_m = 0;
	double angular_momentum_ratio = 0;
	double R_angm = 0;
	bool retrograde_accretion;

	// define supermassive black hole
	auto &smbh = galaxy.smbh;

	// Compute quantities in physical units before proceeding.
	auto m_out = cosmology->comoving_to_physical_mass(delta_mbh);
	auto mbh = cosmology->comoving_to_physical_mass(smbh.mass);
	auto mdot = cosmology->comoving_to_physical_mass(delta_mbh/tau_acc);
	auto spin = smbh.spin;

	if(parameters.accretion_disk_model == AGNFeedbackParameters::PROLONGED){
		auto efficiency = efficiency_luminosity_agn(spin);
		eff = efficiency[0];
		auto r_lso = efficiency[1];
		spin = final_spin(mbh, mbh + m_out, r_lso);
	}
	else{
		do {
			// Initialise angle between accretion disk and total one (assume random orientations)
			auto theta = angle_acc_disk(galaxy);
			auto cos_theta_i = std::cos(theta);

			float r_lso = 0;

			if(mbh > 0){
				auto efficiency = efficiency_luminosity_agn(spin);
				eff = efficiency[0];
				r_lso = efficiency[1];
				mdot_norm = accretion_rate_ratio(mdot, mbh);

				// Self-gravity radius based on King, Pringle, Hofmann 2008 equation 10 but with different constant at the front.
				auto r_sg = 4790.0 * std::pow(parameters.alpha_td, 0.5185) * 
						std::pow(mdot_norm, -0.2962) * 
						std::pow(mbh/1e8, -0.9629); //This is in units of 2r_G
				double m_sg = 0;

				if(mdot_norm < parameters.mdotcrit_adaf){
					// ADAF disk.
					m_sg = mbh;
				}
				else{
					/* Self-gravity mass based on King, Pringle, Hofmann 2008 equation 7
					but with different constant at the front.
					m_sg = M_BH * (H/R) is the expression used.*/
					m_sg = 1.35 * std::pow(parameters.alpha_td, -0.8) * 
						std::pow(mdot_norm, 0.6) * 
						std::pow(mbh / 1e8, 2.2) * 
						std::pow(r_sg, 1.4); // In Msun
				}

				m_in = std::min(m_sg , delta_mbh);

				// Compute new spin based on accretion onto the black hole.
				do{
					if(spin == 0){
						/*If spin equals zero, we break up the spinup into ten smaller pieces to make the
						accretion efficiency more realistic.*/
						auto ichunk = n_accretion_chunks;
						do{
							efficiency = efficiency_luminosity_agn(spin);
							eff = efficiency[0];
							r_lso = efficiency[1];
							mfin = mbh + (1 - eff) * (m_in / n_accretion_chunks);
							spin = final_spin(mbh, mfin, r_lso);
							mbh = mbh + (1- eff) * (m_in / n_accretion_chunks);
							ichunk --;
						}
						while (ichunk > 0);
						m_in = -1; //To stop the loop
					}
					else{
						/* The disk will be inclined with respect to the BH equatorial plane, thus, the BH
						accretion disk system will be subject to Lens-Thirring precession. In this case it is
						necessary to follow the evolution of the misalignment angle phi since it is expected to
						affect the final spin acquired by the BH.*/
						efficiency = efficiency_luminosity_agn(spin);
						eff = efficiency[0];
						r_lso = efficiency[1];

						// Warp radius based on Volonteri 2007 equation 2, with a slightly different constant at the front.
						auto r_warp = 3410 * std::pow( std::abs(spin), 0.625) *
								std::pow(mbh/1e8, 0.125) *
								std::pow(mdot_norm, -0.25) *
								std::pow(parameters.nu2_nu1, -0.625) *
								std::pow(parameters.alpha_td, -0.5); // Units of 2r_G

						// m_warp equal to Sigma R**2, as based on King, Pringle, Hofmann 2008 equation 7, but with a different constant at the front.
						auto m_warp = 1.35 * std::pow(parameters.alpha_td, -0.8) *
								std::pow(mdot_norm, 0.6) *
								std::pow(mbh/1e8, 2.2) *
								std::pow(r_warp, 1.4); //Msun

						if(parameters.accretion_disk_model == AGNFeedbackParameters::WARPEDDISK){
							delta_m = m_warp;

							// check that warped mass isn't larger than the available accreting mass.
							if(delta_m > m_in){
								delta_m = m_in;
								r_warp = r_sg;
							}
							R_angm = r_warp;

						}
						else if(parameters.accretion_disk_model == AGNFeedbackParameters::SELFGRAVITYDISK){
							delta_m = m_in;
							R_angm = r_sg;
						}

						mfin = mbh + (1 - eff) * delta_m;

						//compute angular momentum ratio
						angular_momentum_ratio = (delta_m / (constants::SQRT2 * mbh * std::abs(spin))) * std::sqrt(R_angm);

						auto J_SMBH = std::abs(spin) * std::pow(mbh, 2) * constants::G / constants::c_light_km;
						auto J_disk = 2 * angular_momentum_ratio * J_SMBH;

						// Angle evolution
						auto cos_theta_f = (J_disk + J_SMBH * cos_theta_i)/
								std::sqrt(std::pow(J_SMBH, 2) + std::pow(J_disk, 2) + 2 * J_SMBH * J_disk * cos_theta_i);

						// Making sure cos_theta stays within the correct limits.
						if(cos_theta_f > 1){
							cos_theta_f = 1;
						}
						else if(cos_theta_f < -1){
							cos_theta_f = -1;
						}
						theta = std::acos(cos_theta_f);

						if (cos_theta_f > -1 * angular_momentum_ratio){
							// Prograde accretion
							retrograde_accretion = false;
							efficiency = efficiency_luminosity_agn(spin);
						}
						else{
							// Retrograde accretion
							retrograde_accretion = true;
							efficiency = efficiency_luminosity_agn(-1 * spin);
						}
						eff = efficiency[0];
						r_lso = efficiency[1];
						spin = final_spin(mbh, mfin, r_lso);

						spin = std::abs(spin); // Spin set to magnitude of spin.
						// Set the spin as less than 0 if the last accretion epsiode is retrograde
						if (retrograde_accretion){
							spin = -1 * spin;
						}

						mbh = mbh + (1 - eff) * delta_m; // Not all mass accreted on, some lost as radiation.
						m_in = m_in - delta_m;
						loop_inner += 1;
					}
				}
				while(m_in > 0 && loop_inner < parameters.loop_limit_accretion);

				//check is there is accretion disk mass left (which can happen if
				//we stopped the loop due to the number of iterations hiting the limit. In that case, modify the BH mass accordingly.
				if(m_in > 0){
					mbh = mbh + (1 - eff) * m_in; // Not all mass accreted on, some lost as radiation.
				}

				// End of spin-evolution calculation
				if(parameters.accretion_disk_model == AGNFeedbackParameters::WARPEDDISK){
					m_out = -1;
				}
				else if (parameters.accretion_disk_model == AGNFeedbackParameters::SELFGRAVITYDISK){
					if (m_out > m_sg){
						m_out -= m_sg;
					}
					else{
						m_out = -1;
					}
				}
				loop_outer += 1;
			}
		}
		while(m_out > 0  && loop_outer < parameters.loop_limit_accretion);
	}

	smbh.spin = spin;

	// Check for undefined cases.
	if(spin < -1 || spin > 1 || std::isnan(spin) || std::isinf(spin)){
		std::ostringstream os;
		os << "SMBH in accretion routine has spin not well defined";
		throw invalid_data(os.str());
	}

}

void AGNFeedback::griffin19_spinup_mergers(BlackHole &smbh_primary, const BlackHole &smbh_secondary, const Galaxy &galaxy){

        using namespace constants;

	double t0 = -2.686,
			t2 = -3.454,
			t3 = 2.353,
			s4 = -0.129,
			s5 = -0.384;

	auto m1 = smbh_primary.mass;
	auto m2 = smbh_secondary.mass;
	auto s1 = smbh_primary.spin;
	auto s2 = smbh_secondary.spin;

        // Apply Rezzolla et al. (2008) if at least one of the BHs has a non-zero spin.
	if ( m1 > 0 && m2 > 0 && (s1 !=0 || s2 != 0)){

		// The subscript 1 refers to the spin of black hole 1
		// We are using spherical polar coordinates
		std::default_random_engine generator(exec_params.get_seed(galaxy));
		auto theta_1 = angle_acc_disk(generator);
		auto theta_2 = angle_acc_disk(generator);

		auto cos_theta_1 = std::cos(theta_1);
		auto sin_theta_1 = std::sin(theta_1);

		auto cos_theta_2 = std::cos(theta_2);
		auto sin_theta_2 = std::sin(theta_2);

		auto phi_2 = phi_acc_disk(generator);
		auto cos_phi_2 = std::cos(phi_2);

		// alpha is the cos of the angle between spin 1 and spin 2
		auto alpha = std::abs(sin_theta_1 * sin_theta_2 * cos_phi_2
				+ cos_theta_1 * cos_theta_2);

		//beta is the cos of the angle between the angular momentum of the orbit and spin 1
		auto beta = std::abs(cos_theta_1);
		// gama is the cos of the angle between the angular momentum of the orbit and spin 2
		auto gama = std::abs(cos_theta_2);

		auto symm_mr = m1 * m2 / std::pow(m1 + m2, 2);

		auto mr = m2 / m1;
		auto mrO2 = std::pow(mr, 2);
		auto mrO4 = std::pow(mr, 4);

		if(mr > 1){
			mr = 1 / mr;
			s1 = smbh_secondary.spin;
			s2 = smbh_primary.spin;
		}

		s1 = std::abs(s1);
		s2 = std::abs(s2);

		auto ang_mom = s4 / std::pow(1 + mrO2, 2) * (std::pow(s1, 2) + std::pow(s2, 2) * mrO4 + 2 * s1 * s2 * mrO2 * alpha) +
				(s5 * symm_mr + t0 + 2) / (1 + mrO2) * (s1 * beta + s2 * mrO2 * gama) +
				2 * constants::SQRT3 + t2 * symm_mr + t3 * std::pow(symm_mr, 2);
		ang_mom = std::abs(ang_mom);

		auto param_new_spin = std::pow(s1, 2) + std::pow(s2, 2) * mrO4 + 2 * s1 * s2 * mrO2 *alpha +
				2 * (s2 * beta + s2 * mrO2 * gama) * ang_mom * mr + std::pow(ang_mom, 2) * mrO2;
		auto new_spin =  1.0 / std::pow(1.0 + mr, 2.0) * std::sqrt(param_new_spin);

		if(new_spin > 0.998){
			new_spin = 0.998;
		}

		smbh_primary.spin = new_spin;

		// Check for undefined cases.
		if(new_spin < -1 || new_spin > 1 || std::isnan(new_spin)){
			std::ostringstream os;
			os << " SMBH in merging black holes has spin not well defined";
			throw invalid_data(os.str());
		}
	}
	else if (m1 > 0 && m2 > 0 && s1 ==0 && s2 == 0){
	// when both BHs are spinless, then apply Berti & Volonteri (2008)
		double q = m2 / m1;
		if ( q > 1) q = 1.0 / q;
		
		auto new_spin = 2 * SQRT3 * q / std::pow( 1 + q, 2.0) - 2.029 * std::pow(q, 2.0) / std::pow( 1 + q, 4.0);

		// make sure it doesn't go negative or above 1.
		if ( new_spin < 0) new_spin = 0;
		if ( new_spin > 1) new_spin = 1;

		smbh_primary.spin = new_spin;
	}
	else{
		if(m1 == 0){
			smbh_primary.spin = 0;
		}
	}

}

double AGNFeedback::angle_acc_disk(const Galaxy &galaxy){
	std::default_random_engine generator(exec_params.get_seed(galaxy));
	return angle_acc_disk(generator);
}

double AGNFeedback::angle_acc_disk(std::default_random_engine &generator) {
	return std::acos(distribution(generator));
}

double AGNFeedback::phi_acc_disk(std::default_random_engine &generator){
	return constants::PI2 * (1 + distribution(generator));

}

std::vector<float> AGNFeedback::efficiency_luminosity_agn(float spin){

	std::vector<float> efficiency;

	if(parameters.spin_model == AGNFeedbackParameters::CONSTANTSPIN || parameters.model == AGNFeedbackParameters::CROTON16){
		efficiency.push_back(parameters.nu_smbh);
		efficiency.push_back(0);
		return efficiency;
	}

	double eff, r_lso;

	// Calculate the radiation efficiency in the thin disk approximation

	auto a = std::abs(spin);
	auto a2 = std::pow(a,2);

	auto z1 = 1 + std::pow(1 - a2, 0.333) * (std::pow(1 + a, 0.333) + std::pow(1 - a, 0.333));
	auto z2 = std::sqrt(3 * a2 + std::pow(z1, 2));

	if(a >= 0){
		r_lso = 3 + z2 - std::sqrt((3 - z1) * (3 + z1 + 2 * z2));
	}
	else{
		r_lso = 3 + z2 + std::sqrt((3 - z1) * (3 + z1 + 2 * z2));
	}

	eff = 1 - std::sqrt(1 - 2 / (3 * r_lso));

	if(eff < 0){
		eff = 0.07;
	}
	else if(eff > 0.5){
		eff = 0.5;
	}

	efficiency.push_back(eff);
	efficiency.push_back(r_lso);

	return efficiency;

}

double AGNFeedback::final_spin(const double mbh, const double mfin, const double r_lso){

	double spin = 0;

	auto mrat = mbh / mfin;

	if(mrat > std::sqrt(2 / (3 * r_lso))){
		spin = 0.333 * std::sqrt(r_lso) * mrat * (4 - std::sqrt(3 * r_lso * std::pow(mrat, 2) - 2 ));
		spin = std::min(spin, 0.998);
	}
	else{
		spin = 0.998;
	}

	if (spin < 0 && spin > -0.01) {
		spin = -0.01;
	}
	else if (spin > 0 && spin < 0.01) {
		spin = 0.01;
	}

	return spin;
}


} // namespace shark
