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
 *
 * Stellar feedback classes implementation
 */

#include <cmath>

#include "numerical_constants.h"
#include "stellar_feedback.h"
#include "utils.h"

namespace shark {

StellarFeedbackParameters::StellarFeedbackParameters(const Options &options)
{

	double epsilon_cc = 0, energy=0;

	options.load("stellar_feedback.model", model, true);
	options.load("stellar_feedback.galaxy_scaling", galaxy_scaling);
	options.load("stellar_feedback.radial_feedback", radial_feedback);

	// The parameters below *must* be specified.
	options.load("stellar_feedback.beta_disk", beta_disk, true);
	options.load("stellar_feedback.v_sn", v_sn, true);
	options.load("stellar_feedback.eps_halo", eps_halo);
	options.load("stellar_feedback.eps_disk",eps_disk);
	options.load("stellar_feedback.redshift_power", redshift_power);

	// The parameters below don't need to be specified.
	options.load("stellar_feedback.vkin_sn", vkin_sn);
	options.load("stellar_feedback.e_sn",energy);
	options.load("stellar_feedback.eta_cc",eta_cc);
	options.load("stellar_feedback.epsilon_cc",epsilon_cc);
	options.load("stellar_feedback.beta_halo", beta_halo);
	options.load("stellar_feedback.min_beta", min_beta);

	//convert energy of SNe into Msun (km/s)^2
	e_sn = epsilon_cc * energy *std::pow(constants::MSOLAR_g, -1.0) * std::pow(constants::KILO, -2.0);
}

template <>
StellarFeedbackParameters::StellarFeedbackModel
Options::get<StellarFeedbackParameters::StellarFeedbackModel>(const std::string &name, const std::string &value) const {
	auto lvalue = lower(value);
	if (lvalue == "muratov15") {
		return StellarFeedbackParameters::MURATOV15;
	}
	else if (lvalue == "lacey16") {
		return StellarFeedbackParameters::LACEY16;
	}
	else if (lvalue == "guo11") {
		return StellarFeedbackParameters::GUO11;
	}
	else if (lvalue == "lagos13") {
		return StellarFeedbackParameters::LAGOS13;
	}
	else if (lvalue == "lagos13trunc") {
		return StellarFeedbackParameters::LAGOS13Trunc;
	}
	else if (lvalue == "lacey16reddep") {
		return StellarFeedbackParameters::LACEY16FIRE;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value << ". Supported values are muratov15, lacey16, guo11, lagos13, lagos13trunc and lacey16redrep";
	throw invalid_option(os.str());
}

StellarFeedback::StellarFeedback(StellarFeedbackParameters parameters) :
	parameters(parameters)
{
	// no-op
}

void StellarFeedback::outflow_rate(double sfr, double vsubh, double vgal, double z, double &b1, double &b2, double &bj_1, double &bj_2) {


	double v = vsubh;
	if(parameters.galaxy_scaling && vgal > 0){
		v  = vgal;
	}

	b1 = 0;
	b2 = 0;

	if(sfr <= 0 || v <= 0){
		return;
	}

	double vsn = 1.9 * std::pow(v,1.1);

	double power_index = parameters.beta_disk;
	double const_sn = 0;
	if (parameters.model == StellarFeedbackParameters::MURATOV15){

		if(v > parameters.v_sn){
			power_index = 1;
		}
		const_sn =  std::pow((1+z),parameters.redshift_power) * std::pow(parameters.v_sn/v,power_index);

	}
	else if (parameters.model == StellarFeedbackParameters::LAGOS13){

		double vhot = parameters.v_sn*std::pow(1+z,parameters.redshift_power);
		const_sn =  std::pow(vhot/v,power_index);
	}

	else if (parameters.model == StellarFeedbackParameters::LAGOS13Trunc){
		double vhot = parameters.v_sn*std::pow(1+z,parameters.redshift_power);

		if(v > parameters.v_sn){
			power_index = 1;
		}

		const_sn =  std::pow(vhot/v,power_index);
	}

	else if (parameters.model == StellarFeedbackParameters::LACEY16){

		const_sn = std::pow(parameters.v_sn/v,power_index);
	}
	else if (parameters.model == StellarFeedbackParameters::GUO11){

		const_sn = 0.5 + std::pow(parameters.v_sn/v,power_index);
	}
	else if (parameters.model == StellarFeedbackParameters::LACEY16FIRE){
		const_sn = std::pow((1+z),parameters.redshift_power) * std::pow(parameters.v_sn/v,power_index);
	}

	b1 = parameters.eps_disk * const_sn;

        if(b1 < parameters.min_beta){
		b1 = parameters.min_beta;
		const_sn = b1/parameters.eps_disk;
	}

	double eps_halo = parameters.eps_halo * const_sn *  0.5 * std::pow(vsn,2.0);

	double energ_halo = 0.5 * std::pow(v,2.0);

	double mreheat = b1 * sfr;


	double mejected = eps_halo / energ_halo * sfr - mreheat;

	if(mejected > 0) {
		b2 = mejected/sfr;
		if(b2 > b1){
			b2 = b1;
			b1 += constants::EPS3; //add a small number to b1 to make it strictly larger than b2.
		}
	}
	/*else{
		b1 = eps_halo / energ_halo;
	}*/

	// If no radial feedback is applied, then change in angular momentum reflects that of the mass.
	if(!parameters.radial_feedback){
		bj_1 = b1;
		bj_2 = b2;
	}
	else{
		//TODO: implement this part following Peter Mitchell's model for angular momentum.
	}

}

}  // namespace shark
