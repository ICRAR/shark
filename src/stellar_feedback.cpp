//
// Stellar feedback classes implementation
//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2017
// Copyright by UWA (in the framework of the ICRAR)
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307  USA
//

#include <cmath>

#include "stellar_feedback.h"
#include "numerical_constants.h"

namespace shark {

StellarFeedbackParameters::StellarFeedbackParameters(const Options &options) :
	eps_halo(0),
	vkin_sn(0),
	beta_disk(0),
	beta_halo(0),
	v_sn(0),
	redshift_power(0),
	eps_disk(0),
	eta_cc(0),
	e_sn(0)
{

	double epsilon_cc = 0, energy=0;

	options.load("stellar_feedback.eps_halo", eps_halo, true);
	options.load("stellar_feedback.vkin_sn", vkin_sn, true);
	options.load("stellar_feedback.beta_disk", beta_disk, true);
	options.load("stellar_feedback.beta_disk", beta_halo, true);
	options.load("stellar_feedback.v_sn", v_sn);
	options.load("stellar_feedback.redshift_power", redshift_power);
	options.load("stellar_feedback.eps_disk",eps_disk);
	options.load("stellar_feedback.e_sn",energy);
	options.load("stellar_feedback.eta_cc",eta_cc);
	options.load("stellar_feedback.epsilon_cc",epsilon_cc);

	//convert energy of SNe into Msun (km/s)^2
	e_sn = epsilon_cc * energy *std::pow(constants::MSOLAR_g, -1.0) * std::pow(constants::KILO, -2.0);
}


StellarFeedback::StellarFeedback(StellarFeedbackParameters parameters) :
	parameters(parameters)
{
	// no-op
}

void StellarFeedback::outflow_rate(double sfr, double v, double z, double *b1, double *b2) {
	/*
	 * TODO: add here other models for the outflow rate (e.g. GALFORM models and momentum driven models).
	 */

	*b1 = 0;
	*b2 = 0;

	if(sfr <= 0 || v <= 0){
		return;
	}

	// TODO: add as a different model option
	double power_index = parameters.beta_disk;
	if(v > parameters.v_sn){
		power_index = 1;
	}

	double const_sn =  std::pow((1+z),parameters.redshift_power) * std::pow(parameters.v_sn/v,power_index);
	//4/3 * parameters.eta_cc * parameters.e_sn * std::pow((1+z),parameters.redshift_power) / std::pow(v,2.0) * std::pow(parameters.v_sn/v,power_index);


	*b1 = parameters.eps_disk * const_sn;

	double vsn = 1.9 * std::pow(v,1.1);

	double eps_halo = parameters.eps_halo * const_sn *  0.5 * std::pow(vsn,2.0);

	double energ_halo = 0.5 * std::pow(v,2.0);

	double mreheat = *b1 * sfr;

	double mejected = eps_halo / energ_halo * sfr - mreheat;


	if(mejected > 0) {
		*b2 = mejected/sfr;
		if(*b2 > *b1){
			*b2 = *b1;
		}
	}
	else{
		*b1 = eps_halo / energ_halo;
	}

}

}  // namespace shark
