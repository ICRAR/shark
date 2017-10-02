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
	eta_cc(0),
	e_sn(0),
	epsilon_cc(0),
	beta(0),
	v_sn(0)
{
	options.load("stellar_feedback.eta_cc", eta_cc, true);
	options.load("stellar_feedback.e_sn", e_sn, true);
	options.load("stellar_feedback.epsilon_cc", epsilon_cc, true);
	options.load("stellar_feedback.beta", beta, true);
	options.load("stellar_feedback.v_sn", v_sn);

	//convert energy of SNe into code units.
	e_sn = e_sn *std::pow(constants::MSOLAR_g, -1.0) * std::pow(constants::MPC2CM, -2.0) * std::pow(constants::GYR2S, 2.0);
}


StellarFeedback::StellarFeedback(StellarFeedbackParameters parameters) :
	parameters(parameters)
{
	// no-op
}

double StellarFeedback::outflow_rate(double sfr, double v) {
	/*
	 * TODO: add here other models for the outflow rate (e.g. GALFORM models and momentum driven models).
	 */
	if(sfr <= 0 || v <= 0){
		return 0.0;
	}

	// TODO: add as a different model option
//	double beta  = parameters.epsilon_cc * parameters.e_sn /
//		       std::pow(v, parameters.beta) * parameters.eta_cc * sfr;

	double beta_galform = std::pow(parameters.v_sn/v, parameters.beta);

	return beta_galform;
}

}  // namespace shark
