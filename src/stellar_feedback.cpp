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

namespace shark {

StellarFeedbackParameters::StellarFeedbackParameters(const std::string &filename) :
	Options(filename),
	eta_cc(0),
	e_sn(0),
	epsilon_cc(0),
	beta(0)
{
	load("stellar_feedback.eta_cc", eta_cc);
	load("stellar_feedback.e_sn", e_sn);
	load("stellar_feedback.epsilon_cc", epsilon_cc);
	load("stellar_feedback.beta", beta);
}


StellarFeedback::StellarFeedback(StellarFeedbackParameters parameters) :
	parameters(parameters)
{
	// no-op
}

double StellarFeedback::outflow_rate(double sfr, double v) {
	return parameters.epsilon_cc * parameters.e_sn /
	       std::pow(v, parameters.beta) * parameters.eta_cc * sfr;
}

}  // namespace shark