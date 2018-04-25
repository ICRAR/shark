//
// Stellar feedback classes
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


#ifndef SHARK_STELLAR_FEEDBACK_H_
#define SHARK_STELLAR_FEEDBACK_H_

#include "options.h"

namespace shark {

class StellarFeedbackParameters {

public:
	StellarFeedbackParameters(const Options &options);

	enum StellarFeedbackModel {
		FIRE = 0,
		GUO11,
		LACEY16,
		LAGOS13,
		LAGOS13Trunc,
		LACEY16FIRE
	};

	/**
	 * Note that not all parameters will be used by the stellar feedback model, and which ones and how many are used
	 * depends on the model adopted.
	 * Parameter:
	 * eps_halo: constant in outflow rate from the halo.
	 * vkin_sn: kinetic velocity of SNe remnant.
	 * beta_disk: power-law index of velocity dependence on outflows from the galaxy.
	 * beta_halo: power-law index of velocity dependence on outflows from the halo.
	 * v_sn: normalization of velocity in SNe feedback/
	 * redshift_power: redshift dependence of normalization velocity.
	 * eps_disk: constant in outflow rate from the galaxy.
	 * e_sn: energy released by a single SNe that couples to the ISM.
	 * eta_cc: number of SNe per 1 solar mass of mass formed.
	 * galaxy_scaling: whether we scale SNe outflow rate with the halo or galaxy velocity.
	 * model: adopted SNe feedback model.
	 */
	double eps_halo;
	double vkin_sn;
	double beta_disk;
	double beta_halo;
	double v_sn;
	double redshift_power;
	double eps_disk;
	double e_sn;
	double eta_cc;
	bool galaxy_scaling;
	bool radial_feedback;

	StellarFeedbackModel model;

};


class StellarFeedback {

public:
	StellarFeedback(StellarFeedbackParameters parameters);

	void outflow_rate(double sfr, double vsubh, double vgal, double z, double &b1, double &b2, double &b_1, double &bj_2);

private:
	StellarFeedbackParameters parameters;
};

}  // namespace shark

#endif // SHARK_STELLAR_FEEDBACK_H_
