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
 * Stellar feedback classes
 */

#ifndef SHARK_STELLAR_FEEDBACK_H_
#define SHARK_STELLAR_FEEDBACK_H_

#include "options.h"

namespace shark {

class StellarFeedbackParameters {

public:
	explicit StellarFeedbackParameters(const Options &options);

	enum StellarFeedbackModel {
		MURATOV15 = 0,
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
	 * min_beta: minimum mass loading factor for the galaxy. 
	 */
	double eps_halo = 1;
	double vkin_sn = 0;
	double beta_disk = 0;
	double beta_halo = 0;
	double v_sn = 0;
	double redshift_power = 0;
	double eps_disk = 1;
	double e_sn = 0;
	double eta_cc = 0;
	double min_beta = 0;
	bool galaxy_scaling = false;
	bool radial_feedback = false;

	StellarFeedbackModel model = MURATOV15;

};


class StellarFeedback {

public:
	explicit StellarFeedback(StellarFeedbackParameters parameters);

	void outflow_rate(double sfr, double vsubh, double vgal, double z, double &b1, double &b2, double &bj_1, double &bj_2);

private:
	StellarFeedbackParameters parameters;
};

}  // namespace shark

#endif // SHARK_STELLAR_FEEDBACK_H_
