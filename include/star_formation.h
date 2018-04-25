/*
 * star_formation.h
 *
 *  Created on: 17May,2017
 *      Author: clagos
 */

//
// star formation classes
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


#ifndef INCLUDE_STAR_FORMATION_H_
#define INCLUDE_STAR_FORMATION_H_

#include <memory>

#include "components.h"
#include "cosmology.h"
#include "integrator.h"
#include "options.h"
#include "recycling.h"

namespace shark {

class StarFormationParameters {

public:
	StarFormationParameters(const Options &options);

	/**Parameters:
	 * model: SF model applied
	 * nu_sf: inverse of the molecular gas depletion timescale for normal star-forming disks [Gyr^-1].
	 * Po: pressure normalization in the Blitz & Rosolowsky (2006) SF law [K cm^3].
	 * beta_press: pressure power-law index in the Blitz & Rosolowsky (2006) SF law.
	 * Accuracy_SFeqs: accuracy desired when numerically solving the integration of the SFR surface density.
	 * gas_velocity_dispersion: gas velocity dispersion [km/s].
	 * boost_starburst: boosting factor applied to the SF law in starbursts.
	 * sigma_HI_crit: critical surface density below which gas is assumed to be ionised.
	 * clump_factor_KMT09: clumping factor of the ISM for the Krumholz+ models.
	 * sigma_crit_KMT09: critical gas surface density above which the SF law becomes superlinear in the KMT09 model.
	 * angular_momentum_transfer: boolean parameter indicating whether the user wants to trigger the calculation of angular momentum transfer within the disk.
	 *
	 */
	enum StarFormationModel {
		BR06 = 0,//!< BR06
		GD14,    //!< GD14
		K13,     //!< K13
		KMT09    //!< KMT09
	};

	StarFormationModel model;

	double nu_sf;
	double Po;
	double beta_press;
	double Accuracy_SFeqs;
	double gas_velocity_dispersion;
	double boost_starburst;
	double sigma_HI_crit;
	double clump_factor_KMT09;
	double sigma_crit_KMT09;

	bool angular_momentum_transfer;
};


class StarFormation {

public:
	StarFormation(StarFormationParameters parameters, RecyclingParameters recycleparams, std::shared_ptr<Cosmology> cosmology);

	/**
	 * All input quantities should be in comoving units.
	 */
	double star_formation_rate(double mcold, double mstars, double rgas, double rstars, double zgas, double z,
							   bool burst, double vgal, double &jrate, double jgas);

	double star_formation_rate_surface_density(double r, void * params);

	double fmol(double Sigma_gas, double Sigma_stars, double zgas, double r);

	double midplane_pressure(double Sigma_gas, double Sigma_stars, double r);

	double gd14_sigma_norm(double d_mw, double u_mw);

	unsigned long int get_integration_intervals() {
		return integrator.get_num_intervals();
	}

	void reset_integration_intervals() {
		return integrator.reset_num_intervals();
	}

	double molecular_hydrogen(double mcold, double mstars, double rgas, double rstars, double zgas, double z);

	double molecular_surface_density(double r, void * params);

	void get_molecular_gas(const GalaxyPtr &galaxy, double z, double *m_mol, double *m_atom, double *m_mol_b, double *m_atom_b);

	double ionised_gas_fraction(double mgas, double rgas, double z);

private:
	StarFormationParameters parameters;
	RecyclingParameters recycleparams;
	std::shared_ptr<Cosmology> cosmology;
	Integrator integrator;

};

}  // namespace shark

#endif /* INCLUDE_STAR_FORMATION_H_ */

