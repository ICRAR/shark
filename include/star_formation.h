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

#include "cosmology.h"
#include "integrator.h"
#include "options.h"

namespace shark {

class StarFormationParameters {

public:
	StarFormationParameters(const Options &options);

	int Molecular_BR_Law;
	double nu_sf;
	double Po;
	double beta_press;
	double Accuracy_SFeqs;
	double gas_velocity_dispersion;
};


class StarFormation {

public:
	StarFormation(StarFormationParameters parameters, std::shared_ptr<Cosmology> cosmology);

	/**
	 * All input quantities should be in comoving units.
	 */
	double star_formation_rate(double mcold, double mstars, double rgas, double rstars, double z);

	double star_formation_rate_surface_density(double r, void * params);

	double fmol(double Sigma_gas, double Sigma_stars, double r);

	double midplane_pressure(double Sigma_gas, double Sigma_stars, double r);

	unsigned long int get_integration_intervals() {
		return integrator.get_num_intervals();
	}

	void reset_integration_intervals() {
		return integrator.reset_num_intervals();
	}

	double molecular_hydrogen(double mcold, double mstars, double rgas, double rstars, double z);
	double molecular_surface_density(double r, void * params);

private:
	StarFormationParameters parameters;
	std::shared_ptr<Cosmology> cosmology;
	Integrator integrator;

};

}  // namespace shark

#endif /* INCLUDE_STAR_FORMATION_H_ */

