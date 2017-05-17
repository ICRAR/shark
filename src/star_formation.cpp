/*
 * star_formation.cpp
 *
 *  Created on: 17May,2017
 *      Author: clagos
 */

#include <cmath>

#include "cosmology.h"
#include "star_formation.h"
#include "numerical_constants.h"

namespace shark {

StarFormationParameters::StarFormationParameters(const std::string &filename) :
	Options(filename),
	Molecular_BR_Law(0),
	nu_sf(0),
	Po(0),
	beta_press(0),
	Accuracy_SFeqs(0)

{
	load("star_formation.Molecular_BR_law", Molecular_BR_Law);
	load("star_formation.nu_sf", nu_sf);
	load("star_formation.Po", Po);
	load("star_formation.beta_press", beta_press);
	load("star_formation.Accuracy_SFeqs", Accuracy_SFeqs);

}


StarFormation::StarFormation(StarFormationParameters parameters) :
	parameters(parameters)
{
	// no-op
}

double StarFormation::star_formation_rate(double mcold, double mstar, double rgas, double rstar) {
	return 0;
}

double star_formation_rate_surface_density(double mcold, double mstar, double rgas, double rstar){
//	double reff_gas = rgas * NumericalConstants::MEGA/CosmologicalParameters::Hubble_h;
	return 0;
}
}  // namespace shark


