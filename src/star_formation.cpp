/*
 * star_formation.cpp
 *
 *  Created on: 17May,2017
 *      Author: clagos
 */

#include <cmath>
#include <gsl/gsl_integration.h>

#include "star_formation.h"
#include "numerical_constants.h"

namespace shark {

struct galaxy_properties_for_integration {
	double mcold;
	double mstar;
	double rgas;
	double rstar;
};

StarFormationParameters::StarFormationParameters(const std::string &filename) :
	Options(filename),
	Molecular_BR_Law(0),
	nu_sf(0),
	Po(0),
	beta_press(0),
	Accuracy_SFeqs(0),
    gas_velocity_dispersion(0)
{
	load("star_formation.Molecular_BR_law", Molecular_BR_Law);
	load("star_formation.nu_sf", nu_sf);
	load("star_formation.Po", Po);
	load("star_formation.beta_press", beta_press);
	load("star_formation.Accuracy_SFeqs", Accuracy_SFeqs);
	load("star_formation.gas_velocity_dispersion", gas_velocity_dispersion);
}


StarFormation::StarFormation(StarFormationParameters parameters, std::shared_ptr<Cosmology> cosmology) :
	parameters(parameters),
	cosmology(cosmology)
{
	// no-op
}

double StarFormation::star_formation_rate(double mcold, double mstar, double rgas, double rstar, double z) {

	int smax = 1000;
	gsl_integration_workspace * w
	    = gsl_integration_workspace_alloc (smax);

	gsl_function F;

	/**
	 * All input quantities should be in physical units.
	 */
	galaxy_properties_for_integration props = {
		cosmology->comoving_to_physical_mass(mcold),
		cosmology->comoving_to_physical_mass(mstar),
		cosmology->comoving_to_physical_size(rgas, z),
		cosmology->comoving_to_physical_size(rstar, z)
	};

	struct StarFormationAndProps {
		StarFormation *star_formation;
		galaxy_properties_for_integration *props;
	};

	StarFormationAndProps sf_and_props = {this, &props};
	double (*f)(double, void*) = [](double r, void *ctx) -> double {
		StarFormationAndProps *sf_and_props = reinterpret_cast<StarFormationAndProps *>(ctx);
		return sf_and_props->star_formation->star_formation_rate_surface_density(r, sf_and_props->props);
	};
	F.function = f;
	F.params = &sf_and_props;

	double result, error;
	double rmin = 0;
	double rmax = 10*rgas;

	/**
	 * Here, we integrate the SFR surface density profile out to rmax.
	 */
	gsl_integration_qags (&F, rmin, rmax, 0, 1e-7, smax,
	                        w, &result, &error);

	gsl_integration_workspace_free (w);

	return cosmology->physical_to_comoving_mass(result);

}

double StarFormation::star_formation_rate_surface_density(double r, void * params){

	// apply molecular SF law
	auto props = reinterpret_cast<galaxy_properties_for_integration *>(params);

	double re = props->rgas/constants::RDISK_HALF_SCALE;
	double Sigma_gas = props->mcold/constants::PI2/std::pow(re,2)*std::exp(-r/re);

	double rse = props->rstar/constants::RDISK_HALF_SCALE;
	double Sigma_stars = props->mstar/constants::PI2/std::pow(rse,2)*std::exp(-r/rse);

	return parameters.nu_sf*fmol(Sigma_gas,Sigma_stars,r)*Sigma_gas;
}

double StarFormation::fmol(double Sigma_gas, double Sigma_stars, double r){

	double rmol = std::pow((midplane_pressure(Sigma_gas,Sigma_stars,r)/parameters.Po),parameters.beta_press);

	return rmol/(1+rmol);
}

double StarFormation::midplane_pressure(double Sigma_gas, double Sigma_stars, double r){

	/**
	 * This function calculate the midplane pressure of the disk, and returns it in units of K/cm^-3.
	 */
	//I need to first convert all my quantities to physical quantities. Do this by calling functions from cosmology.!!!

	using namespace constants;

	double hstar = 0.14*r; //scaleheight of the stellar disk; from Kregel et al. (2002).

	double veldisp_star = std::sqrt(PI*G*hstar*Sigma_stars); //stellar velocity dispersion in km/s.

	double pressure = PIO2*G_MPCGYR2*Sigma_gas*(Sigma_gas+(parameters.gas_velocity_dispersion/veldisp_star)*Sigma_stars); //in units of Msun/Mpc/Gyr^2.

	return pressure*Pressure_SimUnits_cgs/k_Boltzmann_erg; //pressure in units of K/cm^3.
}
}  // namespace shark


