/*
 * star_formation.cpp
 *
 *  Created on: 17May,2017
 *      Author: clagos
 */

#include <cmath>

#include "components.h"
#include "logging.h"
#include "star_formation.h"
#include "numerical_constants.h"

namespace shark {

struct galaxy_properties_for_integration {
	double sigma_gas0;
	double sigma_star0;
	double re;
	double rse;
};

StarFormationParameters::StarFormationParameters(const Options &options) :
	Molecular_BR_Law(0),
	nu_sf(0),
	Po(0),
	beta_press(0),
	Accuracy_SFeqs(0),
	gas_velocity_dispersion(0)
{
	options.load("star_formation.Molecular_BR_law", Molecular_BR_Law);
	options.load("star_formation.nu_sf", nu_sf);
	options.load("star_formation.Po", Po);
	options.load("star_formation.beta_press", beta_press);
	options.load("star_formation.Accuracy_SFeqs", Accuracy_SFeqs);
	options.load("star_formation.gas_velocity_dispersion", gas_velocity_dispersion);
}


StarFormation::StarFormation(StarFormationParameters parameters, std::shared_ptr<Cosmology> cosmology) :
	parameters(parameters),
	cosmology(cosmology),
	integrator(1000)
{
	// no-op
}

double StarFormation::star_formation_rate(double mcold, double mstar, double rgas, double rstar, double z) {

	if (mcold <= 0 or rgas <= 0) {
		if(mcold > 0 && rgas <= 0){
			std::ostringstream os;
			os << "Galaxy mcold > 0 and rgas = 0";
			throw invalid_argument(os.str());
		}
		return 0;
	}

	/**
	 * All input quantities should be in physical units.
	 */
	// Define properties that are input for the SFR calculation.

	double re = cosmology->comoving_to_physical_size(rgas / constants::RDISK_HALF_SCALE, z);
	double rse = cosmology->comoving_to_physical_size(rstar / constants::RDISK_HALF_SCALE, z);

	double Sigma_gas = cosmology->comoving_to_physical_mass(mcold) / constants::PI2 / (re * re);
	double Sigma_star = 0;
	if(mstar > 0 and rstar > 0){
		Sigma_star = cosmology->comoving_to_physical_mass(mstar) / constants::PI2 / (rse * rse) ;
	}

	galaxy_properties_for_integration props = {
		Sigma_gas,
		Sigma_star,
		re,
		rse
	};

	struct StarFormationAndProps {
		StarFormation *star_formation;
		galaxy_properties_for_integration *props;
	};

	auto f = [](double r, void *ctx) -> double {
		StarFormationAndProps *sf_and_props = reinterpret_cast<StarFormationAndProps *>(ctx);
		return sf_and_props->star_formation->star_formation_rate_surface_density(r, sf_and_props->props);
	};

	double rmin = 0;
	double rmax = 5.0*re;

	StarFormationAndProps sf_and_props = {this, &props};
	// Adopt 5% accuracy for star formation solution.
	double result = integrator.integrate(f, &sf_and_props, rmin, rmax, 0.0, 0.05);

	/*int bins = 30;
	double integral = 0.0;
	double binr = (rmax-rmin)/bins;
	double rin = 0.0;

	for (int i=0; i<bins; i++){
		rin = rmin+(i+0.5)*binr;
		integral += star_formation_rate_surface_density(rin,sf_and_props.props) * binr;
	}
	double result = integral;*/

	// Avoid negative values.
	if(result < 0){
		result = 0.0;
	}

	if(mcold > 0 && result <= 0){
		std::ostringstream os;
		os << "Galaxy with SFR=0 and mcold " << mcold;
		throw invalid_argument(os.str());
	}

	double sft = mcold/result;

	return cosmology->physical_to_comoving_mass(result);

}

double StarFormation::star_formation_rate_surface_density(double r, void * params){

	using namespace constants;

	// apply molecular SF law
	auto props = reinterpret_cast<galaxy_properties_for_integration *>(params);

	double Sigma_gas = props->sigma_gas0 * std::exp(-r / props->re);

	// Avoid negative numbers.
	if(Sigma_gas < 0){
		Sigma_gas = 0;
	}

	double Sigma_stars = 0;

	// Define Sigma_stars only if stellar mass and radius are positive.
	if(props->rse > 0 and props->sigma_star0 > 0){
		Sigma_stars = props->sigma_star0 * std::exp(-r / props->rse);
	}

	double fracmol = fmol(Sigma_gas, Sigma_stars, r);
	double sfr_density = PI2 * parameters.nu_sf * fracmol * Sigma_gas * r; //Add the 2PI*r to Sigma_SFR to make integration.

	if(props->sigma_gas0 > 0 && sfr_density <= 0){
		std::ostringstream os;
		os << "Galaxy with SFR surface density =0 and cold gas surface density " << props->sigma_gas0;
		throw invalid_argument(os.str());
	}

	return sfr_density;
}

double StarFormation::molecular_surface_density(double r, void * params){

	using namespace constants;

	// apply molecular SF law
	auto props = reinterpret_cast<galaxy_properties_for_integration *>(params);

	double Sigma_gas = props->sigma_gas0 * std::exp(-r / props->re);

	// Avoid negative numbers
	if(Sigma_gas < 0){
		Sigma_gas = 0;
	}

	double Sigma_stars = 0;

	// Define Sigma_stars only if stellar mass and radius are positive.
	if(props->rse > 0 && props->sigma_star0 > 0){
		Sigma_stars = props->sigma_star0 * std::exp(-r / props->rse);
	}

	return PI2 * fmol(Sigma_gas, Sigma_stars, r) * Sigma_gas * r; //Add the 2PI*r to Sigma_SFR to make integration.
}

double StarFormation::fmol(double Sigma_gas, double Sigma_stars, double r){

	double rmol = std::pow((midplane_pressure(Sigma_gas,Sigma_stars,r)/parameters.Po),parameters.beta_press);

	double fmol = rmol/(1+rmol);

	// Avoid calculation errors.
	if(fmol > 1){
		return 1;
	}
	else if(fmol > 0 and fmol < 1){
		return fmol;
	}
	else{
		return 0;
	}
}

double StarFormation::midplane_pressure(double Sigma_gas, double Sigma_stars, double r){

	/**
	 * This function calculate the midplane pressure of the disk, and returns it in units of K/cm^-3.
	 */

	using namespace constants;

	double hstar = 0.14 * r; //scaleheight of the stellar disk; from Kregel et al. (2002).
	double veldisp_star = std::sqrt(PI * G * hstar * Sigma_stars); //stellar velocity dispersion in km/s.

	double star_comp = 0;

	if (Sigma_stars > 0 and veldisp_star > 0) {
		star_comp = (parameters.gas_velocity_dispersion / veldisp_star) * Sigma_stars;
	}

	double pressure = Pressure_Conv * Sigma_gas * (Sigma_gas + star_comp); //in units of K/cm^3.

	return pressure;
}

double StarFormation::molecular_hydrogen(double mcold, double mstar, double rgas, double rstar, double z) {

	if (mcold <= 0 or rgas <= 0) {
		return 0;
	}

	/**
	 * All input quantities should be in physical units.
	 */
	// Define properties that are input for the SFR calculation.

	double re = cosmology->comoving_to_physical_size(rgas / constants::RDISK_HALF_SCALE, z);
	double rse = cosmology->comoving_to_physical_size(rstar / constants::RDISK_HALF_SCALE, z);

	double Sigma_gas = cosmology->comoving_to_physical_mass(mcold) / constants::PI2 / (re * re);
	double Sigma_star = 0;
	if(mstar > 0 and rstar > 0){
		Sigma_star = cosmology->comoving_to_physical_mass(mstar) / constants::PI2 / (rse * rse) ;
	}

	galaxy_properties_for_integration props = {
		Sigma_gas,
		Sigma_star,
		re,
		rse
	};

	struct StarFormationAndProps {
		StarFormation *star_formation;
		galaxy_properties_for_integration *props;
	};

	auto f = [](double r, void *ctx) -> double {
		StarFormationAndProps *sf_and_props = reinterpret_cast<StarFormationAndProps *>(ctx);
		return sf_and_props->star_formation->molecular_surface_density(r, sf_and_props->props);
	};

	double rmin = 0;
	double rmax = 3.0*re;

	StarFormationAndProps sf_and_props = {this, &props};
	double result = integrator.integrate(f, &sf_and_props, rmin, rmax, 0.0, 0.05);

	// Avoid negative values.
	if(result <0){
		result = 0.0;
	}

	return cosmology->physical_to_comoving_mass(result);
}

void StarFormation::get_molecular_gas(const GalaxyPtr &galaxy, double z, double *m_mol, double *m_atom, double *m_mol_b, double *m_atom_b)
{
	*m_mol = 0;
	*m_atom = 0;
	*m_mol_b = 0;
	*m_atom_b = 0;

	if (galaxy->disk_gas.mass > 0) {
		*m_mol = molecular_hydrogen(galaxy->disk_gas.mass,galaxy->disk_stars.mass,galaxy->disk_gas.rscale, galaxy->disk_stars.rscale, z);
		*m_atom = galaxy->disk_gas.mass - *m_mol;
	}
	if (galaxy->bulge_gas.mass > 0) {
		*m_mol_b = molecular_hydrogen(galaxy->bulge_gas.mass,galaxy->bulge_stars.mass,galaxy->bulge_gas.rscale, galaxy->bulge_stars.rscale, z);
		*m_atom_b = galaxy->bulge_gas.mass - *m_mol_b;
	}
}

}  // namespace shark


