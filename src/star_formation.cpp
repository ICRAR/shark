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
	bool burst;
};

StarFormationParameters::StarFormationParameters(const Options &options) :
	model(BR06),
	nu_sf(0),
	Po(0),
	beta_press(0),
	Accuracy_SFeqs(0),
	gas_velocity_dispersion(0),
	sigma_HI_crit(0),
	boost_starburst(1)
{
	options.load("star_formation.model", model);
	options.load("star_formation.nu_sf", nu_sf);
	options.load("star_formation.Po", Po);
	options.load("star_formation.beta_press", beta_press);
	options.load("star_formation.Accuracy_SFeqs", Accuracy_SFeqs);
	options.load("star_formation.gas_velocity_dispersion", gas_velocity_dispersion);
	options.load("star_formation.boost_starburst", boost_starburst);
	options.load("star_formation.sigma_HI_crit", sigma_HI_crit);

	// Convert surface density to internal code units.
	sigma_HI_crit = sigma_HI_crit * std::pow(constants::MEGA,2.0);
}


template <>
StarFormationParameters::StarFormationModel
Options::get<StarFormationParameters::StarFormationModel>(const std::string &name, const std::string &value) const {
	if ( value == "BR06" ) {
		return StarFormationParameters::BR06;
	}
	else if ( value == "GK11" ) {
		return StarFormationParameters::GK11;
	}
	else if (value == "K13"){
		return StarFormationParameters::K13;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value << ". Supported values are BR06, GK11 and K13";
	throw invalid_option(os.str());
}

StarFormation::StarFormation(StarFormationParameters parameters, std::shared_ptr<Cosmology> cosmology) :
	parameters(parameters),
	cosmology(cosmology),
	integrator(1000)
{
	// no-op
}

double StarFormation::star_formation_rate(double mcold, double mstar, double rgas, double rstar, double z, bool burst) {

	if (mcold <= constants::EPS3 or rgas <= 0) {
		if(mcold > constants::EPS3 && rgas <= 0){
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
		rse,
		burst
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

	// If the star formation mode is starburst, then apply boosting in star formation.
	if(props->burst){
		sfr_density = sfr_density * parameters.boost_starburst;
	}

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

	// Avoid negative numbers.
	if(Sigma_gas < 0){
		Sigma_gas = 0;
	}

	// Check for low surface densities..
	if(Sigma_gas < parameters.sigma_HI_crit){
		return 0;
	}

	double Sigma_stars = 0;

	// Define Sigma_stars only if stellar mass and radius are positive.
	if(props->rse > 0 && props->sigma_star0 > 0){
		Sigma_stars = props->sigma_star0 * std::exp(-r / props->rse);
	}

	return PI2 * fmol(Sigma_gas, Sigma_stars, r) * Sigma_gas * r; //Add the 2PI*r to Sigma_SFR to make integration.
}

double StarFormation::fmol(double Sigma_gas, double Sigma_stars, double r){

	double rmol = 0;

	if(parameters.model == StarFormationParameters::BR06){
		rmol = std::pow((midplane_pressure(Sigma_gas,Sigma_stars,r)/parameters.Po),parameters.beta_press);
	}
	else if (parameters.model == StarFormationParameters::GK11){
		//TODO
	}
	else if (parameters.model == StarFormationParameters::K13){
		//TODO
	}

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

double StarFormation::ionised_gas_fraction(double mgas, double rgas, double z){

	double re = cosmology->comoving_to_physical_size(rgas / constants::RDISK_HALF_SCALE, z);

	double sigma0 = mgas/constants::PI2 / (re * re);

	double r_thresh = -re * std::log(parameters.sigma_HI_crit / sigma0);

	double m_in = mgas * ( 1- (1 + r_thresh/re) * std::exp(-r_thresh / re));

	double f_ion = (mgas - m_in) / mgas;

	if(f_ion < 0){
		f_ion = 0;
	}
	else if (f_ion >1){
		std::ostringstream os;
		os << "Galaxy with ionised gas fraction >1! Not possible.";
		throw invalid_argument(os.str());
	}

	return f_ion;

}

void StarFormation::get_molecular_gas(const GalaxyPtr &galaxy, double z, double *m_mol, double *m_atom, double *m_mol_b, double *m_atom_b)
{
	*m_mol = 0;
	*m_atom = 0;
	*m_mol_b = 0;
	*m_atom_b = 0;

	double f_ion = 0;
	double m_neutral = 0;

	// Apply ionised fraction correction only in the case of disks.
	if (galaxy->disk_gas.mass > 0) {
		f_ion = ionised_gas_fraction(galaxy->disk_gas.mass, galaxy->disk_gas.rscale, z);
		m_neutral = (1-f_ion) * galaxy->disk_gas.mass;
		*m_mol = (1-f_ion) * molecular_hydrogen(galaxy->disk_gas.mass,galaxy->disk_stars.mass,galaxy->disk_gas.rscale, galaxy->disk_stars.rscale, z);
		*m_atom = m_neutral - *m_mol;
	}
	if (galaxy->bulge_gas.mass > 0) {
		*m_mol_b = molecular_hydrogen(galaxy->bulge_gas.mass,galaxy->bulge_stars.mass,galaxy->bulge_gas.rscale, galaxy->bulge_stars.rscale, z);
		*m_atom_b = galaxy->bulge_gas.mass - *m_mol_b;
	}
}

}  // namespace shark


