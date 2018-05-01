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
	double zgas;
	bool burst;
};

StarFormationParameters::StarFormationParameters(const Options &options) :
	model(BR06),
	nu_sf(0),
	Po(0),
	beta_press(0),
	Accuracy_SFeqs(0.05),
	gas_velocity_dispersion(0),
	sigma_HI_crit(0),
	boost_starburst(1),
	clump_factor_KMT09(1),
	sigma_crit_KMT09(0),
	angular_momentum_transfer(false)
{
	options.load("star_formation.model", model, true);
	options.load("star_formation.nu_sf", nu_sf, true);
	options.load("star_formation.angular_momentum_transfer", angular_momentum_transfer);

	options.load("star_formation.Accuracy_SFeqs", Accuracy_SFeqs);
	options.load("star_formation.boost_starburst", boost_starburst);
	options.load("star_formation.Po", Po);
	options.load("star_formation.beta_press", beta_press);
	options.load("star_formation.gas_velocity_dispersion", gas_velocity_dispersion);
	options.load("star_formation.sigma_HI_crit", sigma_HI_crit);

	options.load("star_formation.clump_factor_KMT09", clump_factor_KMT09);

	// Convert surface density to internal code units.
	sigma_HI_crit = sigma_HI_crit * std::pow(constants::MEGA,2.0);

	// Define critical density for the normal to starburst SF transition for the KMT09 model in Msun/Mpc^2.
	sigma_crit_KMT09 = 85.0 * std::pow(constants::MEGA , 2.0);
}


template <>
StarFormationParameters::StarFormationModel
Options::get<StarFormationParameters::StarFormationModel>(const std::string &name, const std::string &value) const {
	if ( value == "BR06" ) {
		return StarFormationParameters::BR06;
	}
	else if ( value == "GD14" ) {
		return StarFormationParameters::GD14;
	}
	else if (value == "K13"){
		return StarFormationParameters::K13;
	}
	else if (value == "KMT09"){
		return StarFormationParameters::KMT09;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value << ". Supported values are BR06, GK11, K13 or KMT09";
	throw invalid_option(os.str());
}

StarFormation::StarFormation(StarFormationParameters parameters, RecyclingParameters recycleparams, std::shared_ptr<Cosmology> cosmology) :
	parameters(parameters),
	recycleparams(recycleparams),
	cosmology(cosmology),
	integrator(1000)
{
	// no-op
}

double StarFormation::star_formation_rate(double mcold, double mstar, double rgas, double rstar, double zgas, double z,
								          bool burst, double vgal, double &jrate, double jgas) {

	if (mcold <= constants::EPS3 or rgas <= constants::EPS6 ) {
		if(mcold > constants::EPS3 and rgas <= 0){
			std::ostringstream os;
			os << "Galaxy mcold > 0 and rgas <0";
			throw invalid_argument(os.str());
		}
		if(mcold > constants::EPS3 and rgas <= constants::EPS6){
			std::ostringstream os;
			os << "Galaxy with extremely small size, rgas < 1e-6";
			//throw invalid_argument(os.str());
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
		zgas/recycleparams.zsun,
		burst,
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
	double result = integrator.integrate(f, &sf_and_props, rmin, rmax, 0.0, parameters.Accuracy_SFeqs);

	// Avoid negative values.
	if(result < 0){
		result = 0.0;
	}

	result = cosmology->physical_to_comoving_mass(result);

	//Avoid AM calculation in the case of starbursts.
	if(!burst){
		// Check whether user wishes to calculate angular momentum transfer from gas to stars.
		if(parameters.angular_momentum_transfer){
			auto f_j = [](double r, void *ctx) -> double {
				StarFormationAndProps *sf_and_props = reinterpret_cast<StarFormationAndProps *>(ctx);
				return r * sf_and_props->star_formation->star_formation_rate_surface_density(r, sf_and_props->props);
			};

			StarFormationAndProps sf_and_props = {this, &props};
			double jSFR = integrator.integrate(f_j, &sf_and_props, rmin, rmax, 0.0, parameters.Accuracy_SFeqs);
			jrate = cosmology->physical_to_comoving_mass(jSFR) * vgal; //assumes a flat rotation curve.


			// Avoid negative values.
			if(jrate < 0){
				jrate = 0.0;
			}

			double effecj = jrate / result;

			//Assign maximum value to be jgas.
			if(effecj > jgas and jgas > 0){
				jrate = result * jgas;
			}
		}
		else{
			// case that assumes stellar/gas components have the same sAM.
			jrate = result * jgas;
		}
	}
	else{
		//In the case of starbursts.
		jrate = 0;
	}

	if(mcold > 0 && result <= 0){
		std::ostringstream os;
		os << "Galaxy with SFR=0 and mcold " << mcold;
		throw invalid_argument(os.str());
	}

	return result;

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

	double fracmol = fmol(Sigma_gas, Sigma_stars, props->zgas, r);

	double sfr_density = 0;

	if(parameters.model == StarFormationParameters::BR06 or parameters.model == StarFormationParameters::GD14){
		sfr_density = PI2 * parameters.nu_sf * fracmol * Sigma_gas * r; //Add the 2PI*r to Sigma_SFR to make integration.
	}
	else if (parameters.model == StarFormationParameters::KMT09 or parameters.model == StarFormationParameters::K13){
		double sfr_ff = 0;

		if(Sigma_gas < parameters.sigma_crit_KMT09){
			sfr_ff = std::pow(Sigma_gas/parameters.sigma_crit_KMT09, -0.33);
		}
		else{
			sfr_ff = std::pow(Sigma_gas/parameters.sigma_crit_KMT09, 0.33);
		}

		sfr_density = PI2 * fracmol * sfr_ff * Sigma_gas / 2.6 * r;
	}

	// If the star formation mode is starburst, then apply boosting in star formation.
	if(props->burst){
		sfr_density = sfr_density * parameters.boost_starburst;
	}

	if((props->sigma_gas0 > 0 and fracmol > 0) and sfr_density <= 0){
		std::ostringstream os;
		os << "Galaxy with SFR surface density =0, cold gas surface density " << props->sigma_gas0 << " and fmol > 0";
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

	return PI2 * fmol(Sigma_gas, Sigma_stars, props->zgas, r) * Sigma_gas * r; //Add the 2PI*r to Sigma_SFR to make integration.
}

double StarFormation::fmol(double Sigma_gas, double Sigma_stars, double zgas, double r){

	double rmol = 0;

	if(parameters.model == StarFormationParameters::BR06){
		rmol = std::pow((midplane_pressure(Sigma_gas,Sigma_stars,r)/parameters.Po),parameters.beta_press);
	}
	else if (parameters.model == StarFormationParameters::GD14){
		double d_mw = zgas;
		double u_mw = Sigma_gas / constants::sigma_gas_mw;

		rmol = Sigma_gas / gd14_sigma_norm(d_mw, u_mw);
	}
	else if (parameters.model == StarFormationParameters::K13){
		//TODO
	}
	else if (parameters.model == StarFormationParameters::KMT09){

		double chi   = 0.77 * (1.0 + 3.1 * std::pow(zgas,0.365));
		double s     = std::log(1.0 + 0.6 * chi)/( 0.04 * parameters.clump_factor_KMT09 * Sigma_gas/std::pow(constants::MEGA, 2.0) * zgas);
		double delta = 0.0712 * std::pow(0.1 / s + 0.675, -2.8);
		double func  = std::pow(1.0 + std::pow(0.75 * s / (1.0 + delta), -5.0), -0.2);
		rmol  = (1.0 - func) / func;
		if(rmol < constants::EPS4){
			rmol = constants::EPS4;
		}
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

double StarFormation::gd14_sigma_norm(double d_mw, double u_mw){

	double alpha = 0.5 + 1/(1 + sqrt(u_mw * std::pow(d_mw,2.0)/600.0));

	double g = sqrt(std::pow(d_mw,2.0) + 0.0289);

	double sigma_r1 = 50.0 / g * sqrt(0.01 + u_mw) / (1 + 0.69 * sqrt(0.01 + u_mw)) * std::pow(constants::MEGA,2.0); //In Msun/Mpc^2.

	return sigma_r1;
}

double StarFormation::molecular_hydrogen(double mcold, double mstar, double rgas, double rstar, double zgas, double z) {

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
		rse,
		zgas/recycleparams.zsun
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
	double zgas = 0;

	// Apply ionised fraction correction only in the case of disks.
	if (galaxy->disk_gas.mass > 0) {
		f_ion = ionised_gas_fraction(galaxy->disk_gas.mass, galaxy->disk_gas.rscale, z);
		zgas = galaxy->disk_gas.mass_metals / galaxy->disk_gas.mass;
		m_neutral = (1-f_ion) * galaxy->disk_gas.mass;
		*m_mol = (1-f_ion) * molecular_hydrogen(galaxy->disk_gas.mass,galaxy->disk_stars.mass,galaxy->disk_gas.rscale, galaxy->disk_stars.rscale, zgas, z);
		*m_atom = m_neutral - *m_mol;
	}
	if (galaxy->bulge_gas.mass > 0) {
		zgas = galaxy->bulge_gas.mass_metals / galaxy->bulge_gas.mass;
		*m_mol_b = molecular_hydrogen(galaxy->bulge_gas.mass,galaxy->bulge_stars.mass,galaxy->bulge_gas.rscale, galaxy->bulge_stars.rscale, zgas, z);
		*m_atom_b = galaxy->bulge_gas.mass - *m_mol_b;
	}
}

}  // namespace shark


