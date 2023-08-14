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
 */

#include <cmath>
#include <gsl/gsl_errno.h>

#include "galaxy.h"
#include "logging.h"
#include "numerical_constants.h"
#include "star_formation.h"
#include "utils.h"

namespace shark {

struct galaxy_properties_for_integration {
	double sigma_gas0;
	double sigma_star0;
	double re;
	double rse;
	double zgas;
	double torb;
	bool burst;
};

StarFormationParameters::StarFormationParameters(const Options &options)
{
	options.load("star_formation.model", model, true);
	options.load("star_formation.nu_sf", nu_sf, true);
	options.load("star_formation.angular_momentum_transfer", angular_momentum_transfer);

	options.load("star_formation.accuracy_sf_eqs", Accuracy_SFeqs);
	options.load("star_formation.boost_starburst", boost_starburst);
	options.load("star_formation.po", Po);
	options.load("star_formation.beta_press", beta_press);
	options.load("star_formation.gas_velocity_dispersion", gas_velocity_dispersion);
	options.load("star_formation.sigma_hi_crit", sigma_HI_crit);

	options.load("star_formation.clump_factor_kmt09", clump_factor_KMT09);

	options.load("star_formation.efficiency_sf_kd12", efficiency_sf);
	options.load("star_formation.gmc_surface_density", gmc_surface_density);

	// Convert surface density to internal code units.
	sigma_HI_crit = sigma_HI_crit * std::pow(constants::MEGA,2.0);
        gmc_surface_density = gmc_surface_density * std::pow(constants::MEGA,2.0);

	// Define critical density for the normal to starburst SF transition for the KMT09 model in Msun/Mpc^2.
	sigma_crit_KMT09 = 85.0 * std::pow(constants::MEGA , 2.0);
}


template <>
StarFormationParameters::StarFormationModel
Options::get<StarFormationParameters::StarFormationModel>(const std::string &name, const std::string &value) const {
	auto lvalue = lower(value);
	if (lvalue == "br06") {
		return StarFormationParameters::BR06;
	}
	else if (lvalue == "gd14") {
		return StarFormationParameters::GD14;
	}
	else if (lvalue == "k13") {
		return StarFormationParameters::K13;
	}
	else if (lvalue == "kmt09") {
		return StarFormationParameters::KMT09;
	}
	else if(lvalue == "kd12") {
		return StarFormationParameters::KD12;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value << ". Supported values are br06, gd11, k13, kmt09, kd12";
	throw invalid_option(os.str());
}

StarFormation::StarFormation(StarFormationParameters parameters, RecyclingParameters recycleparams, CosmologyPtr cosmology) :
	parameters(parameters),
	recycleparams(recycleparams),
	cosmology(std::move(cosmology)),
	integrator(1000)
{
	// no-op
}

double StarFormation::star_formation_rate(double mcold, double mstar, double rgas, double rstar, double zgas, double z,
					  bool burst, double vgal, double &jrate, double jgas) {

	if (std::isnan(rgas)) {
		throw invalid_argument("rgas is NaN, cannot calculate star formation rate");
	}

	if (mcold <= constants::EPS3 || rgas <= constants::tolerance) {
		if(mcold > constants::EPS3 && rgas <= 0){
			std::ostringstream os;
			os << "Galaxy mcold > 0 and rgas <0";
			throw invalid_argument(os.str());
		}
		if(mcold > constants::EPS3 && rgas <= constants::tolerance){
			std::ostringstream os;
			os << "Galaxy with extremely small size, rgas < 1e-10";
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
	if(mstar > 0 && rstar > 0){
		Sigma_star = cosmology->comoving_to_physical_mass(mstar) / constants::PI2 / (rse * rse) ;
	}

	//orbital time of the gas in Gyr
	double torb = 2.0 * constants::MPCKM2GYR * rgas / vgal;

	galaxy_properties_for_integration props = {
		Sigma_gas,
		Sigma_star,
		re,
		rse,
		zgas/recycleparams.zsun,
		torb,
		burst,
	};

	struct StarFormationAndProps {
		StarFormation *star_formation;
		galaxy_properties_for_integration *props;
	};

	auto f = [](double r, void *ctx) -> double {
		auto *sf_and_props = static_cast<StarFormationAndProps *>(ctx);
		return sf_and_props->star_formation->star_formation_rate_surface_density(r, sf_and_props->props);
	};

	double rmin = 0;
	double rmax = 5.0*re;

	StarFormationAndProps sf_and_props = {this, &props};

	double result = 0;
	try{
		result = integrator.integrate(f, &sf_and_props, rmin, rmax, 0.0, parameters.Accuracy_SFeqs);
	} catch (gsl_error &e) {
		auto gsl_errno = e.get_gsl_errno();
		std::ostringstream os;
		os << "SFR integration failed with GSL error number " << gsl_errno << ": ";
		os << gsl_strerror(gsl_errno) << ", reason=" << e.get_reason();
		os << ". We'll attempt manual integration now";
		LOG(warning) << os.str();

		// Perform manual integration.
		// TODO: check that error is affordable (i.e., maybe the error is really bad and the
		// program should stop)
		result = manual_integral(f, &sf_and_props, rmin, rmax);
	}

	// Avoid negative values.
	if(result < 0){
		result = 0.0;
	}

	result = cosmology->physical_to_comoving_mass(result);

	//Avoid AM calculation in the case of starbursts.
	if(!burst){
		// Check whether user wishes to calculate angular momentum transfer from gas to stars.
		if(parameters.angular_momentum_transfer){

			// TODO: it would be nice to somehow reuse some of the values from the previous integration
			// in here. At least initially during the first round the integration algorithm will run
			// over the same set of 'r' that it used during the first round of the previous integration,
			// so we could save ourselves lots of calculation by storing those values and reusing them here
			auto f_j = [](double r, void *ctx) -> double {
				auto *sf_and_props = static_cast<StarFormationAndProps *>(ctx);
				return r * sf_and_props->star_formation->star_formation_rate_surface_density(r, sf_and_props->props);
			};

			// React to integration errors by using a way-simpler 4-point manual integration
			double jSFR;
			try{
				jSFR = integrator.integrate(f_j, &sf_and_props, rmin, rmax, 0.0, parameters.Accuracy_SFeqs);
			} catch (gsl_error &e) {
				auto gsl_errno = e.get_gsl_errno();
				std::ostringstream os;
				os << "jSFR integration failed with GSL error number " << gsl_errno << ": ";
				os << gsl_strerror(gsl_errno) << ", reason=" << e.get_reason();
				os << ". We'll attempt manual integration now";
				LOG(warning) << os.str();

				// Perform manual integration.
				// TODO: check that error is affordable (i.e., maybe the error is really bad and the
				// program should stop)
				jSFR = manual_integral(f_j, &sf_and_props, rmin, rmax);
			}

			jrate = cosmology->physical_to_comoving_mass(jSFR) * vgal; //assumes a flat rotation curve.


			// Avoid negative values.
			if(jrate < 0){
				jrate = 0.0;
			}

			double effecj = jrate / result;

			//Assign maximum value to be jgas.
			if(effecj > jgas && jgas > 0){
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

	if(mcold > 0 && result <= 0 && mstar > 0){
		std::ostringstream os;
		os << "Galaxy with SFR=0, mcold " << mcold << ", mstar " << mstar << ", rgas" << rgas << ", rstar" << rstar;
		throw invalid_argument(os.str());
	}

	return result;

}

double StarFormation::star_formation_rate_surface_density(double r, void * params) const {

	using namespace constants;

	// apply molecular SF law
	auto props = static_cast<galaxy_properties_for_integration *>(params);

	double Sigma_gas = props->sigma_gas0 * std::exp(-r / props->re);

	// Avoid negative numbers.
	if(Sigma_gas < 0){
		Sigma_gas = 0;
	}

	double Sigma_stars = 0;

	// Define Sigma_stars only if stellar mass and radius are positive.
	if(props->rse > 0 && props->sigma_star0 > 0){
		Sigma_stars = props->sigma_star0 * std::exp(-r / props->rse);
	}

	double fracmol = fmol(Sigma_gas, Sigma_stars, props->zgas, r);

	double sfr_density = 0;

	if(parameters.model == StarFormationParameters::BR06 ||
			parameters.model == StarFormationParameters::GD14){
		sfr_density = PI2 * parameters.nu_sf * fracmol * Sigma_gas * r; //Add the 2PI*r to Sigma_SFR to make integration.
	}
	else if (parameters.model == StarFormationParameters::KD12){
		double tdep = kd12_taudep(Sigma_gas, props);
		sfr_density = PI2 * fracmol * Sigma_gas / tdep * r;
	}
	else if (parameters.model == StarFormationParameters::KMT09 || parameters.model == StarFormationParameters::K13){
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
	if(props->burst && parameters.model != StarFormationParameters::KD12){
		sfr_density = sfr_density * parameters.boost_starburst;
	}

	if((props->sigma_gas0 > 0 && fracmol > 0) && sfr_density <= 0){
		std::ostringstream os;
		os << "Galaxy with SFR surface density =0, cold gas surface density " << props->sigma_gas0 << " and fmol > 0";
		throw invalid_argument(os.str());
	}

	return sfr_density;
}

double StarFormation::molecular_surface_density(double r, void * params) const {

	using namespace constants;

	// apply molecular SF law
	auto props = static_cast<galaxy_properties_for_integration *>(params);

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

double StarFormation::fmol(double Sigma_gas, double Sigma_stars, double zgas, double r) const {

	double rmol = 0;

	if(parameters.model == StarFormationParameters::BR06 ||
			parameters.model == StarFormationParameters::KD12){
		rmol = std::pow((midplane_pressure(Sigma_gas,Sigma_stars,r)/parameters.Po),parameters.beta_press);
	}
	else if (parameters.model == StarFormationParameters::GD14){
		//Galaxy parameters
		double d_mw = zgas;
		double u_mw = Sigma_gas / constants::sigma_gas_mw;

		double alpha = 0.5 + 1/(1 + sqrt(u_mw * std::pow(d_mw,2.0)/600.0));
		rmol = std::pow(Sigma_gas / gd14_sigma_norm(d_mw, u_mw), alpha);
	}
	else if (parameters.model == StarFormationParameters::K13){

		double func = k13_fmol(zgas, Sigma_gas);
		rmol =  (1.0 - func) / func;
		if(rmol < 1e-4){
			rmol = 1e-4;
		}
	}
	else if (parameters.model == StarFormationParameters::KMT09){

		double func  = kmt09_fmol(zgas, Sigma_gas);
		rmol  = (1.0 - func) / func;
		if(rmol < 1e-4){
			rmol = 1e-4;
		}
	}

	double fmol = rmol/(1+rmol);

	// Avoid calculation errors.
	if(fmol > 1){
		return 1;
	}
	else if(fmol > 0 && fmol < 1){
		return fmol;
	}
	else{
		return 0;
	}
}

double StarFormation::midplane_pressure(double Sigma_gas, double Sigma_stars, double r) const {

	/**
	 * This function calculate the midplane pressure of the disk, and returns it in units of K/cm^-3.
	 */

	using namespace constants;

	double hstar = 0.14 * r; //scaleheight of the stellar disk; from Kregel et al. (2002).
	double veldisp_star = std::sqrt(PI * G * hstar * Sigma_stars); //stellar velocity dispersion in km/s.

	double star_comp = 0;

	if (Sigma_stars > 0 && veldisp_star > 0) {
		star_comp = (parameters.gas_velocity_dispersion / veldisp_star) * Sigma_stars;
	}

	double pressure = Pressure_Conv * Sigma_gas * (Sigma_gas + star_comp); //in units of K/cm^3.

	return pressure;
}

double StarFormation::gd14_sigma_norm(double d_mw, double u_mw) const {

	double g = sqrt(std::pow(d_mw,2.0) + 0.0289);

	double sigma_r1 = 50.0 / g * sqrt(0.01 + u_mw) / (1 + 0.69 * sqrt(0.01 + u_mw)) * std::pow(constants::MEGA,2.0); //In Msun/Mpc^2.

	return sigma_r1;
}

double StarFormation::kmt09_fmol(double zgas, double sigma_gas) const {

	double chi   = 0.77 * (1.0 + 3.1 * std::pow(zgas,0.365));
	double s     = std::log(1.0 + 0.6 * chi)/( 0.04 * parameters.clump_factor_KMT09 * sigma_gas/std::pow(constants::MEGA, 2.0) * zgas);
	double delta = 0.0712 * std::pow(0.1 / s + 0.675, -2.8);
	double func  = std::pow(1.0 + std::pow(0.75 * s / (1.0 + delta), -5.0), -0.2);

	return func;
}

double StarFormation::k13_fmol(double zgas, double sigma_gas) const {

	//Galaxy parameters
	double d_mw = zgas;
	double u_mw = sigma_gas / constants::sigma_gas_mw;
	double Sigma0 = sigma_gas /std::pow(constants::MEGA, 2.0); // gas surface density in Msun/pc^2

	// Calculate cold neutral medium densities in the regimes of hydrostatic and two-phase equilibrium.
	double ncnm_2p    = 23.0 * u_mw * std::pow((1.0 + 3.1 * std::pow(d_mw, 0.365))/4.1, -1) / 10.0; //in units of 10xcm^-3
	double ncnm_hydro = 0.0124068 * std::pow(Sigma0, 2.0) * (1.0 + std::pow(1.0 + 1250.56/std::pow(Sigma0, 2.0), 0.5)) / 10.0; //in units of 10xcm^-3

	// Assign the maximum of the two densities.
	double ncnm = std::max(ncnm_2p, ncnm_hydro);

	double Chi = 7.2 * u_mw/ncnm;

	double Tauc = 0.066 * parameters.clump_factor_KMT09 * d_mw * Sigma0;
	double sfac = std::log10(1 + 0.6 * Chi + 0.01 * std::pow(Chi,2.0)) / (0.6 * Tauc);

	double func = 1 - 0.75 * sfac / (1 + 0.25 * sfac);

	return func;

}

double StarFormation::kd12_taudep(double sigma_gas, void * params) const{

	using namespace constants;

	// apply molecular SF law
	auto props = static_cast<galaxy_properties_for_integration *>(params);

	// in Gyr
	double tdep = 3.5 * props->torb / parameters.efficiency_sf;

	double sigma_gmc = parameters.gmc_surface_density; 
	if(sigma_gmc < sigma_gas){
		sigma_gmc= sigma_gas;
	}

	// in Gyr
	double tgmc = 1.9 / parameters.efficiency_sf * (parameters.gas_velocity_dispersion / 10.0)
			* std::pow(sigma_gmc/1e14, -0.75) * std::pow(sigma_gas/1e13, -0.25);

	double t = tdep;
	if(t > tgmc){
		t = tgmc;
	}

	//Put an upper limit to something similar to the numbers observed in the outkirts of spirals
	if(t > 100){
		t=100;
	}

	return t;
}



double StarFormation::molecular_hydrogen(double mcold, double mstar, double rgas, double rstar, double zgas, double z,
		double &jmol, double jgas, double vgal, bool bulge, bool jcalc) {

	if (mcold <= 0 || rgas <= 0) {
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
	if(mstar > 0 && rstar > 0){
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
		auto *sf_and_props = static_cast<StarFormationAndProps *>(ctx);
		return sf_and_props->star_formation->molecular_surface_density(r, sf_and_props->props);
	};

	double rmin = 0;
	double rmax = 5.0*re;

	StarFormationAndProps sf_and_props = {this, &props};

	// React to integration errors by using a way-simpler 4-point manual integration
	double result = 0;
	try{
		result = integrator.integrate(f, &sf_and_props, rmin, rmax, 0.0, parameters.Accuracy_SFeqs);
	} catch (gsl_error &e) {
		auto gsl_errno = e.get_gsl_errno();
		std::ostringstream os;
		os << "jSFR integration failed with GSL error number " << gsl_errno << ": ";
		os << gsl_strerror(gsl_errno) << ", reason=" << e.get_reason();
		os << ". We'll attempt manual integration now";
		LOG(warning) << os.str();

		// Perform manual integration.
		// TODO: check that error is affordable (i.e., maybe the error is really bad and the
		// program should stop)
		result = manual_integral(f, &sf_and_props, rmin, rmax);
	}

	// Avoid negative values.
	if(result < 0){
		result = 0.0;
	}

	result = cosmology->physical_to_comoving_mass(result);

	//Avoid AM calculation in the case of starbursts.
	if(!bulge && jcalc && result > 0){
		// Check whether user wishes to calculate angular momentum transfer from gas to stars.
		if(parameters.angular_momentum_transfer){

			// TODO: it would be nice to somehow reuse some of the values from the previous integration
			// in here. At least initially during the first round the integration algorithm will run
			// over the same set of 'r' that it used during the first round of the previous integration,
			// so we could save ourselves lots of calculation by storing those values and reusing them here
			auto f_j = [](double r, void *ctx) -> double {
				auto *sf_and_props = static_cast<StarFormationAndProps *>(ctx);
				return r * sf_and_props->star_formation->molecular_surface_density(r, sf_and_props->props);
			};

			// React to integration errors by using a way-simpler 4-point manual integration
			jmol = 0;
			try{
				jmol = integrator.integrate(f_j, &sf_and_props, rmin, rmax, 0.0, parameters.Accuracy_SFeqs);
			} catch (gsl_error &e) {
				auto gsl_errno = e.get_gsl_errno();
				std::ostringstream os;
				os << "jSFR integration failed with GSL error number " << gsl_errno << ": ";
				os << gsl_strerror(gsl_errno) << ", reason=" << e.get_reason();
				os << ". We'll attempt manual integration now";
				LOG(warning) << os.str();

				// Perform manual integration.
				// TODO: check that error is affordable (i.e., maybe the error is really bad and the
				// program should stop)
				jmol = manual_integral(f_j, &sf_and_props, rmin, rmax);
			}

			jmol = cosmology->physical_to_comoving_mass(jmol) * vgal; //assumes a flat rotation curve.

			// Avoid negative values.
			if(jmol < 0){
				jmol = 0.0;
			}

			// Normalize by H2 mass.
			jmol = jmol / result;

			//Assign maximum value to be jgas.
			if(jmol > jgas && jgas > 0){
				jmol = jgas;
			}
		}
		else{
			// case that assumes stellar/gas components have the same sAM.
			jmol = jgas;
		}
	}
	else{
		//In the case of bulges or the case where we do not need to calculate j (jcalc == false) or the case where mmol == 0.
		jmol = 0;
	}

	return result;
}

double StarFormation::ionised_gas_fraction(double mgas, double rgas, double z) const
{
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

StarFormation::molecular_gas StarFormation::get_molecular_gas(const Galaxy &galaxy, double z, bool jcalc)
{
	double m_mol    = 0;
	double m_atom   = 0;
	double m_mol_b  = 0;
	double m_atom_b = 0;
	double j_mol     = 0;
	double j_atom    = 0;

	double f_ion = 0;
	double m_neutral = 0;
	double zgas = 0;
	double jgas = 0;
	double vgal = galaxy.vmax;

	// Apply ionised fraction correction only in the case of disks.
	if (galaxy.disk_gas.mass > 0) {
		jgas = galaxy.disk_gas.sAM;

		f_ion = ionised_gas_fraction(galaxy.disk_gas.mass, galaxy.disk_gas.rscale, z);
		zgas = galaxy.disk_gas.mass_metals / galaxy.disk_gas.mass;
		m_neutral = (1-f_ion) * galaxy.disk_gas.mass;

		m_mol = (1-f_ion) * molecular_hydrogen(galaxy.disk_gas.mass,galaxy.disk_stars.mass,galaxy.disk_gas.rscale, galaxy.disk_stars.rscale, zgas, z, j_mol, jgas, vgal, false, jcalc);
		m_atom = m_neutral - m_mol;

		if(jcalc){
			// Calculate specific AM of atomic gas.
			j_atom = (jgas * m_neutral - j_mol * m_mol) / m_atom;
		}
	}
	if (galaxy.bulge_gas.mass > 0) {
		double dummy_jmol;
		zgas = galaxy.bulge_gas.mass_metals / galaxy.bulge_gas.mass;
		m_mol_b = molecular_hydrogen(galaxy.bulge_gas.mass,galaxy.bulge_stars.mass,galaxy.bulge_gas.rscale, galaxy.bulge_stars.rscale, zgas, z, dummy_jmol, jgas, vgal, true, jcalc);
		m_atom_b = galaxy.bulge_gas.mass - m_mol_b;
	}

	return molecular_gas {m_mol, m_atom, m_mol_b, m_atom_b, j_mol, j_atom};
}

double StarFormation::manual_integral(func_t f, void * params, double rmin, double rmax) const
{
	double integral = 0;

	int nbins = 30;

	// Perform integral in bins of log(r+1).
	double rminl = std::log10(rmin+1);
	double rmaxl = std::log10(rmax+1);

	double rbin = (rmaxl - rminl) / nbins;

	for (int bin=0; bin<nbins-1; bin++){
		double ri = rminl + rbin * bin;
		double rf = rminl + rbin * (bin+1);

		ri = std::pow(10.0,ri) - 1.0;
		rf = std::pow(10.0,rf) - 1.0;
		double rx =(rf + ri) * 0.5 ;

		integral += f(rx, &params) * (rf - ri);
	}

	// Avoid negative numbers.
	if(integral< 0){
		integral = 0;
	}

	return integral;
}

}  // namespace shark


