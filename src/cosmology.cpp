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

#include <cerrno>
#include <cmath>
#include <fstream>
#include <sstream>
#include <map>
#include <tuple>
#include <stdexcept>

#include "cosmology.h"
#include "data.h"
#include <gsl/gsl_integration.h>
#include "logging.h"
#include "numerical_constants.h"
#include "utils.h"

namespace shark {

enum power_spectrum_type
{
	MILLENIUM = 0,
	MILLENIUM_GAS,
	PLANCK14,
	PLANCK15,
	PLANCK15extended,
	WDM_DOVE
};


template <>
power_spectrum_type Options::get<power_spectrum_type>(const std::string &name, const std::string &value) const
{
	auto lvalue = lower(value);
	if (lvalue == "millenium") {
		return power_spectrum_type::MILLENIUM;
	}
	else if (lvalue == "millenium_gas") {
		return power_spectrum_type::MILLENIUM_GAS;
	}
	else if (lvalue == "planck14") {
		return power_spectrum_type::PLANCK14;
	}
	else if (lvalue == "planck15") {
		return power_spectrum_type::PLANCK15;
	}
	else if (lvalue == "planck15extended") {
		return power_spectrum_type::PLANCK15extended;
	}
	else if (lvalue == "wdm_dove") {
		return power_spectrum_type::WDM_DOVE;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value << ". ";
	os << "Supported values are millenium, millenium_gas, planck14, planck15, planck15extended and wdm_dove";
	throw invalid_option(os.str());
}

static
std::map<power_spectrum_type, std::string> power_spectrum_files {
	{MILLENIUM, "pk_Mill.dat"},
	{MILLENIUM_GAS, "pk_MillGas_norm.dat"},
	{PLANCK14, "pk_Planck_norm.dat"},
	{PLANCK15, "pk_Planck_SURFS_norm.dat"},
	{PLANCK15extended, "pk_planck_surfs_extended.dat"},
	{WDM_DOVE, "pk_WDMDove.dat"}
};

CosmologicalParameters::CosmologicalParameters(const Options &options)
{
	power_spectrum_type ps_type = PLANCK15;
	options.load("cosmology.omega_m", OmegaM);
	options.load("cosmology.omega_b", OmegaB);
	options.load("cosmology.omega_l", OmegaL);
	options.load("cosmology.n_s", n_s);
	options.load("cosmology.sigma8", sigma8);
	options.load("cosmology.hubble_h", Hubble_h);
	options.load("cosmology.power_spectrum", ps_type);
	load_tables(get_static_data_filepath(std::string("Power_Spec/") + power_spectrum_files[ps_type]));
}

void CosmologicalParameters::load_tables(const std::string &power_spec_file)
{

	LOG(debug) << "Reading table " << power_spec_file ;

	std::ifstream f = open_file(power_spec_file);
	std::string line;
	while ( std::getline(f, line) ) {

		trim(line);
		if (empty_or_comment(line)) {
			continue;
		}

		double k,p;
		std::istringstream iss(line);
		iss >> k >> p;

		power_spectrum.k.push_back(k);
		power_spectrum.p.push_back(p);

	}
	f.close();

}

Cosmology::Cosmology(CosmologicalParameters parameters) :
	parameters(std::move(parameters))
{
    initialise_lookup_tables();
}

void Cosmology::initialise_lookup_tables()
{
    // Precompute sigma(M) over a logarithmic mass grid. 
	// This is done to speed up the code later on.
    // Masses are in Msun/h.

    constexpr double log_m_min = 4.0;
    constexpr double log_m_max = 17.0;
    constexpr double dlogm = 0.01;

    const int n = static_cast<int>(
        std::round((log_m_max - log_m_min) / dlogm)) + 1;

    log_mass_grid.reserve(n);
    log_sigma_grid.reserve(n);

    for (int i = 0; i < n; ++i) {

        const double logm = log_m_min + i * dlogm;
        const double mass = std::pow(10.0, logm);

        const double sigma = sigma_mass_direct(mass);

        if (!std::isfinite(sigma) || sigma <= 0.0) {
            throw std::runtime_error(
                "Invalid sigma(M) while initialising lookup table");
        }

        log_mass_grid.push_back(logm);
        log_sigma_grid.push_back(std::log10(sigma));
    }


	// ------------------------------------------------------------
	// Growth-factor lookup table
	// ------------------------------------------------------------

	// Tabulate D(a) on a uniform logarithmic scale-factor grid.
	// Covers a = 0.001 to 1, corresponding to z = 999 to 0.

	constexpr double log_a_min = -3.0;
	constexpr double log_a_max = 0.0;
	constexpr int n_growth = 1001;

	const double dloga =
    	(log_a_max - log_a_min) / (n_growth - 1);

	log_a_grid.reserve(n_growth);
	log_growth_grid.reserve(n_growth);

	for (int i = 0; i < n_growth; ++i) {

    	const double loga = log_a_min + i * dloga;
	    const double a = std::pow(10.0, loga);
    	const double z = 1.0 / a - 1.0;

    	const double D = growth_factor_direct(z);

    	if (!std::isfinite(D) || D <= 0.0) {
        	throw std::runtime_error(
            	"Invalid growth factor while initialising lookup table");
	    }

    	log_a_grid.push_back(loga);
    	log_growth_grid.push_back(std::log10(D));
	}	

}

double Cosmology::comoving_to_physical_angularmomentum(double r, double z) const
{
	return r / std::pow(parameters.Hubble_h, 2) / (1 + z);
}

double Cosmology::comoving_to_physical_size(double r, double z) const {
	// We DO NOT need to divide by (1+z) because the sizes of galaxies and halos are calculated
	// from their mass and velocity, and the latter is in physical units from the VELOCIraptor catalogue.
	return r / parameters.Hubble_h; ////(1+z);
}

double Cosmology::physical_to_comoving_size(double r, double z) const {
	// We DO NOT need to multiply by (1+z) because the sizes of galaxies and halos are calculated
	// from their mass and velocity, and the latter is in physical units from the VELOCIraptor catalogue.
	return r * parameters.Hubble_h; ////(1+z);
}

double Cosmology::comoving_to_physical_velocity(double v, double z) const {
	return v;
}

double Cosmology::comoving_to_physical_mass(double m) const {
	return m/parameters.Hubble_h;
}

double Cosmology::physical_to_comoving_mass(double m) const {
	return m*parameters.Hubble_h;
}

double Cosmology::convert_redshift_to_age(double z) const {

	/**
	 * Function that calculates an age of the universe from a redshift.
	 */

	using namespace constants;

	double Hubble_Time=1.0/H0100PGYR; //The Hubble time for H_0=100km/s/Mpc.

	double err = 10e-5;
	double t;

	double a = 1/(1+z);

	if(std::abs(1-parameters.OmegaM)< err && parameters.OmegaL == 0 ){//Einstein-de Sitter universe.
		t = Hubble_Time*2*a*std::sqrt(a)/(3*parameters.Hubble_h);
	}
	else if (parameters.OmegaM < 1 && parameters.OmegaL == 0){// Open universe with no cosmological constant.
		double zplus1 = 1/a;
		t = Hubble_Time*(parameters.OmegaM/(2.0*parameters.Hubble_h* std::pow((1-parameters.OmegaM),1.5)))*(2.0*std::sqrt(1.0-parameters.OmegaM)*
				std::sqrt(parameters.OmegaM*(zplus1-1.0)+1.0)/(parameters.OmegaM
	            *zplus1)-std::acosh((parameters.OmegaM*(zplus1-1.0)-parameters.OmegaM+2.0)/(parameters.OmegaM*zplus1)));
	}
	else if(std::abs(1 - (parameters.OmegaM + parameters.OmegaL)) < err){//Flat with non-zero lambda.
		t = Hubble_Time*(2/(3*parameters.Hubble_h*std::sqrt(1-parameters.OmegaM)))*std::asinh(std::sqrt((1.0/parameters.OmegaM-1.0)*a)*a);
	}
	else{
		std::ostringstream os;
		os << "Error in age of the universe calculation -- not coded for this cosmology" << std::strerror(errno);
		throw std::runtime_error(os.str());
	}

	return t;
}


double Cosmology::convert_age_to_redshift_lcdm(double t) const {

	/**
	 * Function that calculates the redshift given an age. However, this function only applies to a flat universe with non-zero lambda.
	 */

	using namespace constants;

	double Hubble_Time=1.0/H0100PGYR; //The Hubble time for H_0=100km/s/Mpc.

	double a = std::pow(std::sinh(t / (Hubble_Time*(2/(3*parameters.Hubble_h*std::sqrt(1-parameters.OmegaM))))) / std::sqrt((1.0/parameters.OmegaM-1.0)) , 2.0/3.0);

	double z = 1.0 / a - 1;

	return z;
}


double Cosmology::expansion_factor(double t) const {

	using namespace constants;

	double err = std::pow(10,-5);
	double Hubble_Time=1.0/H0100PGYR; //The Hubble time for H_0=100km/s/Mpc.

	double a;

	if(std::abs(1-parameters.OmegaM)< err && parameters.OmegaL == 0){//Einstein-de Sitter universe.
		a = std::pow((1.5*t*parameters.Hubble_h/Hubble_Time),2/3);
	}
	else if(std::abs(1-parameters.OmegaM+parameters.OmegaL) < err){//Flat with non-zero lambda.
		double y = 1.5*t*parameters.Hubble_h*std::sqrt(1.0-parameters.OmegaM)/Hubble_Time;
		double sinhy = 0.5*(std::exp(y)-std::exp(-y));
		a = std::pow((sinhy*(1.0-parameters.OmegaM)*std::pow(parameters.OmegaM,2)),(2.0/3.0))/(parameters.OmegaM*(1.0-parameters.OmegaM));
	}
	else{
		std::ostringstream os;
		os << "Error in expansion factor calculation -- not coded for this cosmology" << std::strerror(errno);
		throw std::runtime_error(os.str());
	}

	return a;
}

double Cosmology::hubble_parameter (double z) const {
	double H2 = (parameters.OmegaM * std::pow(1.0 + z, 3.0) + parameters.OmegaL);
	return parameters.Hubble_h * 100.0 * std::sqrt(H2);
}

double Cosmology::critical_density (double z) const {

	// Function returns the critical density in units of Msun/cMpc^3
	
	auto h = hubble_parameter(z) / 100.0; // we want h not H.

	return 2.7754e11 * std::pow(h, 2);
}

double Cosmology::growth_factor_direct(double z) const
{
	/*
	* Linear growth factor for a flat Lambda-CDM cosmology.
	* Normalised such that D(z=0) = 1.
	*/

	const double Om = parameters.OmegaM;
	const double Ol = parameters.OmegaL;

	if (z < 0.0 || Om <= 0.0 || Ol < 0.0 ||
		std::abs(Om + Ol - 1.0) > 1.e-4) {
			throw std::runtime_error(
				"growth_factor requires z >= 0 and a flat LCDM cosmology");
	}

	if (z == 0.0)
		return 1.0;

	const double a = 1.0 / (1.0 + z);

	// Einstein-de Sitter case.
	if (Ol == 0.0)
		return a;

	// Integrand of the growth-factor expression.
	auto integrand = [](double x, void *params) -> double {

		const double *omega = static_cast<double *>(params);

		const double Om = omega[0];
		const double Ol = omega[1];

		return x * std::sqrt(x) /
				std::pow(Om + Ol * x * x * x, 1.5);
	};

	double omega[2] = {Om, Ol};

	gsl_function F;
	F.function = integrand;
	F.params = omega;

	gsl_integration_workspace *workspace =
			gsl_integration_workspace_alloc(1000);

	auto integral = [&](double upper) {

		double result, error;

		gsl_integration_qag(
			&F, 0.0, upper,
			0.0, 1.e-7,
			1000, 6,
			workspace,
			&result, &error);

		return result;
	};

	const double I_a = integral(a);
	const double I_1 = integral(1.0);

	gsl_integration_workspace_free(workspace);

	const double E_a = std::sqrt(Om / (a * a * a) + Ol);

	return E_a * I_a / I_1;
}

double Cosmology::growth_factor(double z) const
{
    if (!std::isfinite(z) || z < 0.0) {
        throw std::runtime_error(
            "growth_factor requires a finite, non-negative redshift");
    }

    if (z == 0.0)
        return 1.0;

    // Convert redshift to logarithmic scale factor.
    const double loga = -std::log10(1.0 + z);

    // Outside the lookup table, use direct integration.
    if (loga < log_a_grid.front() ||
        loga > log_a_grid.back()) {
        return growth_factor_direct(z);
    }

    // Uniform logarithmic scale-factor grid.
    const double dloga = log_a_grid[1] - log_a_grid[0];

    const std::size_t i = std::min(
        static_cast<std::size_t>(
            (loga - log_a_grid.front()) / dloga),
        log_a_grid.size() - 2
    );

    // Linear interpolation in log(a)-log(D).
    const double t =
        (loga - log_a_grid[i]) /
        (log_a_grid[i + 1] - log_a_grid[i]);

    const double logD =
        (1.0 - t) * log_growth_grid[i] +
        t * log_growth_grid[i + 1];

    return std::pow(10.0, logD);
}

double Cosmology::sigma_mass_direct(double mass) const
{
	/*TODO: check that the maximum tabulated wavenumber captures well
	 * the whole range of k_max*R required for Shark (especially at high-z).*/

	/*
	 * Calculate the RMS linear density fluctuation sigma(M)
	 * using the linear matter power spectrum at z=0.
	 *
	 * Units:
	 *   mass : Msun/h
	 *   k    : h/Mpc
	 *   P(k) : (Mpc/h)^3
	 *   R    : Mpc/h
	 */

	if (mass <= 0.0) {
		throw std::runtime_error("sigma_mass_direct requires a positive mass");
	}

	const auto &k = parameters.power_spectrum.k;
	const auto &pk = parameters.power_spectrum.p;

	if (k.size() < 2 || k.size() != pk.size()) {
		throw std::runtime_error("Invalid power spectrum table in sigma_mass_direct");
	}

	// Mean matter density in (Msun/h)/(Mpc/h)^3.
	const double rho_m = 2.7754e11 * parameters.OmegaM;

	// Top-hat smoothing radius in Mpc/h.
	const double R = std::cbrt(
			3.0 * mass / (4.0 * M_PI * rho_m));

	// Fourier transform of a spherical top-hat window.
	auto window = [](double x) {

		if (std::abs(x) < 0.01) {
			const double x2 = x * x;
			return 1.0 - x2 / 10.0
				   + x2 * x2 / 280.0;
		}

		return 3.0 *
			(std::sin(x) - x * std::cos(x))
			/ (x * x * x);
	};

	double sigma2 = 0.0;

	// Integrate over ln(k) using the trapezoidal rule.
	for (std::size_t i = 0; i < k.size() - 1; ++i) {

		const double k1 = k[i];
		const double k2 = k[i + 1];

		const double W1 = window(k1 * R);
		const double W2 = window(k2 * R);

		const double f1 = k1 * k1 * k1 * pk[i] * W1 * W1;
		const double f2 = k2 * k2 * k2 * pk[i + 1] * W2 * W2;

		sigma2 += 0.5 * (f1 + f2) * std::log(k2 / k1);
	}

	sigma2 /= (2.0 * M_PI * M_PI);

	return std::sqrt(sigma2);
}

double Cosmology::sigma_mass(double mass) const
{
    if (!std::isfinite(mass) || mass <= 0.0) {
        throw std::runtime_error(
            "sigma_mass requires a positive, finite mass");
    }

    const double logm = std::log10(mass);

    // Outside the lookup table, use direct integration.
    if (logm < log_mass_grid.front() ||
        logm > log_mass_grid.back()) {
        return sigma_mass_direct(mass);
    }

    // Uniform logarithmic mass grid.
    const double dlogm = log_mass_grid[1] - log_mass_grid[0];

    const std::size_t i = std::min(
        static_cast<std::size_t>(
            (logm - log_mass_grid.front()) / dlogm),
        log_mass_grid.size() - 2
    );

    // Linear interpolation in log(M)-log(sigma).
    const double t =
        (logm - log_mass_grid[i]) /
        (log_mass_grid[i + 1] - log_mass_grid[i]);

    const double logsigma =
        (1.0 - t) * log_sigma_grid[i] +
        t * log_sigma_grid[i + 1];

    return std::pow(10.0, logsigma);
}

} // namespace shark
