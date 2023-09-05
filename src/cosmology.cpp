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
	else if (lvalue == "wdm_dove") {
		return power_spectrum_type::WDM_DOVE;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value << ". ";
	os << "Supported values are millenium, millenium_gas, planck14, planck15 and wdm_dove";
	throw invalid_option(os.str());
}

static
std::map<power_spectrum_type, std::string> power_spectrum_files {
	{MILLENIUM, "pk_Mill.dat"},
	{MILLENIUM_GAS, "pk_MillGas_norm.dat"},
	{PLANCK14, "pk_Planck_norm.dat"},
	{PLANCK15, "pk_Planck_SURFS_norm.dat"},
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
	// no-op
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

} // namespace shark
