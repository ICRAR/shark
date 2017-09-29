/*
 * cosmology.c
 *
 *  Created on: 13Jun.,2017
 *      Author: clagos
 */


#include <cmath>
#include <fstream>
#include <sstream>
#include <map>
#include <tuple>
#include <errno.h>
#include <stdexcept>

#include "cosmology.h"
#include "logging.h"
#include "numerical_constants.h"
#include "components.h"
#include "utils.h"

namespace shark {

CosmologicalParameters::CosmologicalParameters(const Options &options) :
	OmegaM(0),
	OmegaB(0),
	OmegaL(0),
	n_s(0),
	sigma8(0),
	Hubble_h(0),
	power_spectrum()
{

	std::string power_spec_file;

	options.load("cosmology.OmegaM", OmegaM);
	options.load("cosmology.OmegaB", OmegaB);
	options.load("cosmology.OmegaL", OmegaL);
	options.load("cosmology.n_s", n_s);
	options.load("cosmology.sigma8", sigma8);
	options.load("cosmology.Hubble_h", Hubble_h);
	options.load("cosmology.power_spectrum_file", power_spec_file, true);

	load_tables(power_spec_file);
}

void CosmologicalParameters::load_tables(const std::string &power_spec_file)
{

	using namespace std;

	LOG(debug) << "Reading table " << power_spec_file ;

	ifstream f = open_file(power_spec_file);
	string line;
	while ( getline(f, line) ) {

		trim(line);
		if (empty_or_comment(line)) {
			continue;
		}

		double k,p;
		istringstream iss(line);
		iss >> k >> p;

		power_spectrum.k.push_back(k);
		power_spectrum.p.push_back(p);

	}
	f.close();

}

Cosmology::Cosmology(CosmologicalParameters parameters) :
	parameters(parameters)
{
	// no-op
}

double Cosmology::comoving_to_physical_size(double r, double z){
	return r/parameters.Hubble_h;
}

double Cosmology::comoving_to_physical_velocity(double v, double z){
	return v;
}

double Cosmology::comoving_to_physical_mass(double m){
	return m/parameters.Hubble_h;
}

double Cosmology::physical_to_comoving_mass(double m){
	return m*parameters.Hubble_h;
}

double Cosmology::convert_redshift_to_age(double z){

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

double Cosmology::expansion_factor(double t){

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

double Cosmology::hubble_parameter (double z){
	double H2 = (parameters.OmegaM * std::pow(1.0 + z, 3.0) + parameters.OmegaL);
	return parameters.Hubble_h * 100.0 * std::sqrt(H2);
}

}
