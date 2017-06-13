/*
 * cosmology.c
 *
 *  Created on: 13Jun.,2017
 *      Author: clagos
 */


#include <cmath>
#include <fstream>
#include <map>
#include <tuple>

#include "cosmology.h"
#include "logging.h"
#include "numerical_constants.h"
#include "components.h"

namespace shark {

CosmologicalParameters::CosmologicalParameters(const std::string &filename) :
	Options(filename),
	OmegaM(0),
	OmegaB(0),
	OmegaL(0),
	n_s(0),
	sigma8(0),
	Hubble_h(0),
	power_spectrum()
{

	std::string power_spec_file;

	load("cosmology.OmegaM", OmegaM);
	load("cosmology.OmegaB", OmegaB);
	load("cosmology.OmegaL", OmegaL);
	load("cosmology.n_s", n_s);
	load("cosmology.sigma8", sigma8);
	load("cosmology.Hubble_h", Hubble_h);
	load("cosmology.power_spectrum_file", power_spec_file, true);

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
		if (is_skipable(line)) {
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

Cosmology::Cosmology(CosmologicalParameters) :
	parameters(parameters)
{
	// no-op
}

double Cosmology::comoving_to_physical_size(double r, double z){
	return r/parameters.Hubble_h/(1+z);
}

double Cosmology::comoving_to_physical_velocity(double v, double z){
	return v/(1+z);
}

double Cosmology::comoving_to_physical_mass(double m){
	return m/parameters.Hubble_h;
}

double Cosmology::physical_to_comoving_mass(double m){
	return m*parameters.Hubble_h;
}

}
