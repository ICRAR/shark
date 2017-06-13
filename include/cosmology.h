//
// Cosmological parameters used as inputs for SHArk
//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2017
// Copyright by UWA (in the framework of the ICRAR)
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307  USA
//
/*
 * Cosmology.h
 *
 *  Created on: 10Apr.,2017
 *      Author: clagos
 */

#ifndef SHARK_COSMOLOGY_H_
#define SHARK_COSMOLOGY_H_

#include <vector>
#include <string>

#include "options.h"

namespace shark {

/**
 * structure to save power spectrum
 */

struct PowerSpectrumTable {
	std::vector<double> k;
	std::vector<double> p;
};


/**
 * A set of cosmological parameters
 */
class CosmologicalParameters : public Options {

public:
	CosmologicalParameters(const std::string &filename);

	typedef std::map<double, std::string> tables_idx;

	float OmegaM;
	float OmegaB;
	float OmegaL;
	float n_s;
	float sigma8;
	float Hubble_h;
	PowerSpectrumTable power_spectrum;

private:
	void load_tables(const std::string &power_spec_file);
};


/**
 * Cosmology class that will contain all cosmological parameters.
 */
class Cosmology {

public:
	Cosmology(CosmologicalParameters parameters);

	double comoving_to_physical_size(double r, double z);
	double comoving_to_physical_velocity(double v, double z);
	double comoving_to_physical_mass(double m);
	double physical_to_comoving_mass(double m);

private:
	CosmologicalParameters parameters;

};

}  // namespace shark

#endif // SHARK_COSMOLOGY_H_
