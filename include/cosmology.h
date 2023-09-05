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
 *
 * Cosmological parameters used as inputs for SHArk
 */

#ifndef SHARK_COSMOLOGY_H_
#define SHARK_COSMOLOGY_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

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
class CosmologicalParameters {

public:
	explicit CosmologicalParameters(const Options &options);

	float OmegaM = 0;
	float OmegaB = 0;
	float OmegaL = 0;
	float n_s = 0;
	float sigma8 = 0;
	float Hubble_h = 0;
	PowerSpectrumTable power_spectrum;

private:
	void load_tables(const std::string &power_spec_file);
};


/**
 * Cosmology class that will contain all cosmological parameters.
 */
class Cosmology {

public:
	explicit Cosmology(CosmologicalParameters parameters);

	double comoving_to_physical_angularmomentum(double r, double z) const;
	double comoving_to_physical_size(double r, double z) const;
	double comoving_to_physical_velocity(double v, double z) const;
	double comoving_to_physical_mass(double m) const;
	double physical_to_comoving_mass(double m) const;
	double physical_to_comoving_size(double r, double z) const;
	double convert_redshift_to_age(double z) const;
	double convert_age_to_redshift_lcdm(double t) const;
	double expansion_factor(double t) const;

	/**
	 * universal_baryon_fraction: calculates the baryon density with respect to the total matter density.
	 * @return
	 */
	double universal_baryon_fraction() const {
		return parameters.OmegaB/parameters.OmegaM;
	};

	double universal_baryon_fraction_relative_to_dm() const {
		return parameters.OmegaB/(parameters.OmegaM - parameters.OmegaB);
	};

	double hubble_parameter (double z) const;
	double critical_density (double z) const;

	CosmologicalParameters parameters;

};

/// Type to be used by users handling pointers to this class
using CosmologyPtr = std::shared_ptr<Cosmology>;

template <typename ...Ts>
CosmologyPtr make_cosmology(Ts&&...ts)
{
	return std::make_shared<Cosmology>(std::forward<Ts>(ts)...);
}

}  // namespace shark

#endif // SHARK_COSMOLOGY_H_
