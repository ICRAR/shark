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
 * Star formation classes
 */

#ifndef INCLUDE_STAR_FORMATION_H_
#define INCLUDE_STAR_FORMATION_H_

#include <memory>

#include "cosmology.h"
#include "integrator.h"
#include "options.h"
#include "recycling.h"

namespace shark {

class StarFormationParameters {

public:
	explicit StarFormationParameters(const Options &options);

	/**Parameters:
	 * model: SF model applied
	 * nu_sf: inverse of the molecular gas depletion timescale for normal star-forming disks [Gyr^-1].
	 * Po: pressure normalization in the Blitz & Rosolowsky (2006) SF law [K cm^3].
	 * beta_press: pressure power-law index in the Blitz & Rosolowsky (2006) SF law.
	 * Accuracy_SFeqs: accuracy desired when numerically solving the integration of the SFR surface density.
	 * gas_velocity_dispersion: gas velocity dispersion [km/s].
	 * boost_starburst: boosting factor applied to the SF law in starbursts.
	 * sigma_HI_crit: critical surface density below which gas is assumed to be ionised.
	 * clump_factor_KMT09: clumping factor of the ISM for the Krumholz+ models.
	 * sigma_crit_KMT09: critical gas surface density above which the SF law becomes superlinear in the KMT09 model.
	 * angular_momentum_transfer: boolean parameter indicating whether the user wants to trigger the calculation of angular momentum transfer within the disk.
	 *
	 */
	enum StarFormationModel {
		BR06 = 0,//!< BR06
		GD14,    //!< GD14
		K13,     //!< K13
		KMT09,   //!< KMT09
		KD12     //!< KD12
	};

	StarFormationModel model = BR06;

	double nu_sf = 0;
	double Po = 0;
	double beta_press = 0;
	double Accuracy_SFeqs = 0.05;
	double gas_velocity_dispersion = 0;
	double boost_starburst = 1;
	double sigma_HI_crit = 0;
	double clump_factor_KMT09 = 1;
	double sigma_crit_KMT09 = 0;
	double efficiency_sf = 1;
        double gmc_surface_density = 85.0; //in Msun/pc^2

	bool angular_momentum_transfer = false;
};


class StarFormation {

public:

	struct molecular_gas {
		double m_mol;
		double m_atom;
		double m_mol_b;
		double m_atom_b;
		double j_mol;
		double j_atom;
	};

	StarFormation(StarFormationParameters parameters, RecyclingParameters recycleparams, CosmologyPtr cosmology);

	using func_t = double (*)(double x, void *);

	/**
	 * All input quantities should be in comoving units.
	 */
	double star_formation_rate(double mcold, double mstars, double rgas, double rstars, double zgas, double z,
							   bool burst, double vgal, double &jrate, double jgas);

	double star_formation_rate_surface_density(double r, void * params) const;

	double manual_integral(func_t f, void * params, double rmin, double rmax) const;

	double fmol(double Sigma_gas, double Sigma_stars, double zgas, double r) const;

	double midplane_pressure(double Sigma_gas, double Sigma_stars, double r) const;

	double gd14_sigma_norm(double d_mw, double u_mw) const;

	double kmt09_fmol(double zgas, double sigma_gas) const;

	double k13_fmol(double zgas, double sigma_gas) const;

	double kd12_taudep(double sigma_gas, void * params) const;

	std::size_t get_integration_intervals() {
		return integrator.get_num_intervals();
	}

	void reset_integration_intervals() {
		return integrator.reset_num_intervals();
	}

	double molecular_hydrogen(double mcold, double mstars, double rgas, double rstars, double zgas, double z, double &jmol,  double jgas, double vgal, bool bulge, bool jcalc);

	double molecular_surface_density(double r, void * params) const;

	molecular_gas get_molecular_gas(const Galaxy &galaxy, double z, bool jcalc);

	double ionised_gas_fraction(double mgas, double rgas, double z) const;

private:
	StarFormationParameters parameters;
	RecyclingParameters recycleparams;
	CosmologyPtr cosmology;
	Integrator integrator;

};

/// A collection of galaxy-indexed molecular gas objects
using molgas_per_galaxy = std::map<Galaxy::id_t, StarFormation::molecular_gas>;

}  // namespace shark

#endif /* INCLUDE_STAR_FORMATION_H_ */

