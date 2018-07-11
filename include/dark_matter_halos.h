//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2018
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

/**
 * @file
 */

#ifndef INCLUDE_DARK_MATTER_HALOS_H_
#define INCLUDE_DARK_MATTER_HALOS_H_

#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include <gsl/gsl_sf_lambert.h>

#include "mixins.h"
#include "numerical_constants.h"
#include "components.h"
#include "cosmology.h"
#include "simulation.h"

namespace shark {

class DarkMatterHaloParameters {

public:
	DarkMatterHaloParameters(const Options &options);

	enum DarkMatterProfile {
		NFW = 0,
		EINASTO
	};

	enum SizeModel{
		MO98 = 0,
		COLE00
	};

	DarkMatterProfile haloprofile = NFW;
	SizeModel sizemodel = MO98;
	bool random_lambda = false;

};

class DarkMatterHalos {

public:
	DarkMatterHalos(const DarkMatterHaloParameters &params, const CosmologyPtr &cosmology, SimulationParameters &sim_params);
	virtual ~DarkMatterHalos() {};

	virtual double grav_potential_halo(double r, double c) const = 0;

	double energy_circular (double r, double c);

	virtual double enclosed_mass(double r, double c) const = 0;

	double halo_dynamical_time (HaloPtr &halo, double z);

	double halo_virial_radius(Subhalo &subhalo);

	double halo_virial_velocity (double mvir, double redshift);

	float halo_lambda (float lambda, double z);

	double disk_size_theory (Subhalo &subhalo, double z);

	double halo_concentration (HaloPtr &halo);

	double nfw_concentration(double mvir, double z);

	//double mmw98_nfw_concentration(double mvir, double vmax, double rvir);

	void cooling_gas_sAM(Subhalo &subhalo, double z);

	void disk_sAM(Subhalo &subhalo, Galaxy &galaxy);

	void bulge_sAM(Subhalo &subhalo, Galaxy &galaxy);

	void transfer_bulge_am(SubhaloPtr &subhalo, GalaxyPtr &galaxy, double z);

	double v2halo (double x, double m, double c, double r);
	double v2disk (double x, double m, double c, double r);
	double v2bulge (double x, double m, double c, double r);

	void generate_random_orbits(xyz<float> &pos, xyz<float> &v, xyz<float> &L, double total_am, const HaloPtr &halo);

protected:
	DarkMatterHaloParameters params;
	CosmologyPtr cosmology;
	SimulationParameters sim_params;
	std::default_random_engine generator;
	std::lognormal_distribution<double> distribution;
	std::uniform_real_distribution<double> flat_distribution;

};

/// Type used by users to keep track o
typedef std::shared_ptr<DarkMatterHalos> DarkMatterHalosPtr;

class NFWDarkMatterHalos : public DarkMatterHalos {

public:
	using DarkMatterHalos::DarkMatterHalos;

	double grav_potential_halo(double r, double c) const override;
	double enclosed_mass(double r, double c) const override;
};

class EinastoDarkMatterHalos : public DarkMatterHalos {

public:
	using DarkMatterHalos::DarkMatterHalos;

	double grav_potential_halo(double r, double c) const override;
	double enclosed_mass(double r, double c) const override;
};


/// Factory of DarkMatterHaloPtrs
template <typename ...Ts>
DarkMatterHalosPtr make_dark_matter_halos(const DarkMatterHaloParameters &dmh_parameters, Ts&&...ts)
{
	if (dmh_parameters.haloprofile == DarkMatterHaloParameters::NFW) {
		return std::make_shared<NFWDarkMatterHalos>(dmh_parameters, std::forward<Ts>(ts)...);
	}
	else if (dmh_parameters.haloprofile == DarkMatterHaloParameters::EINASTO) {
		return std::make_shared<EinastoDarkMatterHalos>(dmh_parameters, std::forward<Ts>(ts)...);
	}

	std::ostringstream os;
	os << "Dark Matter halo profile " << dmh_parameters.haloprofile
	   << " not currently supported";
	throw invalid_argument(os.str());
}

} // namespace shark


#endif /* INCLUDE_DARK_MATTER_HALOS_H_ */