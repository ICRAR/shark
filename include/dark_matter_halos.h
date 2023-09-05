//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2018
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

#ifndef INCLUDE_DARK_MATTER_HALOS_H_
#define INCLUDE_DARK_MATTER_HALOS_H_

#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "mixins.h"
#include "numerical_constants.h"
#include "components.h"
#include "cosmology.h"
#include "simulation.h"
#include "execution.h"

namespace shark {

class DarkMatterHaloParameters {

public:
	explicit DarkMatterHaloParameters(const Options &options);

	enum DarkMatterProfile {
		NFW = 0,
		EINASTO
	};

	enum SizeModel{
		MO98 = 0,
		COLE00
	};

	enum ConcentrationModel{
        	DUFFY08 = 0,
		DUTTON14 
	};

	DarkMatterProfile haloprofile = NFW;
	SizeModel sizemodel = MO98;
        ConcentrationModel concentrationmodel = DUFFY08;

	/**
	Note that if random_lambda = true, the values of lambda will be drawn from a log-normal distribution randomly. However, if at the same time 
	use_converged_lambda_catalog = true, the latter will be done only for halos that have a number of particles below min_part_convergence.
	**/
	bool random_lambda = false;
	bool use_converged_lambda_catalog = false; 
	int  min_part_convergence = 100;

};

class DarkMatterHalos {

public:
	DarkMatterHalos(
		const DarkMatterHaloParameters &params,
		CosmologyPtr cosmology,
		SimulationParameters &sim_params,
		ExecutionParameters exec_params);
	virtual ~DarkMatterHalos() = default;

	virtual double grav_potential_halo(double r, double c) const = 0;

	double energy_circular (double r, double c);

	virtual double enclosed_mass(double r, double c) const = 0;

	double halo_dynamical_time (HaloPtr &halo, double z);

	double subhalo_dynamical_time (Subhalo &subhalo, double z);

	double halo_virial_radius(const HaloPtr &halo, double z);

	double halo_virial_velocity (double mvir, double redshift);

	float halo_lambda (const Subhalo &subhalo, float m, double z, double npart);

	double disk_size_theory (Subhalo &subhalo, double z);

	double halo_concentration (HaloPtr &halo);

	double nfw_concentration(double mvir, double z);

	//double mmw98_nfw_concentration(double mvir, double vmax, double rvir);

	void cooling_gas_sAM(Subhalo &subhalo, double z);

	float enclosed_total_mass(const Subhalo &subhalo, double z, float r);

	void disk_sAM(Subhalo &subhalo, Galaxy &galaxy, double z);

	void bulge_sAM(Subhalo &subhalo, Galaxy &galaxy, double z);

	void transfer_bulge_am(SubhaloPtr &subhalo, Galaxy &galaxy, double z);

	double v2halo (double x, double m, double c, double r);
	double v2disk (double x, double m, double c, double r);
	double v2bulge (double x, double m, double c, double r);

	void generate_random_orbits(xyz<float> &pos, xyz<float> &v, xyz<float> &L, double total_am, const HaloPtr &halo, const Galaxy &galaxy);

protected:
	DarkMatterHaloParameters params;
	CosmologyPtr cosmology;
	SimulationParameters sim_params;
	ExecutionParameters exec_params;

private:
	xyz<float> random_point_in_sphere(float r, std::default_random_engine &generator);
};

/// Type used by users to keep track o
using DarkMatterHalosPtr= std::shared_ptr<DarkMatterHalos>;

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
