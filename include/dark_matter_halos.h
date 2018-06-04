/*
 * dark_matter_halos.h
 *
 *  Created on: 4Aug.,2017
 *      Author: clagos
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

class NfwDistribution {

private:
	std::uniform_real_distribution<double> uniform;
	const double a;
	const double norm;

public:
	NfwDistribution(const double c) :
	    uniform(0, 1),
	    a(1 / c),
	    norm(std::log((a + 1) / a) - 1 / (a + 1))
	{
	}

	template <typename G>
	double operator()(G &g) {
		auto p    = uniform(g);
		double m  = p * norm;
		double wm = - 1.0 / (std::exp( m + 1));
		double w  = gsl_sf_lambert_W0(wm);

		double x = a * (std::exp(w + m + 1) - 1);
		return x;
	}

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