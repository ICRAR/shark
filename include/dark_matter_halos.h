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
#include <vector>

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
	DarkMatterHalos(DarkMatterHaloParameters &params, std::shared_ptr<Cosmology> cosmology, SimulationParameters &sim_params);
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

	void cooling_gas_sAM(Subhalo &subhalo, Galaxy &galaxy);

	void disk_sAM(Subhalo &subhalo, Galaxy &galaxy);

	void bulge_sAM(Subhalo &subhalo, Galaxy &galaxy);

	void transfer_bulge_am(SubhaloPtr &subhalo, GalaxyPtr &galaxy, double z);

	double v2halo (double x, double m, double c, double r);
	double v2disk (double x, double m, double c, double r);
	double v2bulge (double x, double m, double c, double r);

protected:
	DarkMatterHaloParameters params;
	std::shared_ptr<Cosmology> cosmology;
	SimulationParameters sim_params;
	std::default_random_engine generator;
	std::lognormal_distribution<double> distribution;

};

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

} // namespace shark




#endif /* INCLUDE_DARK_MATTER_HALOS_H_ */
