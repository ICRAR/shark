/*
 * dark_matter_halos.h
 *
 *  Created on: 4Aug.,2017
 *      Author: clagos
 */

#ifndef INCLUDE_DARK_MATTER_HALOS_H_
#define INCLUDE_DARK_MATTER_HALOS_H_


#include <memory>
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

	DarkMatterProfile haloprofile;

};


class DarkMatterHalos {

public:
	DarkMatterHalos(std::shared_ptr<Cosmology> cosmology, SimulationParameters &sim_params);
	virtual ~DarkMatterHalos() {};

	virtual double grav_potential_halo(double r, double c) const = 0;

	double energy_circular (double r, double c);

	virtual double enclosed_mass(double r, double c) const = 0;

	double halo_dynamical_time (HaloPtr &halo, double z);

	double halo_virial_radius(HaloPtr &halo);

	double halo_virial_velocity (double mvir, double redshift);

	double halo_lambda (xyz<float> L, double mvir, double redshift);

	double disk_size_theory (Subhalo &subhalo);

	double halo_concentration (HaloPtr &halo);

	double gao_nfw_concentration(double mvir, double a_form);

	double mmw98_nfw_concentration(double mvir, double vmax, double rvir);

	void galaxy_velocity(Subhalo &subhalo);

	double v2halo (double x, double m, double c, double r);
	double v2disk (double x, double m, double c, double r);
	double v2bulge (double x, double m, double c, double r);

protected:
	std::shared_ptr<Cosmology> cosmology;
	SimulationParameters sim_params;

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
