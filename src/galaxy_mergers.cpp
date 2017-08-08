/*
 * galaxy_mergers.cpp
 *
 *  Created on: 4Aug.,2017
 *      Author: clagos
 */

#include <cmath>
#include <memory>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "components.h"
#include "galaxy_mergers.h"
#include "numerical_constants.h"

namespace shark {

GalaxyMergerParameters::GalaxyMergerParameters(const std::string &filename) :
	Options(filename),
	major_merger_ratio(0),
	minor_merger_burst_ratio(0),
	gas_fraction_minor_merger(0),
	merger_random_seed(-1),
	jiang08()
	{

	load("galaxy_mergers.major_merger_ratio", major_merger_ratio);
	load("galaxy_mergers.minor_merger_burst_ratio", minor_merger_burst_ratio);
	load("galaxy_mergers.gas_fraction_minorm_merger", gas_fraction_minor_merger);
	load("galaxy_mergers.merger_random_seed", merger_random_seed);
	load("galaxy_mergers.jiang08_a", jiang08[0]);
	load("galaxy_mergers.jiang08_b", jiang08[1]);
	load("galaxy_mergers.jiang08_c", jiang08[2]);
	load("galaxy_mergers.jiang08_d", jiang08[3]);

}

GalaxyMergers::GalaxyMergers(GalaxyMergerParameters parameters, std::shared_ptr<DarkMatterHalos> darkmatterhalo) :
	parameters(parameters),
	darkmatterhalo(darkmatterhalo)
{
	// no-op
}


void GalaxyMergers::orbital_parameters(double vr, double vt, double f){

	//
	using namespace constants;

	//double f2 = 1+1/std::max(mass_ratio,1+tolerance);

	//Now generate two random numbers between 0 and 3.
	std::default_random_engine generator;
	std::uniform_int_distribution<double> distribution(0,3);

	auto dice = std::bind ( distribution, generator );

	vr = dice();
	vt = dice();

}

double GalaxyMergers::merging_timescale_orbital(double vr, double vt, double f, double c){

	/**
	 * Input variables:
	 * (vr,vt): radial and tangential velocities of satellite galaxies.
	 * f: output of function mass_ratio_function.
	 * c: concentration parameter of halo.
	 */

	//Calculate energy of the orbit.

	double rVirial = 1; //Virial radius in scaled unit system.

	double E = 0.5 * (std::pow(vr,2)+std::pow(vt,2)) / f + darkmatterhalo->grav_potential_halo(rVirial , c);


	//TODO: check that I'm using the right input values for the merging timescale.

	int status;
	int iter = 0, max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r = 0, r_expected = sqrt (5.0);
	double x_hi = 1, x_lo = E;
	gsl_function F;
	struct quadratic_params params = {1.0, 0.0, -5.0};

	F.function = darkmatterhalo->energy_circular;
	F.params = &params;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc (T);
	gsl_root_fsolver_set (s, &F, x_lo, x_hi);

	do {
		iter++;
		status = gsl_root_fsolver_iterate (s);
		r = gsl_root_fsolver_root (s);
		x_lo = gsl_root_fsolver_x_lower (s);
		x_hi = gsl_root_fsolver_x_upper (s);
		status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);

		if (status == GSL_SUCCESS)
			printf ("Converged:\n");

	}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free (s);

	double rc_to_rvir = r;

	double vc = std::sqrt(f * darkmatterhalo->enclosed_mass(rc_to_rvir, c)) / rc_to_rvir;

	double eta = (vt/vc) / rc_to_rvir;

	//Apply Jiang et al. (2008) merging timescale.

	return std::sqrt(rc_to_rvir)*(parameters.jiang08[0] * std::pow(eta, parameters.jiang08[1]) + parameters.jiang08[3]) /2 / parameters.jiang08[2];
}

double GalaxyMergers::mass_ratio_function(double mp, double ms){

	return 1+1/std::max(ms/mp, 1+constants::tolerance);
}

double GalaxyMergers::merging_timescale_mass(double mp, double ms){

	double mass_ratio = mp/ms;

	return mass_ratio/std::log(1+mass_ratio);
}

double GalaxyMergers::merging_timescale(std::shared_ptr<Subhalo> &primary, std::shared_ptr<Subhalo> &secondary){

	/**
	 * Function calculates the dynamical friction timescale for the subhalo secondary to merge into the subhalo primary.
	 * This should be calculated only in the snapshot before the secondary mergers onto the primary (i.e. disappears from merger tree).
	 */
	double vt,vr;

	double ms = secondary->Mvir;

	double mp = primary->Mvir;

	double c = primary->host_halo->concentration;

	auto halo = primary->host_halo;

	double tau_dyn = darkmatterhalo->halo_dynamical_time(halo);

	double f = mass_ratio_function(mp, ms);

	//Draw orbital parameters from PDF in Benson et al. (2005).
	orbital_parameters(vr, vt, f);

	double tau_orbits = merging_timescale_orbital(vr, vt, f, c);

	double tau_mass = merging_timescale_mass(mp, ms);

	return tau_orbits * tau_mass * tau_dyn;

}

void GalaxyMergers::merging_subhalos(std::shared_ptr<Halo> &halo){

	/**
	 * This function evaluates whether subhalos in each timestep are disappearing from the merger tree, and if they are
	 * it passes that satellite galaxy onto the central subhalo and calculates a dynamical friction timescale that the
	 * galaxy will save in case it's a satellite.
	 *
	 * Function receives as input a halo.
	 */

	auto central_subhalo = halo->central_subhalo;

	for(std::shared_ptr<Subhalo> &subhalo: halo->satellite_subhalos) {
		//Identify which subhalos will disappear in the next snapshot

		if(subhalo->last_snapshot_identified == subhalo->snapshot){

			auto satellite_subhalo = subhalo;

			//Calculate dynamical friction timescale.
			double tau_fric = merging_timescale(central_subhalo, satellite_subhalo);

			//Assign tau_fric to all satellite galaxies in the subhalo that will disappear.
			for (std::shared_ptr<Galaxy> &galaxies: satellite_subhalo->galaxies){
				galaxies->tmerge = tau_fric;
			}

			//Now transfer the galaxies in this subhalo to the central subhalo.
			central_subhalo->galaxies.insert(central_subhalo->galaxies.end(), satellite_subhalo->galaxies.begin(), satellite_subhalo->galaxies.end());

			//Delete galaxies from satellite_subhalo
			satellite_subhalo->galaxies.clear();
		}
	}


}

void GalaxyMergers::merging_galaxies(std::shared_ptr<Halo> &halo){

	/**
	 * This function determines which galaxies are merging in this snapshot by comparing tmerge with the duration of the snapshot.
	 */

	//First define central subhalo.

	Subhalo central_subhalo = halo->central_subhalo;


}




}  // namespace shark
