/*
 * galaxy_mergers.cpp
 *
 *  Created on: 4Aug.,2017
 *      Author: clagos
 */

#include <cmath>
#include <random>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "agn_feedback.h"
#include "components.h"
#include "galaxy_mergers.h"
#include "numerical_constants.h"
#include "physical_model.h"

namespace shark {

GalaxyMergerParameters::GalaxyMergerParameters(const Options &options) :
	major_merger_ratio(0),
	minor_merger_burst_ratio(0),
	merger_random_seed(-1),
	jiang08(4),
	f_orbit(1),
	cgal(0.5)
	{

	options.load("galaxy_mergers.major_merger_ratio", major_merger_ratio, true);
	options.load("galaxy_mergers.minor_merger_burst_ratio", minor_merger_burst_ratio, true);
	options.load("galaxy_mergers.merger_random_seed", merger_random_seed);

	options.load("galaxy_mergers.jiang08_a", jiang08[0], true);
	options.load("galaxy_mergers.jiang08_b", jiang08[1], true);
	options.load("galaxy_mergers.jiang08_c", jiang08[2], true);
	options.load("galaxy_mergers.jiang08_d", jiang08[3], true);

	options.load("galaxy_mergers.f_orbit", f_orbit);
	options.load("galaxy_mergers.cgal", cgal);

}

GalaxyMergers::GalaxyMergers(GalaxyMergerParameters parameters, std::shared_ptr<DarkMatterHalos> darkmatterhalo, std::shared_ptr<BasicPhysicalModel> physicalmodel, std::shared_ptr<AGNFeedback> agnfeedback) :
	parameters(parameters),
	darkmatterhalo(darkmatterhalo),
	physicalmodel(physicalmodel),
	agnfeedback(agnfeedback)
{
	// no-op
}


void GalaxyMergers::orbital_parameters(double &vr, double &vt, double f){

	//double f2 = 1+1/std::max(mass_ratio,1+tolerance);

	//Now generate two random numbers between 0 and 3.
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0, 3);

	vr = distribution(generator);
	vt = distribution(generator);
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
	double r = 0;
	double x_hi = 1, x_lo = E;
	gsl_function F;

	// Structured passed as void * to GSL
	struct root_solver_pars {
		double c;
		std::shared_ptr<DarkMatterHalos> dark_matter_halo;
	};

	root_solver_pars pars {c, darkmatterhalo};
	F.function = [](double x, void * params) -> double {
		auto pars = static_cast<root_solver_pars *>(params);
		return pars->dark_matter_halo->energy_circular(x, pars->c);
	};
	F.params = &pars;

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

	/**
	 * Input variables:
	 * mp: mass  of primary galaxy.
	 * ms: mass of secondary galaxy.
	 */

	return 1+1/std::max(ms / mp, 1. + constants::tolerance);
}

double GalaxyMergers::merging_timescale_mass(double mp, double ms){

	/**
	 * Input variables:
	 * mp: mass  of primary galaxy.
	 * ms: mass of secondary galaxy.
	 */

	double mass_ratio = mp/ms;

	return mass_ratio/std::log(1+mass_ratio);
}

void GalaxyMergers::merging_timescale(SubhaloPtr &primary, SubhaloPtr &secondary){

	/**
	 * Function calculates the dynamical friction timescale for the subhalo secondary to merge into the subhalo primary.
	 * This should be calculated only in the snapshot before the secondary mergers onto the primary (i.e. disappears from merger tree).
	 * Inputs:
	 * a primary and secondary galaxy. The primary is the central galaxy.
	 */

	auto halo = primary->host_halo;

	double tau_dyn = darkmatterhalo->halo_dynamical_time(halo);

	double mp = primary->Mvir + primary->central_galaxy()->baryon_mass();

	for (auto &galaxy: secondary->galaxies){

		double mgal = galaxy->baryon_mass();
		double ms = secondary->Mvir + mgal;
		double tau_mass = merging_timescale_mass(mp, ms);

		// Define merging timescale and redefine type of galaxy.

		galaxy->tmerge = tau_mass * tau_dyn;
		galaxy->galaxy_type = Galaxy::TYPE2;

	}
//	double vt,vr;

//	double c = primary->host_halo->concentration;

//	double f = mass_ratio_function(mp, ms);

	//Calculate well the part of orbital parameters.
	//Draw orbital parameters from PDF in Benson et al. (2005).
	//orbital_parameters(vr, vt, f);

	//double tau_orbits = merging_timescale_orbital(vr, vt, f, c);

	//return tau_orbits * tau_mass * tau_dyn;


}

void GalaxyMergers::merging_subhalos(HaloPtr &halo){

	/**
	 * This function evaluates whether subhalos in each timestep are disappearing from the merger tree, and if they are
	 * it passes that satellite galaxy onto the central subhalo and calculates a dynamical friction timescale that the
	 * galaxy will save in case it's a satellite.
	 *
	 * Function receives as input a halo.
	 */

	auto central_subhalo = halo->central_subhalo;

	// Assign halo concentration.
	halo->concentration = halo->central_subhalo->concentration;

	for(auto &subhalo: halo->satellite_subhalos) {
		//Identify which subhalos will disappear in the next snapshot

		if(subhalo->last_snapshot_identified == subhalo->snapshot){

			auto satellite_subhalo = subhalo;

			//Calculate dynamical friction timescale for all galaxies in satellite_subhalo.
			merging_timescale(central_subhalo, satellite_subhalo);

			//transfer all mass from the satellite_subhalo to the central_subhalo.
			transfer_baryon_mass(satellite_subhalo, central_subhalo);

			//Now transfer the galaxies in this subhalo to the central subhalo.
			satellite_subhalo->transfer_galaxies_to(central_subhalo);
		}
	}

}

void GalaxyMergers::merging_galaxies(HaloPtr &halo, double z, double delta_t){

	/**
	 * This function determines which galaxies are merging in this snapshot by comparing tmerge with the duration of the snapshot.
	 * Inputs:
	 * halo: halo in which the two galaxies merging live.
	 * z: current redshift.
	 * delta_t: time interval between the current snapshots.
	 */

	//First define central subhalo.

	auto &central_subhalo = halo->central_subhalo;
	if (!central_subhalo) {
		std::ostringstream os;
		os << halo << " has no central subhalo";
		throw exception(os.str());
	}

	/**
	 * First find central galaxy of central subhalo.
	 */
	GalaxyPtr central_galaxy = central_subhalo->central_galaxy();
	if(!central_galaxy){
		std::ostringstream os;
		os << central_subhalo << " has no central galaxy";
		throw exception(os.str());
	}

	std::vector<GalaxyPtr> all_sats_to_delete;

	for (auto &galaxy: central_subhalo->galaxies){
		if(galaxy->galaxy_type == Galaxy::TYPE2){
			/**
			 * If merging timescale is less than the duration of this snapshot, then proceed to merge. Otherwise, update merging timescale.
			 */
			if(galaxy->tmerge < delta_t){
				create_merger(central_galaxy, galaxy, halo, z, delta_t);

				// Accummulate all satellites that we need to delete at the end.
				all_sats_to_delete.push_back(galaxy);
			}
			else{
				galaxy->tmerge = galaxy->tmerge - delta_t;
			}
		}
	}

	// Now destroy and remove satellite galaxy.
	// TODO: Maybe not most efficiently, but it will do for now
	for(auto &galaxy: all_sats_to_delete) {
		// we know "it" will point to something known because all_stas_to_delete
		// contain galaxies from central_subhalo->galaxies, which hasn't been
		// modified in the meanwhile. Thus, skip any check on "it" and use it blindly
		auto it = std::find(central_subhalo->galaxies.begin(), central_subhalo->galaxies.end(), galaxy);
		central_subhalo->galaxies.erase(it);
	}

}

void GalaxyMergers::create_merger(GalaxyPtr &central, GalaxyPtr &satellite, HaloPtr &halo, double z, double delta_t){

	/**
	 * This function classifies the merger and computes the starburst in case it takes place.
	 * Inputs:
	 * central: central galaxy.
	 * satellite: satellite galaxy.
	 * halo: halo in which the two galaxies merging live.
	 * z: current redshift.
	 * delta_t: time interval between the current snapshots.
	 */


	//First define central subhalo.

	auto &central_subhalo = halo->central_subhalo;

	double mbar_central = central->baryon_mass();

	double mbar_satellite = satellite->baryon_mass();

	double mass_ratio = mbar_satellite/mbar_central;

	//If mass ratio>1 is because satellite galaxy is more massive than central, so redefine mass ratio accordingly.
	if(mass_ratio > 1){
		mass_ratio = 1 / mass_ratio;
	}

	/**
	 * First, calculate remnant galaxy's bulge size based on merger properties.
	 */
	central->bulge_stars.rscale = bulge_size_merger(mass_ratio, central, satellite, halo);

	central->bulge_gas.rscale = central->bulge_stars.rscale;

	// Black holes merge regardless of the merger type.
	central->smbh.mass += satellite->smbh.mass;

	central->smbh.mass_metals += satellite->smbh.mass_metals;

	/**
	 * Evaluate major mergers
	 */
	if(mass_ratio >= parameters.major_merger_ratio){

		/**
		 * Transfer mass. In the case of major mergers, transfer all stars and gas to the bulge.
		 */

		central->bulge_stars.mass += central->disk_stars.mass + satellite->stellar_mass();
		central->bulge_stars.mass_metals += central->disk_stars.mass_metals + satellite->stellar_mass_metals();
		central->bulge_gas.mass += central->disk_gas.mass + satellite->gas_mass();
		central->bulge_gas.mass_metals +=  central->disk_gas.mass_metals + satellite->gas_mass_metals();

		// calculate black hole growth due to starburst.
		double delta_mbh = agnfeedback->smbh_growth_starburst(central->bulge_gas.mass);
		double delta_mzbh = delta_mbh/central->bulge_gas.mass * central->bulge_gas.mass_metals;

		// Grow SMBH.
		central->smbh.mass += delta_mbh;
		central->smbh.mass_metals += delta_mzbh;

		// Reduce gas available for star formation due to black hole growth.
		central->bulge_gas.mass -= delta_mbh;
		central->bulge_gas.mass_metals -= delta_mzbh;

		//Make all disk values 0.

		central->disk_stars.mass = 0;
		central->disk_stars.mass_metals = 0;
		central->disk_gas.mass = 0;
		central->disk_gas.mass_metals = 0;

		/**
		 * Triger starburst with available gas.
		 */
		physicalmodel->evolve_galaxy_starburst(*central_subhalo, *central, z, delta_t);

	}
	else{//minor mergers

		/**
		 * Transfer mass. In the case of minor mergers, transfer the satellite' stars to the central bulge, and the satellite' gas to the central's disk.
		 */

		central->bulge_stars.mass += satellite->stellar_mass();

		central->bulge_stars.mass_metals += satellite->stellar_mass_metals();

		central->disk_gas.mass += satellite->gas_mass();

		central->disk_gas.mass_metals +=  satellite->gas_mass_metals();

		if(mass_ratio >= parameters.minor_merger_burst_ratio){
			/**
			 * Triger starburst with available gas.
			 */

			central->bulge_gas.mass += central->disk_gas.mass;
			central->bulge_gas.mass_metals +=  central->disk_gas.mass_metals;

			// calculate black hole growth due to starburst.
			double delta_mbh = agnfeedback->smbh_growth_starburst(central->bulge_gas.mass);
			double delta_mzbh = delta_mbh/central->bulge_gas.mass * central->bulge_gas.mass_metals;

			// Grow SMBH.
			central->smbh.mass += delta_mbh;
			central->smbh.mass_metals += delta_mzbh;

			// Reduce gas available for star formation due to black hole growth.
			central->bulge_gas.mass -= delta_mbh;
			central->bulge_gas.mass_metals -= delta_mzbh;

			//Make gas disk values 0.
			central->disk_gas.mass = 0;
			central->disk_gas.mass_metals = 0;

			physicalmodel->evolve_galaxy_starburst(*central_subhalo, *central, z, delta_t);

		}

	}


}

double GalaxyMergers::bulge_size_merger(double mass_ratio, GalaxyPtr &central, GalaxyPtr &satellite, HaloPtr &halo){

	/**
	 * This function calculates the bulge sizes resulting from a galaxy mergers following Cole et al. (2000). This assumes
	 * that the internal energy of the remnant spheroid just after the mergers is equal to the sum of the internal and relative
	 * orbital energies of the two merging galaxies (neglecting any energy dissipation and mass loss during the merger).
	 * Inputs:
	 * central: central galaxy.
	 * satellite: satellite galaxy.
	 * halo: halo in which the two galaxies merging live.
	 */

	double mtotal_central = 0;
	double rcentral = 0;

	//Define central properties depending on whether merger is major or minor.
	if(mass_ratio >= parameters.major_merger_ratio){

 		double mbar_central = central->baryon_mass();

		rcentral = central->composite_size();

		double enc_mass = darkmatterhalo->enclosed_mass(rcentral/darkmatterhalo->halo_virial_radius(halo), halo->concentration);

		//Because central part of the DM halo behaves like the baryons, the mass of the central galaxy includes
		//the DM mass enclosed by rcentral.
		mtotal_central = mbar_central + halo->Mvir * enc_mass;
	}
	else{
		mtotal_central = central->bulge_mass();

		rcentral = central->bulge_stars.rscale;
	}

	double mbar_satellite = satellite->baryon_mass();

	double rsatellite = satellite->composite_size();

	return r_remnant(mtotal_central, mbar_satellite, rcentral, rsatellite);

}


double GalaxyMergers::r_remnant(double mc, double ms, double rc, double rs){

	/**
	 * Input variables:
	 * mc: mass central.
	 * ms: mass satellite.
	 * rc: radius central.
	 * rs: radius satellite.
	 */

	double factor1  = std::pow(mc,2)/rc;

	double factor2 = std::pow(ms,2)/rs;

	double factor3 = parameters.f_orbit/parameters.cgal *  mc * ms / (rc + rs);

	return std::pow((mc + ms),2)/ (factor1 + factor2 + factor3);
}

void GalaxyMergers::transfer_baryon_mass(SubhaloPtr satellite, SubhaloPtr central){

	central->hot_halo_gas.mass += satellite->hot_halo_gas.mass;
	central->hot_halo_gas.mass_metals += satellite->hot_halo_gas.mass_metals;

	central->cold_halo_gas.mass += satellite->cold_halo_gas.mass;
	central->cold_halo_gas.mass_metals += satellite->cold_halo_gas.mass_metals;

	central->ejected_galaxy_gas.mass += satellite->ejected_galaxy_gas.mass;
	central->ejected_galaxy_gas.mass_metals += satellite->ejected_galaxy_gas.mass_metals;

}



}  // namespace shark
