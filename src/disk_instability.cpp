/*
 * disk_instability.cpp
 *
 *  Created on: 29Sep.,2017
 *      Author: clagos
 */

#include <cmath>
#include <fstream>
#include <map>
#include <stdexcept>
#include <tuple>

#include "components.h"
#include "disk_instability.h"
#include "numerical_constants.h"

namespace shark {

DiskInstabilityParameters::DiskInstabilityParameters(const Options &options) :
	stable(0),
	fint(0)
	{

	options.load("disk_instability.stable", stable);
	options.load("disk_instability.fint", fint);


}

DiskInstability::DiskInstability(DiskInstabilityParameters parameters,
		GalaxyMergerParameters merger_params,
		SimulationParameters simparams,
		std::shared_ptr<DarkMatterHalos> darkmatterhalo,
		std::shared_ptr<BasicPhysicalModel> physicalmodel,
		std::shared_ptr<AGNFeedback> agnfeedback) :
	parameters(parameters),
	merger_params(merger_params),
	simparams(simparams),
	darkmatterhalo(darkmatterhalo),
	physicalmodel(physicalmodel),
	agnfeedback(agnfeedback)
{
	// no-op
}

void DiskInstability::evaluate_disk_instability (HaloPtr &halo, int snapshot, double delta_t){


	double z = simparams.redshifts[snapshot];

	for (auto &subhalo: halo->all_subhalos()){
		for (auto &galaxy: subhalo->galaxies){
			double f = toomre_parameter(galaxy, subhalo);
			if(f < parameters.stable){

				/**
				 * Estimate new bulge size.
				 */
				galaxy->bulge_gas.rscale = bulge_size(galaxy);
				galaxy->bulge_stars.rscale = galaxy->bulge_gas.rscale;

				/**
				 * Transfer all stars and gas to the bulge.
				 */
				galaxy->bulge_stars.mass += galaxy->disk_stars.mass;
				galaxy->bulge_stars.mass_metals += galaxy->disk_stars.mass_metals;
				galaxy->bulge_gas.mass += galaxy->disk_gas.mass;
				galaxy->bulge_gas.mass_metals +=  galaxy->disk_gas.mass_metals;

				//Make all disk values 0.
				galaxy->disk_stars.mass = 0;
				galaxy->disk_stars.mass_metals = 0;
				galaxy->disk_gas.mass = 0;
				galaxy->disk_gas.mass_metals = 0;

				transfer_history_disk_to_bulge(galaxy, snapshot);

				/*
				galaxy->disk_gas.rscale = 0;
				galaxy->disk_stars.rscale = 0;
				galaxy->disk_gas.sAM = 0;
				galaxy->disk_stars.sAM = 0;*/

				create_starburst(subhalo, galaxy, z, delta_t);
			}
		}
	}

}

double DiskInstability::toomre_parameter(GalaxyPtr &galaxy, SubhaloPtr &subhalo){

	double vc = subhalo->Vcirc;
	double md =  galaxy->disk_mass();
	double rd = galaxy->disk_gas.rscale;

	if(md <= 0 or rd <= 0){
		return 100;
	}

	double denom = 1.68 * constants::G * md / rd;

	double t = vc / std::sqrt(denom);

	return t;
}

double DiskInstability::bulge_size(GalaxyPtr &galaxy){

	double md = galaxy->disk_mass();
	double mb = galaxy->bulge_mass();

	double rd = galaxy->disk_stars.rscale;
	double rb = galaxy->bulge_stars.rscale;

	double c = merger_params.cgal;

	double bc = 0;
	if(mb > 0 and rb > 0){
		bc = c * std::pow(mb,2.0)/rb;
	}

	double dc = c * std::pow(md,2.0)/rd;

	double combined_c = parameters.fint * (md * mb)/(rd + rb);

	double rnew = c * std::pow((md + mb),2.0) / (bc + dc + combined_c);


	if(std::isnan(rnew) or rnew <= 0 or rnew > 0.1){
		std::ostringstream os;
		os << galaxy << " has a bulge size not well defined in disk instabilities.";
		throw invalid_data(os.str());
	}


	return rnew;

}

void DiskInstability::create_starburst(SubhaloPtr &subhalo, GalaxyPtr &galaxy, double z, double delta_t){

	// Trigger starburst only in case there is gas in the bulge.
	if(galaxy->bulge_gas.mass > constants::tolerance_mass){

		// Calculate black hole growth due to starburst.
		double delta_mbh = agnfeedback->smbh_growth_starburst(galaxy->bulge_gas.mass, subhalo->Vvir);
		double delta_mzbh = 0;
		if(galaxy->bulge_gas.mass > 0){
			delta_mzbh = delta_mbh/galaxy->bulge_gas.mass * galaxy->bulge_gas.mass_metals;
		}

		double tdyn = agnfeedback->smbh_accretion_timescale(*galaxy, z);

		// Define accretion rate.
		galaxy->smbh.macc_sb = delta_mbh/tdyn;

		// Grow SMBH.
		galaxy->smbh.mass += delta_mbh;
		galaxy->smbh.mass_metals += delta_mzbh;

		// Reduce gas available for star formation due to black hole growth.
		galaxy->bulge_gas.mass -= delta_mbh;
		galaxy->bulge_gas.mass_metals -= delta_mzbh;

		// Trigger starburst.
		physicalmodel->evolve_galaxy_starburst(*subhalo, *galaxy, z, delta_t);

		// Check for small gas reservoirs left in the bulge.
		if(galaxy->bulge_gas.mass < constants::tolerance_mass){
			galaxy->disk_gas.mass += galaxy->bulge_gas.mass;
			galaxy->disk_gas.mass_metals += galaxy->bulge_gas.mass_metals;
			galaxy->bulge_gas.mass = 0;
			galaxy->bulge_gas.mass_metals = 0;

			/*// Calculate disk size.
			galaxy->disk_gas.rscale = darkmatterhalo->disk_size_theory(*subhalo);
			galaxy->disk_stars.rscale = galaxy->disk_gas.rscale;*/
		}
	}
}

void DiskInstability::transfer_history_disk_to_bulge(GalaxyPtr &galaxy, int snapshot){

	/**
	 * Function transfers the disk stellar mass history to bulge of the central galaxy.
	 */

	//Transfer history of stellar mass growth until the previous snapshot.
	for(int s=simparams.min_snapshot; s <= snapshot-1; s++) {

		auto it = std::find_if(galaxy->history.begin(), galaxy->history.end(), [s](const HistoryItem &hitem) {
			return hitem.snapshot == s;
		});

		if (it == galaxy->history.end()){ //galaxy didn't exist.
			//no-opt.
		}
		else { // both galaxies exist at this snapshot
			auto &hist = *it;

			//tranfer disk information to bulge.
			hist.sfr_bulge += hist.sfr_disk;
			hist.stellar_bulge.mass += hist.stellar_disk.mass;
			hist.stellar_bulge.mass_metals += hist.stellar_disk.mass_metals;

			//make disk properties = 0;
			hist.sfr_disk = 0;
			hist.stellar_disk.mass = 0;
			hist.stellar_disk.mass_metals= 0;
		}
	}

}


}//end namespace shark
