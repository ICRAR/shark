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
 */

#include <cmath>
#include <fstream>
#include <map>
#include <stdexcept>
#include <tuple>

#include "disk_instability.h"
#include "galaxy.h"
#include "halo.h"
#include "numerical_constants.h"
#include "subhalo.h"

namespace shark {

DiskInstabilityParameters::DiskInstabilityParameters(const Options &options)
{
	options.load("disk_instability.stable", stable, true);
	options.load("disk_instability.fint", fint);
}

DiskInstability::DiskInstability(DiskInstabilityParameters parameters,
		GalaxyMergerParameters merger_params,
		SimulationParameters simparams,
		DarkMatterHalosPtr darkmatterhalo,
		std::shared_ptr<BasicPhysicalModel> physicalmodel,
		AGNFeedbackPtr agnfeedback) :
	parameters(parameters),
	merger_params(merger_params),
	simparams(std::move(simparams)),
	darkmatterhalo(std::move(darkmatterhalo)),
	physicalmodel(std::move(physicalmodel)),
	agnfeedback(std::move(agnfeedback))
{
	// no-op
}

void DiskInstability::evaluate_disk_instability (HaloPtr &halo, int snapshot, double delta_t){

	double z = simparams.redshifts[snapshot];

	for (auto &subhalo: halo->all_subhalos()){
		for (auto &galaxy: subhalo->galaxies){
			double f = toomre_parameter(galaxy);
			if(f < parameters.stable){
				/**
 				* Count number of disk instability episodes.
		 		*/ 
				galaxy.interaction.disk_instabilities += 1;

				/**
				 * Estimate new bulge size.
				 */
				galaxy.bulge_gas.rscale = bulge_size(galaxy);
				galaxy.bulge_stars.rscale = galaxy.bulge_gas.rscale;

				/**
				 * Transfer all stars and gas to the bulge.
				 */
				galaxy.bulge_stars.mass += galaxy.disk_stars.mass;
				galaxy.bulge_stars.mass_metals += galaxy.disk_stars.mass_metals;
				galaxy.bulge_gas.mass += galaxy.disk_gas.mass;
				galaxy.bulge_gas.mass_metals +=  galaxy.disk_gas.mass_metals;

				// Keep track of bulge mass that comes from the transfer of stars from the disk to the bulge.
				galaxy.diskinstabilities_assembly_stars.mass += galaxy.disk_stars.mass;
				galaxy.diskinstabilities_assembly_stars.mass_metals += galaxy.disk_stars.mass_metals;

				/**Assume both stars and gas mix up well during mergers.
				 * And calculate a pseudo specific AM as in mergers.*/
				//
				//calculate bulge specific angular momentum based on assuming conservation.
				//
				//effective_angular_momentum(galaxy);
				if(galaxy.bulge_mass() > 0){
					double v_pseudo = std::sqrt(constants::G * galaxy.bulge_mass() / galaxy.bulge_gas.rscale);
					galaxy.bulge_gas.sAM   = galaxy.bulge_gas.rscale * v_pseudo;
					galaxy.bulge_stars.sAM = galaxy.bulge_gas.sAM;
				}

				//Make all disk values 0.
				galaxy.disk_stars.restore_baryon();
				galaxy.disk_gas.restore_baryon();

				transfer_history_disk_to_bulge(galaxy, snapshot);

				create_starburst(subhalo, galaxy, z, snapshot, delta_t);
			}
		}
	}

}

double DiskInstability::toomre_parameter(const Galaxy &galaxy) const
{
	double vc = galaxy.vmax;
	double md =  galaxy.disk_mass();
	//double rd = galaxy->disk_gas.rscale;
	double rd = galaxy.disk_size();

	if(md <= 0 || rd <= 0){
		return 100;
	}

	double denom = 1.68 * constants::G * md / rd;

	double t = vc / std::sqrt(denom);

	return t;
}

double DiskInstability::bulge_size(const Galaxy &galaxy) const
{
	double md = galaxy.disk_mass();
	double mb = galaxy.bulge_mass();

	double rd = 0, rb = 0;
	if(md > 0){
		rd = (galaxy.disk_stars.rscale * galaxy.disk_stars.mass + galaxy.disk_gas.rscale * galaxy.disk_gas.mass) / md;
	}
	if(mb > 0){
		rb = galaxy.bulge_gas.rscale;
	}
	double c = merger_params.cgal;

	double bc = 0;
	if(mb > 0 && rb > 0){
		bc = c * std::pow(mb,2.0)/rb;
	}

	double dc = c * std::pow(md,2.0)/rd;

	double combined_c = parameters.fint * (md * mb)/(rd + rb);

	double rnew = c * std::pow((md + mb),2.0) / (bc + dc + combined_c);


	if(std::isnan(rnew) || rnew <= 0 || rnew >= 3){
		std::ostringstream os;
		os << galaxy << " has a bulge size not well defined in disk instabilities: " << rnew;
		throw invalid_data(os.str());
	}

	if(rnew <= constants::EPS6){
		std::ostringstream os;
		os << "Galaxy with extremely small size, rbulge_gas < 1e-6, in disk instabilities";
		//throw invalid_argument(os.str());
	}

	return rnew;

}

void DiskInstability::create_starburst(SubhaloPtr &subhalo, Galaxy &galaxy, double z, int snapshot, double delta_t){

	// Trigger starburst only in case there is gas in the bulge.
	if(galaxy.bulge_gas.mass > merger_params.mass_min){

		// Calculate black hole growth due to starburst.
		double tdyn = agnfeedback->smbh_accretion_timescale(galaxy, z);
		double delta_mbh = agnfeedback->smbh_growth_starburst(galaxy.bulge_gas.mass, subhalo->Vvir, tdyn, galaxy);
		double delta_mzbh = 0;
		if(galaxy.bulge_gas.mass > 0){
			delta_mzbh = delta_mbh/galaxy.bulge_gas.mass * galaxy.bulge_gas.mass_metals;
		}


		// Define accretion rate.
		galaxy.smbh.macc_sb += delta_mbh/tdyn;

		// Grow SMBH.
		galaxy.smbh.mass += delta_mbh;
		galaxy.smbh.mass_metals += delta_mzbh;

		// Reduce gas available for star formation due to black hole growth.
		galaxy.bulge_gas.mass -= delta_mbh;
		galaxy.bulge_gas.mass_metals -= delta_mzbh;

		// Trigger starburst.
		physicalmodel->evolve_galaxy_starburst(*subhalo, galaxy, z, delta_t, false);

		// Check for small gas reservoirs left in the bulge.
		if(galaxy.bulge_gas.mass > 0 && galaxy.bulge_gas.mass < merger_params.mass_min){

			galaxy.disk_gas        += galaxy.bulge_gas;

			if(galaxy.disk_gas.rscale == 0){
				galaxy.disk_gas.rscale = galaxy.bulge_gas.rscale;
				galaxy.disk_gas.sAM    = galaxy.bulge_gas.sAM;

				if (std::isnan(galaxy.disk_gas.sAM) || std::isnan(galaxy.disk_gas.rscale)) {
					throw invalid_argument("rgas or sAM are NaN, cannot continue at disk_instabilities");
				}
			}

			galaxy.bulge_gas.restore_baryon();

		}
	}
}

void DiskInstability::transfer_history_disk_to_bulge(Galaxy &galaxy, int snapshot){

	/**
	 * Function transfers the disk stellar mass history to bulge of the central galaxy.
	 */

	//Transfer history of stellar mass growth until the previous snapshot.
	for(int s=simparams.min_snapshot; s <= snapshot-1; s++) {

		auto it = std::find_if(galaxy.history.begin(), galaxy.history.end(), [s](const HistoryItem &hitem) {
			return hitem.snapshot == s;
		});

		if (it == galaxy.history.end()){ //galaxy didn't exist.
			//no-opt.
		}
		else { // both galaxies exist at this snapshot
			auto &hist = *it;

			//transfer disk information to bulge formed via disk instabilites
			hist.sfr_bulge_diskins   += hist.sfr_disk;
			hist.sfr_z_bulge_diskins += hist.sfr_z_disk;

			//make disk properties = 0;
			hist.sfr_disk   = 0;
			hist.sfr_z_disk = 0;
		}
	}

}

void DiskInstability::effective_angular_momentum(Galaxy &galaxy){

	double AM_disk_stars  = galaxy.disk_stars.angular_momentum();
	double AM_bulge_stars = galaxy.bulge_stars.angular_momentum();
	double AM_disk_gas    = galaxy.disk_gas.angular_momentum();
	double AM_bulge_gas   = galaxy.bulge_gas.angular_momentum();

	double mgas  = galaxy.gas_mass();
	double mstar = galaxy.stellar_mass();

	// Calculate effective specific angular momentum only if masses are > 0. Assume gas and stars mix up well.
	if(mgas + mstar > 0){
		double j = (AM_disk_stars + AM_bulge_stars + AM_disk_gas + AM_bulge_gas) / (mgas + mstar);
		galaxy.bulge_gas.sAM   = j;
		galaxy.bulge_stars.sAM = j;
	}

}


}//end namespace shark
