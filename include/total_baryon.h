//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2019
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
 * Total baryon-related classes and functionality
 */

#ifndef INCLUDE_TOTAL_BARYON_H_
#define INCLUDE_TOTAL_BARYON_H_

#include <map>
#include <vector>

#include "baryon.h"

namespace shark {

/**
 * Class to track all the mass componets in the simulation.
 */
class TotalBaryon {

public:

	/**
	 * mcold: total cold gas mass in disk+bulge.
	 * mstars: total stellar mass in disk+bulge.
	 * mstars_burst_galaxymergers: total stellar mass formed via bursts driven by galaxy mergers.
	 * mstars_burst_diskinstabilities: total stellar mass formed via bursts driven by disk instabilities.
	 * mhot_halo: total hot halo gas.
	 * mcold_halo: total cold halo gas (that is cooling during the current snapshot).
	 * mejected_halo: total hot gas that has been ejected from galaxies due to stellar feedback and that has not been reincorporated yet onto the hot halo gas reservoir.
	 * mlost_halo: total gas mass that has been ejected from galaxies and halos due to QSO feedback.
	 * mBH: total mass locked in black holes.
	 * mHI: total mass in the form of atomic gas.
	 * mH2: total mass in the form of molecular gas.
	 * mDM: total mass in the form of dark matter.
	 * SFR: integrated SFR of all galaxies over a snapshot for the disk and the bulge.
	 * baryon_total_created: keeps track of the baryons deposited in DM halos to ensure mass convervations.
	 * max_BH: maximum mass of the SMBHs in this snapshot.
	 */

	std::vector<BaryonBase> mcold;
	std::vector<BaryonBase> mstars;
	std::vector<BaryonBase> mstars_burst_galaxymergers;
	std::vector<BaryonBase> mstars_burst_diskinstabilities;
	std::vector<BaryonBase> mhot_halo;
	std::vector<BaryonBase> mcold_halo;
	std::vector<BaryonBase> mejected_halo;
	std::vector<BaryonBase> mlost_halo;
	std::vector<BaryonBase> mBH;
	std::vector<BaryonBase> mHI;
	std::vector<BaryonBase> mH2;
	std::vector<BaryonBase> mDM;

	std::vector<double> SFR_disk;
	std::vector<double> SFR_bulge;
	std::vector<double> max_BH;

	/**
	 * Vectors of integers that keep track of number of mergers of galaxies and disk instability episodes in each snapshot.
	 * major_mergers: counts number of major mergers per snapshot.
	 * minor_mergers: counts number of minor mergers per snapshot.
	 * disk_instabil: counts number of disk instability episodes per snapshot.
	*/
	std::vector<int> major_mergers;
	std::vector<int> minor_mergers;
	std::vector<int> disk_instabil;

	std::map<int, double> baryon_total_created;
	std::map<int, double> baryon_total_lost;

	std::vector<double> get_masses (const std::vector<BaryonBase> &B) const;
	std::vector<double> get_metals (const std::vector<BaryonBase> &B) const;

};

}  // namespace shark

#endif /* INCLUDE_TOTAL_BARYON_H_ */
