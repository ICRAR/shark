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

#ifndef INCLUDE_DISK_INSTABILITY_H_
#define INCLUDE_DISK_INSTABILITY_H_

#include "agn_feedback.h"
#include "components.h"
#include "galaxy_mergers.h"
#include "physical_model.h"
#include "simulation.h"

namespace shark {

class DiskInstabilityParameters {

public:
	explicit DiskInstabilityParameters(const Options &options);

	float stable = 0;
	float fint = 2;
};

class DiskInstability{

public:
	DiskInstability (DiskInstabilityParameters parameters,
			GalaxyMergerParameters merger_params,
			SimulationParameters simparams,
			DarkMatterHalosPtr darkmatterhalo,
			std::shared_ptr<BasicPhysicalModel> physicalmodel,
			AGNFeedbackPtr agnfeedback);

	double bulge_size(const Galaxy &galaxy) const;

	double toomre_parameter(const Galaxy &galaxy) const;

	void evaluate_disk_instability (HaloPtr &halo, int snapshot, double delta_t);

	void create_starburst(SubhaloPtr &subhalo, Galaxy &galaxy, double z, int snapshot, double delta_t);

	void transfer_history_disk_to_bulge(Galaxy &galaxy, int snapshot);

	void effective_angular_momentum(Galaxy &galaxy);

private:
	DiskInstabilityParameters parameters;
	GalaxyMergerParameters merger_params;
	SimulationParameters simparams;
	DarkMatterHalosPtr darkmatterhalo;
	std::shared_ptr<BasicPhysicalModel> physicalmodel;
	AGNFeedbackPtr agnfeedback;

};

} // namespace shark


#endif /* INCLUDE_DISK_INSTABILITY_H_ */
