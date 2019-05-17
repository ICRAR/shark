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
 * Galaxy-related classes and functionality
 */

#ifndef INCLUDE_GALAXY_H_
#define INCLUDE_GALAXY_H_

#include <ostream>
#include <vector>

#include "baryon.h"
#include "components.h"
#include "mixins.h"

namespace shark {

/**
 * Structure that saves the history of relevant baryon components
 * needed for SED calculation later on.
 */
struct HistoryItem {
	float sfr_disk;
	float sfr_bulge_mergers;
	float sfr_bulge_diskins;
	float sfr_z_disk;
	float sfr_z_bulge_mergers;
	float sfr_z_bulge_diskins;
	int snapshot;
};

// TODO: add documentation
struct InteractionItem {
	int major_mergers = 0;
	int minor_mergers = 0;
	int disk_instabilities = 0;

	// TODO: rename to simply "restore"
	void restore_interaction_item()
	{
		major_mergers = 0;
		minor_mergers = 0;
		disk_instabilities = 0;
	}
};

/**
 * A basic galaxy.
 *
 * Galaxies have at least a bulge, a disk and a black hole. The most basic
 * galaxy allowed in SHArk contains a disk and a bulge, with either component
 * having gas and stars. A Galaxy also requires to have a SMBH.
 */

class Galaxy : public Identifiable<int> {
public:

	using Identifiable::Identifiable;

	/**
	 * An enumeration of types of galaxies
	 */
	enum galaxy_type_t {
		CENTRAL = 0,
		TYPE1,
		TYPE2,
		FLYBY
	};


	/**
	 * The ID of the descendant of this galaxy.
	 */
	id_t descendant_id = -1;

	/**
	 * The type of galaxy
	 */
	galaxy_type_t galaxy_type = CENTRAL;

	Baryon bulge_stars;
	Baryon bulge_gas;
	Baryon disk_stars;
	Baryon disk_gas;
	Baryon galaxymergers_burst_stars;
	Baryon galaxymergers_assembly_stars;
	Baryon diskinstabilities_burst_stars;
	Baryon diskinstabilities_assembly_stars;
	BlackHole smbh;

	//save average star formation rates and metallicities of the newly formed stars.
	float sfr_disk = 0;
	float sfr_bulge_mergers = 0;
	float sfr_bulge_diskins = 0;
	float sfr_z_disk = 0;
	float sfr_z_bulge_mergers = 0;
	float sfr_z_bulge_diskins = 0;

	/**
	 * Keep track of mean stellar age using:
	 *  mean_stellar_age: stellar mass formed times the mean age at which they formed.
	 *  total_stellar_mass_ever_formed: total stellar mass ever formed (without including the effects of stellar populations).
	 */
	float mean_stellar_age = 0;
	float total_stellar_mass_ever_formed = 0;

	/// maximum circular velocity.
	float vmax = 0;

	/// star formation and gas history of this galaxy across snapshots
	std::vector<HistoryItem>  history;

	/// interactions of this galaxy during this snapshot
	InteractionItem interaction;

	/// dynamical friction timescale, defined only if galaxy is satellite.
	float tmerge = 0;
	/// concentration of the subhalo this galaxy was before becoming type 2, only relevant for type 2 galaxies
	float concentration_type2 = 0;
	/// subhalo mass of this galaxy before it became type 2, only relevant for type 2 galaxies
	float msubhalo_type2 = 0;
	/// subhalo virial velocity of this galaxy before it became type 2, only relevant for type 2 galaxies
	float vvir_type2 = 0;
	/// subhalo spin parameter of this galaxy before it became type 2, only relevant for type 2 galaxies
	float lambda_type2 = 0;

	/**
	 * Define functions to calculate total mass and metals of various components.
	 */

	double disk_mass() const
	{
		return disk_gas.mass + disk_stars.mass;
	}

	double disk_mass_metals() const
	{
		return disk_gas.mass_metals + disk_stars.mass_metals;
	}

	double bulge_mass() const
	{
		return bulge_gas.mass + bulge_stars.mass;
	}

	double bulge_mass_metals() const
	{
		return bulge_gas.mass_metals + bulge_stars.mass_metals;
	}

	double baryon_mass() const
	{
		return disk_gas.mass + disk_stars.mass + bulge_gas.mass + bulge_stars.mass;
	}

	double stellar_mass() const
	{
		return disk_stars.mass + bulge_stars.mass;
	}

	double stellar_mass_metals() const
	{
		return disk_stars.mass_metals + bulge_stars.mass_metals;
	}

	double gas_mass() const
	{
		return disk_gas.mass + bulge_gas.mass;
	}

	double gas_mass_metals() const
	{
		return disk_gas.mass_metals + bulge_gas.mass_metals;
	}

	double disk_size() const
	{
		double rgas  = 0;
		double rstar = 0;

		if(disk_gas.mass > 0){
			rgas = disk_gas.rscale;
		}
		if(disk_stars.mass > 0){
			rstar = disk_stars.rscale;
		}

		double rcomp = 0.0;

		// Define rcomp only if disk has mass.
		if(disk_mass() > 0){
			rcomp = (disk_stars.mass * rstar + disk_gas.mass * rgas) / disk_mass();
		}

		return rcomp;
	}

	double bulge_size() const
	{
		double rgas  = 0;
		double rstar = 0;

		if(bulge_gas.mass > 0){
			rgas = bulge_gas.rscale;
		}
		if(bulge_stars.mass > 0){
			rstar = bulge_stars.rscale;
		}

		double rcomp = 0.0;

		// Define rcomp only if bulge has mass.
		if(bulge_mass() > 0){
			rcomp = (bulge_stars.mass * rstar + bulge_gas.mass * rgas) / bulge_mass();
		}

		return rcomp;
	}

	double composite_size() const
	{
		double rdisk = disk_size();
		double rbulge = bulge_stars.rscale;

		double rcomp = 0.0;

		// Define rcomp only if galaxy has mass.
		if(baryon_mass() > 0){
			rcomp = (disk_mass() * rdisk + bulge_mass() * rbulge) / baryon_mass();
		}

		return rcomp;
	}

	double stellar_size() const
	{
		double rdisk = disk_size();
		double rbulge = bulge_size();

		double rcomp = 0.0;

		// Define rcomp only if galaxy has mass.
		if(baryon_mass() > 0){
			rcomp = (disk_stars.mass * rdisk + bulge_mass() * rbulge) / baryon_mass();
		}

		return rcomp;
	}

	double angular_momentum() const
	{
		return disk_gas.angular_momentum() + disk_stars.angular_momentum() + bulge_gas.angular_momentum() + bulge_stars.angular_momentum();
	}

};

template <typename T>
std::basic_ostream<T> &operator<<(std::basic_ostream<T> &stream, const Galaxy &galaxy)
{
	stream << "<Galaxy " << galaxy.id << ", type=" << galaxy.galaxy_type << '>';
	return stream;
}

template <typename T>
std::basic_ostream<T> &operator<<(std::basic_ostream<T> &stream, const GalaxyPtr &galaxy)
{
	stream << *galaxy;
	return stream;
}

}  // namespace shark

#endif /* INCLUDE_GALAXY_H_ */