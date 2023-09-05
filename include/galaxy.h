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
#include "numerical_constants.h"

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

struct BHHistoryItem {
	float macc_hh;
	float macc_sb;
	float massembly;
	float mbh;
	float spin;
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

class Galaxy : public Identifiable<galaxy_id_t> {
public:

	using Identifiable::Identifiable;

	/// Deleted copy constructor to avoid copying Galaxy objects
	Galaxy(const Galaxy &other) = delete;
	/// Move constructor
	Galaxy(Galaxy &&other) = default;

	/// Deleted copy assignment operator to avoid copying Galaxy objects
	Galaxy &operator=(const Galaxy &other) = delete;
	/// Move assignment operator
	Galaxy &operator=(Galaxy &&other) = default;

	/**
	 * An enumeration of types of galaxies
	 */
	enum galaxy_type_t {
		CENTRAL = 0,
		TYPE1,
		TYPE2,
		FLYBY
	};

	Baryon bulge_stars;
	Baryon bulge_gas;
	Baryon disk_stars;
	Baryon disk_gas;
	BaryonBase galaxymergers_burst_stars;
	BaryonBase galaxymergers_assembly_stars;
	BaryonBase diskinstabilities_burst_stars;
	BaryonBase diskinstabilities_assembly_stars;
	BlackHole smbh;

	// save first snapshot in which galaxy appears.
	int birth_snapshot;

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

	// ram-pressure stripping radius
	float r_rps = 0;
	BaryonBase ram_pressure_stripped_gas;

	/// star formation and gas history of this galaxy across snapshots
	std::vector<HistoryItem> history;

	/// black hole history of this galaxy across snapshots
	std::vector<BHHistoryItem> bh_history;

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
	/// tracking of mass lost due to tidal stripping
	BaryonBase stars_tidal_stripped;
	/// The ID of the descendant of this galaxy, -1 if no descendant is defined
	id_t descendant_id = -1;
	/// The type of this galaxy
	galaxy_type_t galaxy_type = CENTRAL;


        float mheat_ratio = 0;

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

	double enclosed_mass_disk(float r) const
	{

		return enclosed_mass_exponential(r, disk_gas.mass, disk_gas.rscale) +
				enclosed_mass_exponential(r, disk_stars.mass, disk_stars.rscale);
	}

	double enclosed_bulge_mass(float r) const
	{

		return enclosed_mass_exponential(r, bulge_gas.mass, bulge_gas.rscale) +
				enclosed_mass_plummer(r, bulge_stars.mass, bulge_stars.rscale);

	}

	double enclosed_mass_gas(float r) const
	{

		return enclosed_mass_exponential(r, disk_gas.mass, disk_gas.rscale) +
				enclosed_mass_exponential(r, bulge_gas.mass, bulge_gas.rscale);
	}

	double enclosed_mass_exponential(float r, float m, float r50) const
	{
		if(m == 0){
			return 0;
		}
		else{
			auto rnorm = r/(r50/1.67);
			return m * (1 - (1 + rnorm) * std::exp(-rnorm));
		}
	}

	double enclosed_mass_plummer(float r, float m, float r50) const
	{
		if(m ==0){
			return 0;
		}
		else{
			auto re = r50/1.3;
			return m * std::pow(r, 3) / std::pow( std::pow(r,2) + std::pow(re,2), 1.5);
		}
	}

	double surface_density_disk(float r){
		return surface_density_exponential(r, disk_gas.mass, disk_gas.rscale) +
				surface_density_exponential(r, disk_stars.mass, disk_stars.rscale);
	}

	double surface_density_bulge(float r){
		return surface_density_exponential(r, bulge_gas.mass, bulge_gas.rscale) +
				surface_density_plummer(r, bulge_stars.mass, bulge_stars.rscale);
	}

	double surface_density_gas(float r){
		return surface_density_exponential(r, disk_gas.mass, disk_gas.rscale) +
				surface_density_exponential(r, bulge_gas.mass, bulge_gas.rscale);
	}

	double surface_density_stars(float r){
		return surface_density_exponential(r, disk_stars.mass, disk_stars.rscale) +
				surface_density_plummer(r, bulge_stars.mass, bulge_stars.rscale);
	}

	double surface_density_exponential(float r, float m, float r50) const
	{
		if(m == 0){
			return 0;
		}
		else{
			auto re = r50/1.67;
			auto rnorm = r/re;
			return m / (constants::PI2 * std::pow(re,2)) * std::exp(-rnorm);
		}
	}

	double surface_density_plummer(float r, float m, float r50) const
	{
		if(m ==0){
			return 0;
		}
		else{
			auto re = r50/1.3;
			return m * std::pow(re, 2) / (constants::PI * std::pow( std::pow(r,2) + std::pow(re,2), 2));
		}
	}

};

// Support for less-based comparison of Galaxy objects. This allows them to be
// put into a set or be keys for a std::map
inline
bool operator<(const Galaxy &a, const Galaxy &b)
{
	return a.id < b.id;
}

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

template <typename T>
std::basic_ostream<T> &operator<<(std::basic_ostream<T> &stream, const ConstGalaxyPtr &galaxy)
{
	stream << *galaxy;
	return stream;
}

}  // namespace shark

#endif /* INCLUDE_GALAXY_H_ */
