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
 * Subhalo-related classes and functionality
 */

#ifndef INCLUDE_SUBHALO_H_
#define INCLUDE_SUBHALO_H_

#include <cstdint>
#include <iosfwd>
#include <memory>
#include <vector>

#include "baryon.h"
#include "mixins.h"

namespace shark {

// Forward definitions
class Halo;
class Galaxy;
class Subhalo;
using HaloPtr = std::shared_ptr<Halo>;
using GalaxyPtr = std::shared_ptr<Galaxy>;
using SubhaloPtr = std::shared_ptr<Subhalo>;
using galaxies_size_type = std::vector<GalaxyPtr>::size_type;

/**
 * This structure keeps track of the properties of the halo gas,
 * which are necessary to implement a more sophisticated cooling model.
 */
struct CoolingSubhaloTracking {
	std::vector<double> deltat;
	std::vector<double> temp;
	std::vector<double> mass;
	std::vector<double> tcooling;
	double rheat {0};
};


/** This class defines what a subhalo is. In an abstract sense, a subhalo is the quantum units of how dark matter clusters. Subhalos can
 * coexist with other subhalos in the same halo. A subhalo con also host 0, 1 or more galaxies depending on how we allow galaxies to populate
 * subhalos.
 */
class Subhalo : public Identifiable<std::int64_t>, public Spatial<float> {

public:

	/**
	 * Initialize values in zero.
	 */
	Subhalo(id_t id, int snapshot):
		Identifiable(id),
		snapshot(snapshot)
	{
		//no-op
	}

	/**
	 * An enumeration of types of subhalos
	 */
	enum subhalo_type_t {
		CENTRAL = 0,
		SATELLITE,
		FLYBY
	};

	/**
	 * The snapshot at which this subhalo is found
	 */
	int snapshot;

	/**
	 * Whether this subhalo has a descendant or not
	 */
	bool has_descendant = false;

	/**
	 * Boolean property indicating if subhalo is a main progenitor of its descendant.
	 */
	bool main_progenitor = false;

	/**
	 * Boolean property indicating if subhalo is the result of an interpolation in snapshots were descendants were missing. In this case Dhalos puts a subhalo in those snapshots
	 * to ensure continuation of the merger tree.
	 */
	bool IsInterpolated = false;

	/**
	 * The ID of the descendant of this subhalo.
	 * Valid only if has_descendant is \code{true}
	 */
	id_t descendant_id = 0;

	/**
	 * The ID of the Halo containing the descendant of this subhalo
	 */
	id_t descendant_halo_id = 0;

	/**
	 * The snapshot at which the descendant of this subhalo can be found
	 */
	int descendant_snapshot = -1;

    /**
     * Integer that shows if this subhalo will disappear from the tree in the next snapshot.
     * last_snapshot_identified = 1 if disappears in the next snapshot, =0 otherwise.
     */

    int last_snapshot_identified = -1;

	/**
	 * A pointer to the descendant of this subhalo.
	 * If this pointer is set then descendant_id and descendant_subhalo are
	 * meaningless.
	 */
	SubhaloPtr descendant;

	/**
	 * The list of galaxies in this subhalo.
	 */
	std::vector<GalaxyPtr> galaxies;

	/** Returns a pointer to the central galaxy. If no central galaxy is found
	 in this Subhalo, then an empty pointer is returned.
	 */
	GalaxyPtr central_galaxy() const;

	/** Returns a pointer to the type1 galaxy. If no type1 galaxy is found
	 in this Subhalo, then an empty pointer is returned.
	 */
	GalaxyPtr type1_galaxy() const;

	/**
	 * Returns all the type 2 satellites of this subhalo.
	 */
	std::vector<GalaxyPtr> all_type2_galaxies() const;

	/**
	 * The subhalo type
	 */
	subhalo_type_t subhalo_type = CENTRAL;

	/**
	 * haloID: ID of the Halo this Subhalo belong to
	 */
	id_t haloID = 0;

	/** Vvir: virial velocity of the subhalo [km/s]
	 * Mvir: dark matter mass of the subhalo [Msun/h]
	 * Mgas: gas mass in the subhalo [Msun/h]. This is different than 0 if the input simulation is a hydrodynamical simulation.
	 * L: angular momentum of subhalo [Msun/h km/s Mpc/h]
	 * Vcirc: maximum circular velocity of the subhalo [km/s]
	 * concentration: NFW concentration parameter of subhalo
	 * lambda: spin parameter of subhalo
	 * infall_t: redshift at which the subhalo became a type > 0.
	 *  */
	float Vvir = 0;
	float Mvir = 0;
	float Mgas = 0;
	xyz<float> L {0, 0, 0};
	float Vcirc = 0;
	float concentration = 0;
	float lambda = 0;
	float infall_t = 0;

	/**
	 * cooling_subhalo_tracking: saves que information of the virial temperature, total halo gas and cooling time history.
	 */
	CoolingSubhaloTracking cooling_subhalo_tracking;


	/**
	 * Hot gas component of the halo and outside the galaxies that is
	 * allowed to cool down and/or fall onto the galaxy.
	 */
	Baryon hot_halo_gas;

	/**
	 * Cold gas component of the halo and outside the galaxies that has
	 * cooled down.
	 */
	Baryon cold_halo_gas;

	/**
	 * Hot gas component of the halo and outside galaxies that tracks
	 * the ejected outflowing gas from the galaxy and that is not
	 * available for cooling yet.
	 */
	Baryon ejected_galaxy_gas;

	/**
	 * Lost gas reservoir which tracks the gas that is outflowing due to
	 * QSO feedback and that has escaped the halo.
	 */
	Baryon lost_galaxy_gas {};

	/**
	 * A list of pointers to the ascendants of this subhalo, sorted by mass in
	 * descending order
	 */
	std::vector<SubhaloPtr> ascendants;

	/**
	 * The accreted baryonic mass onto the subhalo. This information comes from the merger tree.
	 */
	float accreted_mass = 0;

	/**
	 * The halo that holds this subhalo.
	 */
	HaloPtr host_halo;

	/**
	 * @return The main progenitor of this Subhalo
	 */
	SubhaloPtr main() const;

	/**
	 * Copies the galaxies from this Subhalo into @a target
	 *
	 * @param target The subhalo where galaxies will be copied to
	 * @param gals The galaxies to copy to the target subhalo; defaults to all our galaxies
	 */
	void copy_galaxies_to(SubhaloPtr &target, const std::vector<GalaxyPtr> &gals) const;

	/**
	 * Transfers (i.e., moves) the galaxies from this Subhalo into @a target
	 *
	 * @param target The subhalo where galaxies will be transferred to
	 */
	void transfer_galaxies_to(SubhaloPtr &target);

	/**
	 * Transfers (i.e., moves) the galaxies that are type=2 from this Subhalo into @a target
	 *
	 * @param target The subhalo where type 2 galaxies will be transferred to
	 */
	void transfer_type2galaxies_to(SubhaloPtr &target);

	/**
	 * Removes galaxies from this Subhalo
	 *
	 * @param to_remove A vector of galaxies to remove.
	 */
	void remove_galaxies(const std::vector<GalaxyPtr> &to_remove);

	/// Returns the number of galaxies contained in this Halo
	galaxies_size_type galaxy_count() const
	{
		return galaxies.size();
	}

	/**
	 * @return The total baryon mass contained in this Subhalo
	 */
	double total_baryon_mass() const;

	/**
	 * Checks whether this subhalo has the correct composition of galaxies
	 * depending on its type
	 */
	void check_subhalo_galaxy_composition() const;

	/**
	 * Like check_subhalo_galaxy_composition, but fails if this subhalo is not
	 * a central subhalo
	 */
	void check_central_subhalo_galaxy_composition() const;

	/**
	 * Like check_subhalo_galaxy_composition, but fails if this subhalo is not
	 * a satellite subhalo
	 */
	void check_satellite_subhalo_galaxy_composition() const;

private:
	void do_check_satellite_subhalo_galaxy_composition() const;
	void do_check_central_subhalo_galaxy_composition() const;

};

template <typename T>
std::basic_ostream<T> &operator<<(std::basic_ostream<T> &stream, const Subhalo &subhalo)
{
	stream << "<Subhalo " << subhalo.id;
	if (subhalo.host_halo) {
		stream << " @ " << subhalo.host_halo;
	}
	stream << ">";
	return stream;
}

template <typename T>
std::basic_ostream<T> &operator<<(std::basic_ostream<T> &stream, const SubhaloPtr &subhalo)
{
	stream << *subhalo;
	return stream;
}

}  // namespace shark

#endif /* INCLUDE_SUBHALO_H_ */