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
 * Halo-related classes and functionality
 */

#ifndef INCLUDE_HALO_H_
#define INCLUDE_HALO_H_

#include <iosfwd>
#include <memory>
#include <vector>

#include "components.h"
#include "mixins.h"

namespace shark {

/**
 * A halo.
 *
 * Halos are the largest gravitationally bound structures in the universe. They
 * must contain at least one subhalo inside.
 */
class Halo : public Identifiable<halo_id_t>, public Spatial<float> {

public:

	Halo(id_t halo_id, int snapshot) :
		Identifiable(halo_id),
		snapshot(snapshot)
	{
		// no-op
	}

	/**
	 * @return the total number of subhalos contained in this halo
	 */
	std::size_t subhalo_count() const
	{
		std::size_t count = (central_subhalo ? 1 : 0);
		return count + satellite_subhalos.size();
	}

	/**
	 * Returns a new vector containing pointers to all subhalos contained in
	 * this halo (i.e., the central and satellite subhalos).
	 *
	 * @return A vector with all subhalos
	 */
	std::vector<SubhaloPtr> all_subhalos() const;

	/**
	 * Removes @a subhalo from this Halo. If the subhalo is not part of this
	 * Halo, a subhalo_not_found exception is thrown.
	 *
	 * @param subhalo The subhalo to remove
	 */
	void remove_subhalo(const SubhaloPtr &subhalo);

	/**
	 * Adds @a subhalo to this Halo.
	 *
	 * @param subhalo The subhalo to add
	 */
	void add_subhalo(SubhaloPtr &&subhalo);

	/**
	 * Returns the z=0 halo in which this halo ends up.
	 * @return
	 */
	HaloPtr final_halo() const;

	///
	/// Returns the number of galaxies contained in this Halo
	///
	galaxies_size_type galaxy_count() const;

	/**
	 * @return The total baryon mass contained in this Halo
	 */
	double total_baryon_mass() const;

	/**
	 * @return The baryon mass contained inside this Halo (exclusing ejected and lost mass)
	 */
	double inside_halo_baryon_mass() const;

	/**
	 * @return total mass of ascendants of this Halo
	 */
	double total_mass_ascendants() const;


	/**
	 * @return The main progenitor of this halo
	 */
	HaloPtr main_progenitor() const;

	/// The ascendant Halos of this Halo
	std::vector<HaloPtr> ascendants;
	/// The subhalos contained in this Halo
	std::vector<SubhaloPtr> satellite_subhalos;
	/// The central subhalo of this Halo
	SubhaloPtr central_subhalo;
	/// The descendant Halo of this Halo, if any
	HaloPtr descendant;
	/// The merger tree that holds this halo
	MergerTreePtr merger_tree;
	/**
	 * The mass contained in the subhalos.
	 * This quantity should be =1 for classic SAMs, but with Rodrigo Canas work
	 * on VELOCIraptor, this quantity could be less than 1.
	 */
	float mass_fraction_subhalos = -1;
	/// virial velocity of the halo [km/s]
	float Vvir = 0;
	/// dark matter mass of the halo [Msun/h]
	float Mvir = 0;
	/// gas mass in the halo [Msun/h]. This is different than 0 if the input simulation is a hydrodynamical simulation
	float Mgas = 0;
	/// NFW concentration parameter of halo
	float concentration = 0;
	/// spin parameter of halo
	float lambda = 0;
	/// redshift at which the halo had 80% of its mass in place
	float age_80 = 0;
	/// redshift at which the halo had 50% of its mass in place
	float age_50 = 0;
	/// cumulative excess jet power from satellite subhalos (this is >0 only if satellite completely shut off their own cooling and have some left over jet power)
	float excess_jetfeedback = 0;
	/// The snapshot at which this halo is found
	int snapshot;
	/**
	 * in the case the user runs the code with the option of ignoring late,
	 * massive forming halos (which are usually issues in the merger tree
	 * builder), this boolean parameter indicates whether this halo has been
	 * flagged as having this issue.
	 */
	bool ignore_gal_formation = false;

        /**
         * Boolean parameter tracking whether this halo is in quasi hydrostatic equilibrium
         */

        bool hydrostatic_eq = false;


};

template <typename T>
std::basic_ostream<T> &operator<<(std::basic_ostream<T> &stream, const Halo &halo)
{
	stream << "<Halo " << halo.id;
	if (halo.merger_tree) {
		stream << " @ " << halo.merger_tree;
	}
	stream << ">";
	return stream;
}

template <typename T>
std::basic_ostream<T> &operator<<(std::basic_ostream<T> &stream, const HaloPtr &halo)
{
	stream << *halo;
	return stream;
}

/**
 * Adds @p parent as the parent halo of @p halo, and likewise add @p halo as a
 * descendant of @p parent. If a descendant has already been set in @p parent
 * and is different from @p descendant then an error is thrown.
 *
 * @param halo The halo of which @p parent is a parent (i.e., a descendant of
 *        @p parent).
 * @param parent A parent of @p halo
 */
void add_parent(const HaloPtr &halo, const HaloPtr &parent);

/**
 * Utility class to be use as a unary predicate to filter HaloPtr objects belonging
 * to a given snapshot
 */
class in_snapshot {
public:
	in_snapshot() {}
	explicit in_snapshot(int snapshot)
	 : m_snapshot(snapshot)
	{}

	bool operator()(const HaloPtr &halo) const
	{
		return halo->snapshot == m_snapshot;
	}

private:
	int m_snapshot = -1;
};

/**
 * A binary predicate class that compares Halo objects according to their
 * snapshots.
 */
class by_snapshot {
public:
	bool operator()(const HaloPtr &halo, int snapshot) const
	{
		return halo->snapshot < snapshot;
	}

	bool operator()(int snapshot, const HaloPtr &halo) const
	{
		return snapshot < halo->snapshot;
	}

	bool operator()(const HaloPtr &a, const HaloPtr &b) const
	{
		return a->snapshot < b->snapshot;
	}
};

}  // namespace shark

#endif /* INCLUDE_HALO_H_ */
