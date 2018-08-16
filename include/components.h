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
 *
 * Base classes used create an evolving system.
 */

#ifndef SHARK_COMPONENTS_H_
#define SHARK_COMPONENTS_H_

#include <cassert>
#include <iterator>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <vector>

#include "mixins.h"

namespace shark {

// Forward-defines
class Galaxy;
class Subhalo;
class Halo;
class MergerTree;

typedef std::shared_ptr<Galaxy> GalaxyPtr;
typedef std::shared_ptr<Subhalo> SubhaloPtr;
typedef std::shared_ptr<Halo> HaloPtr;
typedef std::shared_ptr<MergerTree> MergerTreePtr;

/// Type used by galaxy_count() methods
typedef typename std::vector<GalaxyPtr>::size_type galaxies_size_type;

/**
 * The common base for all baryon component types.
 */
class BaryonBase {

public:

	/**
	 * Mass content of the baryon component
	 */
	float mass = 0;

	/**
	 * Metallicity of the baryon component
	 */
	float mass_metals = 0;

	BaryonBase &operator+=(const BaryonBase &b) {
		mass += b.mass;
		mass_metals += b.mass_metals;
		return *this;
	}

	friend BaryonBase operator+(BaryonBase &lhs, const BaryonBase &rhs) {
		lhs += rhs;
		return lhs;
	}

	void restore_baryon(){
		mass = 0;
		mass_metals = 0;
	}

};

/**
 * A common baryon component.
 * Note that black holes are not baryon components as they use their own class.
 */
class Baryon : public BaryonBase {
public:

	/**
	 * A scale radius
	 */
	float rscale = 0;

	/**
	 * Specific angular momentum
	 */
	float sAM = 0;

	friend Baryon operator+(Baryon &lhs, const Baryon &rhs) {
		lhs += rhs;
		return lhs;
	}

	void restore_baryon(){
		BaryonBase::restore_baryon();
		rscale = 0;
		sAM = 0;
	}

	float angular_momentum(){
		return mass * sAM;
	}

};

/**
 * Black hole baryon component.
 *
 * Because of the no hair theorem, black holes are only allowed to have a mass,
 * an accretion rate and a spin. This implies that to extend black holes from
 * the basic galaxy_component class, one needs to include at most an accretion
 * rate and a spin.
 */
class BlackHole : public BaryonBase {

public:

	/** accretion rate onto the black hole during hot halo mode. */
	float macc_hh = 0;

	/** macc_sb: accretion rate onto the black hole during starbursts. */
	float macc_sb = 0;
};



/**
 * Structure that saves the history of relevant baryon components needed for SED calculation later on.
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

struct InteractionItem{
	int major_mergers = 0;
	int minor_mergers = 0;
	int disk_instabilities = 0;

        void restore_interaction_item(){
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

	Baryon bulge_stars {};
	Baryon bulge_gas {};
	Baryon disk_stars {};
	Baryon disk_gas {};
	Baryon galaxymergers_burst_stars {};
	Baryon galaxymergers_assembly_stars {};
	Baryon diskinstabilities_burst_stars {};
	Baryon diskinstabilities_assembly_stars {};
	BlackHole smbh {};

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

	//save maximum circular velocity.
	float vmax = 0;

	//save star formation and gas history
	std::vector<HistoryItem>  history {};

	//save interactions of this galaxy during this snapshot.
	InteractionItem interaction {};

	/**
	 * tmerge: dynamical friction timescale, which is defined only if galaxy is satellite.
	 * concentration_type2: concentration of the subhalo this galaxy was before becoming type 2 (only relevant for type 2 galaxies).
	 * msubhalo_type2: subhalo mass of this galaxy before it became type 2 (only relevant for type 2 galaxies).
	 * vvir_type2: subhalo virial velocity of this galaxy before it became type 2 (only relevant for type 2 galaxies).
	 * lambda_type2: subhalo spin parameter of this galaxy before it became type 2 (only relevant for type 2 galaxies).
	 */
	float tmerge = 0;
	float concentration_type2 = 0;
	float msubhalo_type2 = 0;
	float vvir_type2 = 0;
	float lambda_type2 = 0;

	/**
	 * Define functions to calculate total mass and metals of various components.
	 */

	double disk_mass(){
		return disk_gas.mass + disk_stars.mass;
	}

	double disk_mass_metals(){
		return disk_gas.mass_metals + disk_stars.mass_metals;
	}

	double bulge_mass(){
		return bulge_gas.mass + bulge_stars.mass;
	}

	double bulge_mass_metals(){
		return bulge_gas.mass_metals + bulge_stars.mass_metals;
	}

	double baryon_mass(){
		return disk_gas.mass + disk_stars.mass + bulge_gas.mass + bulge_stars.mass;
	}

	double stellar_mass(){
		return disk_stars.mass + bulge_stars.mass;
	}

	double stellar_mass_metals(){
		return disk_stars.mass_metals + bulge_stars.mass_metals;
	}

	double gas_mass(){
		return disk_gas.mass + bulge_gas.mass;
	}

	double gas_mass_metals(){
		return disk_gas.mass_metals + bulge_gas.mass_metals;
	}

	double disk_size(){

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

	double bulge_size(){


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

	double composite_size(){

		double rdisk = disk_size();
		double rbulge = bulge_stars.rscale;

		double rcomp = 0.0;

		// Define rcomp only if galaxy has mass.
		if(baryon_mass() > 0){
			rcomp = (disk_mass() * rdisk + bulge_mass() * rbulge) / baryon_mass();
		}

		return rcomp;
	}

	double stellar_size(){

		double rdisk = disk_size();
		double rbulge = bulge_size();

		double rcomp = 0.0;

		// Define rcomp only if galaxy has mass.
		if(baryon_mass() > 0){
			rcomp = (disk_stars.mass * rdisk + bulge_mass() * rbulge) / baryon_mass();
		}

		return rcomp;
	}

	double angular_momentum(){

		return disk_gas.angular_momentum() + disk_stars.angular_momentum() + bulge_gas.angular_momentum() + bulge_stars.angular_momentum();
	}

};


/**
 * This structure keeps track of the properties of the halo gas, which are necessary to implement a more sophisticated cooling model.
 */
struct CoolingSubhaloTracking {
	/**
	 * Initialize values in zero.
	 */
	CoolingSubhaloTracking():
		deltat(),
		temp(),
		mass(),
		tcooling(),
		rheat(0)
	{
		//no=op
	};
	std::vector<double> deltat;
	std::vector<double> temp;
	std::vector<double> mass;
	std::vector<double> tcooling;
	double rheat;
};


/** This class defines what a subhalo is. In an abstract sense, a subhalo is the quantum units of how dark matter clusters. Subhalos can
 * coexist with other subhalos in the same halo. A subhalo con also host 0, 1 or more galaxies depending on how we allow galaxies to populate
 * subhalos.
 */
class Subhalo : public Identifiable<long>, public Spatial<float> {

public:

	/**
	 * Initialize values in zero.
	 */
	Subhalo(long id, int snapshot):
		Identifiable(id),
		Spatial(),
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
	SubhaloPtr descendant {};

	/**
	 * The list of galaxies in this subhalo.
	 */
	std::vector<GalaxyPtr> galaxies {};

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
	 * Mvir: virial mass of the subhalo [Msun/h]
	 * L: angular momentum of subhalo [Msun/h km/s Mpc/h]
	 * Vcirc: maximum circular velocity of the subhalo [km/s]
	 * concentration: NFW concentration parameter of subhalo
	 * lambda: spin parameter of subhalo
	 *  */
	float Vvir = 0;
	float Mvir = 0;
	xyz<float> L {0, 0, 0};
	float Vcirc = 0;
	float concentration = 0;
	float lambda = 0;

	/**
	 * cooling_subhalo_tracking: saves que information of the virial temperature, total halo gas and cooling time history.
	 */
	CoolingSubhaloTracking cooling_subhalo_tracking {};


	/**
	 * Hot gas component of the halo and outside the galaxies that is
	 * allowed to cool down and/or fall onto the galaxy.
	 */
	Baryon hot_halo_gas {};

	/**
	 * Cold gas component of the halo and outside the galaxies that has
	 * cooled down.
	 */
	Baryon cold_halo_gas {};

	/**
	 * Hot gas component of the halo and outside galaxies that tracks
	 * the ejected outflowing gas from the galaxy and that is not
	 * available for cooling yet.
	 */
	Baryon ejected_galaxy_gas {};

	/**
	 * A list of pointers to the ascendants of this subhalo, sorted by mass in
	 * descending order
	 */
	std::vector<SubhaloPtr> ascendants {};

	/**
	 * The accreted baryonic mass onto the subhalo. This information comes from the merger tree.
	 */
	float accreted_mass = 0;

	/**
	 * The halo that holds this subhalo.
	 */
	HaloPtr host_halo {};

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

namespace detail {

template <bool constant>
struct subhalos_iterator_traits {
};

template <>
struct subhalos_iterator_traits<false> {
	typedef Halo * halo_pointer;
	typedef Halo & halo_reference;
	typedef SubhaloPtr & subhalo_reference;
	typedef SubhaloPtr * subhalo_pointer;
	typedef std::vector<SubhaloPtr>::iterator satellite_iterator;
};

template <>
struct subhalos_iterator_traits<true> {
	typedef const Halo * halo_pointer;
	typedef const Halo & halo_reference;
	typedef const SubhaloPtr & subhalo_reference;
	typedef const SubhaloPtr * subhalo_pointer;
	typedef std::vector<SubhaloPtr>::const_iterator satellite_iterator;
};

}  // namespace detail

/**
 * A halo.
 *
 * Halos are the largest gravitationally bound structures in the universe. They
 * must contain at least one subhalo inside.
 */
class Halo : public Identifiable<long>, public Spatial<float> {

public:

	template <bool constant>
	class subhalos_iterator {

	private:
		bool central;
		typename detail::subhalos_iterator_traits<constant>::halo_pointer halo;
		typename detail::subhalos_iterator_traits<constant>::satellite_iterator satellite_it {};

	public:
		typedef ptrdiff_t difference_type;
		typedef SubhaloPtr value_type;
		typedef std::forward_iterator_tag iterator_category;
		typedef typename detail::subhalos_iterator_traits<constant>::subhalo_pointer pointer;
		typedef typename detail::subhalos_iterator_traits<constant>::subhalo_reference reference;
		typedef typename detail::subhalos_iterator_traits<constant>::halo_reference halo_reference;

		subhalos_iterator(halo_reference halo) : central(false), halo(&halo)
		{
			if (!halo.central_subhalo && halo.satellite_subhalos.empty()) {
				this->halo = nullptr;
			}
			else if (halo.central_subhalo) {
				central = true;
			}
			else {
				satellite_it = halo.satellite_subhalos.begin();
			}
		}
		subhalos_iterator() : central(false), halo(nullptr) {};

		reference operator*() const
		{
			if (central)
				return halo->central_subhalo;
			return *satellite_it;
		}

		subhalos_iterator &operator++()
		{
			if (central) {
				central = false;
				satellite_it = halo->satellite_subhalos.begin();
			}
			else {
				satellite_it++;
			}
			if (satellite_it == halo->satellite_subhalos.end()) {
				halo = nullptr;
			}
			return *this;
		}

		bool operator ==(const subhalos_iterator &lhs) const
		{
			if (bool(halo) == bool(lhs.halo))
				return true;
			if (central && lhs.central)
				return true;
			return central == lhs.central && satellite_it == lhs.satellite_it;
		}

		bool operator !=(const subhalos_iterator &lhs) const
		{
			return !(*this == lhs);
		}
	};

	template <typename HaloT>
	class subhalos_view {

	public:
		typedef subhalos_iterator<false> iterator;
		typedef subhalos_iterator<true> const_iterator;

		subhalos_view(HaloT &halo) : halo(halo) {};

		const_iterator begin() const
		{
			return const_iterator(halo);
		}

		iterator begin()
		{
			return iterator(halo);
		}

		const_iterator end() const
		{
			return const_iterator();
		}

		iterator end()
		{
			return iterator();
		}

	private:
		HaloT &halo;
	};

	Halo(id_t halo_id, int snapshot) :
		Identifiable(halo_id),
		Spatial(),
		snapshot(snapshot)
	{
		// no-op
	}

	/**
	 * The central subhalo
	 */
	SubhaloPtr central_subhalo {};

	/**
	 * The subhalos contained in this halo
	 */
	std::vector<SubhaloPtr> satellite_subhalos {};

	/**
	 * @return the total number of subhalos contained in this halo
	 */
	unsigned long subhalo_count() const
	{
		unsigned int count = (central_subhalo ? 1 : 0);
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
	 * The mass contained in the subhalos.
	 * This quantity should be =1 for classic SAMs, but with Rodrigo Canas work
	 * on VELOCIraptor, this quantity could be less than 1.
	 */
	float mass_fraction_subhalos = -1;

	/** TODO: document these */
	float Vvir = 0;
	float Mvir = 0;
	float concentration = 0;
	float lambda = 0;
	float cooling_rate = 0;

	/**
	 * The snapshot at which this halo is found
	 */
	int snapshot;

	bool main_progenitor = false;

	HaloPtr descendant {};
	std::set<HaloPtr> ascendants {};

	/**
	 * The merger tree that holds this halo.
	 */
	MergerTreePtr merger_tree {};

	/**
	 * Adds @a subhalo to this Halo.
	 *
	 * @param subhalo The subhalo to add
	 */
	void add_subhalo(const SubhaloPtr &&subhalo);

	///
	/// Returns the number of galaxies contained in this Halo
	///
	galaxies_size_type galaxy_count() const;

	/**
	 * @return The total baryon mass contained in this Halo
	 */
	double total_baryon_mass() const;

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
 * A merger tree.
 *
 * A merger tree contains halos, which are indexed by snapshot,
 * and an ID to identify it.
 */
class MergerTree : public Identifiable<long> {
public:

	using Identifiable::Identifiable;

	/**
	 * All halos contained in this merger tree, indexed by snapshot number
	 */
	std::map<int, std::vector<HaloPtr>> halos;

	void add_halo(const HaloPtr &halo) {
		halos[halo->snapshot].push_back(halo);
	}

	std::vector<HaloPtr> &halos_at(int snapshot)
	{
		static std::vector<HaloPtr> empty;
		if (halos.find(snapshot) == halos.end()) {
			return empty;
		}
		return halos.at(snapshot);
	}
};

class TotalBaryon {

public:

	/**
	 * Class to track all the mass componets in the simulation.
	 */

	TotalBaryon() :
		SFR_disk(0),
		SFR_bulge(0)
	{
		// no-op
	}

	/**
	 * mcold: total cold gas mass in disk+bulge.
	 * mstars: total stellar mass in disk+bulge.
	 * mstars_burst_galaxymergers: total stellar mass formed via bursts driven by galaxy mergers.
	 * mstars_burst_diskinstabilities: total stellar mass formed via bursts driven by disk instabilities.
	 * mhot_halo: total hot halo gas.
	 * mcold_halo: total cold halo gas (that is cooling during the current snapshot).
	 * mejected_halo: total hot gas that has been ejected from galaxies due to feeback and that has not been reincorporated yet onto the hot halo gas reservoir.
	 * mBH: total mass locked in black holes.
	 * mHI: total mass in the form of atomic gas.
	 * mH2: total mass in the form of molecular gas.
	 * mDM: total mass in the form of dark matter.
	 * SFR: integrated SFR of all galaxies over a snapshot.
	 * baryon_total_created: keeps track of the baryons deposited in DM halos to ensure mass convervations.
	 */

	std::vector<BaryonBase> mcold;
	std::vector<BaryonBase> mstars;
	std::vector<BaryonBase> mstars_burst_galaxymergers;
	std::vector<BaryonBase> mstars_burst_diskinstabilities;
	std::vector<BaryonBase> mhot_halo;
	std::vector<BaryonBase> mcold_halo;
	std::vector<BaryonBase> mejected_halo;
	std::vector<BaryonBase> mBH;
	std::vector<BaryonBase> mHI;
	std::vector<BaryonBase> mH2;
	std::vector<BaryonBase> mDM;

	std::vector<double> SFR_disk;
	std::vector<double> SFR_bulge;

	/**
	 * Vectors of integers that keep track of number of mergers of galaxies and disk instability episodes in each snapshot.
	 * major_mergers: counts number of major mergers per snapshot.
	 * minor_mergers: counts number of minor mergers per snapshot.
	 * disk_instabil: counts number of disk instability episodes per snapshot.
 	*/ 
	std::vector<int> major_mergers;
	std::vector<int> minor_mergers;
	std::vector<int> disk_instabil;

	std::map<int,double> baryon_total_created;
	std::map<int,double> baryon_total_lost;

	std::vector<double> get_masses (const std::vector<BaryonBase> &B) const;
	std::vector<double> get_metals (const std::vector<BaryonBase> &B) const;

};

template <typename T>
std::basic_ostream<T> &operator<<(std::basic_ostream<T> &stream, const MergerTree &merger_tree)
{
	stream << "<MergerTree " << merger_tree.id << ">";
	return stream;
}

template <typename T>
std::basic_ostream<T> &operator<<(std::basic_ostream<T> &stream, const MergerTreePtr &merger_tree)
{
	stream << *merger_tree;
	return stream;
}

}  // namespace shark

#endif // SHARK_COMPONENTS_H_
