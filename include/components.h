//
// Base classes that make up a SHArk solving system.
//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2017
// Copyright by UWA (in the framework of the ICRAR)
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307  USA
//

#ifndef SHARK_COMPONENTS_H_
#define SHARK_COMPONENTS_H_

#include <algorithm>
#include <cassert>
#include <map>
#include <memory>
#include <numeric>
#include <ostream>
#include <vector>

#include "logging.h"
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


/**
 * The common base for all baryon component types.
 */
class BaryonBase {

public:

	/**
	 * Initialize values in zero.
	 */
	BaryonBase():
		mass(0),
		mass_metals(0)
	{
		// no-op
	}
	/**
	 * Mass content of the baryon component
	 */
	float mass;

	/**
	 * Metallicity of the baryon component
	 */
	float mass_metals;
};

/**
 * A common baryon component.
 * Note that black holes are not baryon components as they use their own class.
 */
class Baryon : public BaryonBase {
public:
	/**
	 * Initialize values in zero.
	 */
	Baryon():
		rscale(0),
		sAM(0)
	{
		// no-op
	}
	/**
	 * A scale radius
	 */
	float rscale;

	/**
	 * Specific angular momentum
	 */
	float sAM;
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

	/**
	 * Initialize values in zero.
	 */
	BlackHole():
		macc_hh(0),
		macc_sb(0)
	{
		// no-op
	}
	/**
	 * macc_hh: accretion rate onto the black hole during hot halo mode.
	 * macc_sb: accretion rate onto the black hole during starbursts.
	 */

	float macc_hh;
	float macc_sb;
};



/**
 * Structure that saves the history of relevant baryon components needed for SED calculation later on.
 */
struct HistoryItem {

	float sfr_disk;
	float sfr_bulge;
	Baryon stellar_disk;
	Baryon stellar_bulge;
	Baryon gas_disk;
	Baryon gas_bulge;
	int snapshot;

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

	/**
	 * Initialize values in zero.
	 */
	Galaxy():
		galaxy_type(),
		bulge_stars(),
		bulge_gas(),
		disk_stars(),
		disk_gas(),
		smbh(),
		sfr_disk(0),
		sfr_bulge(0),
		tmerge(0)
	{
		//no-op
	}

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
	 * The type of galaxy
	 */
	galaxy_type_t galaxy_type;

	Baryon bulge_stars;
	Baryon bulge_gas;
	Baryon disk_stars;
	Baryon disk_gas;
	Baryon burst_stars;
	BlackHole smbh;

	//save average star formation rates.
	float sfr_disk;
	float sfr_bulge;

	//save star formation and gas history
	std::vector<HistoryItem>  history;

	/**
	 * dynamical friction timescale, which is defined only is galaxy is satellite.
	 */
	float tmerge;

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


	double composite_size(){


		double rdisk = 0;
		double rbulge = 0;

		if(disk_mass() > 0){
			rdisk = (disk_stars.mass * disk_stars.rscale + disk_gas.mass * disk_gas.rscale) / disk_mass();
		}
		if(bulge_mass() > 0){
			rbulge = (bulge_stars.mass * bulge_stars.rscale + bulge_gas.mass * bulge_gas.rscale) / bulge_mass();
		}

		double rcomp = 0.0;

		// Define rcomp only if galaxy has mass.
		if(baryon_mass() > 0){
			rcomp = (disk_mass() * rdisk + bulge_mass() * rbulge) / baryon_mass();
		}

		return rcomp;
	}
};

/** This class extends the galaxy to include spatial information.*/
class SpatialGalaxy : public Galaxy, public Spatial<float> {
};

/** Extend galaxy to be satellite by including a merging timescale. */
class SatelliteGalaxy : public Galaxy {
public:
	float tmerge;
};

/** This class extends the satellite galaxy to include spatial information.*/
class SpatialSatelliteGalaxy : public SatelliteGalaxy, public Spatial<float> {
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
	Subhalo():
		snapshot(-1),
		has_descendant(false),
		main_progenitor(false),
		IsInterpolated(false),
		descendant_id(0),
		descendant_halo_id(0),
		descendant_snapshot(-1),
		last_snapshot_identified(-1),
		descendant(0),
		galaxies(),
		subhalo_type(),
		haloID(0),
		Vvir(0),
		Mvir(0),
		L{0, 0, 0},
		Vcirc(0),
		concentration(0),
		lambda(0),
		cooling_subhalo_tracking(),
		hot_halo_gas(),
		cold_halo_gas(),
		ejected_galaxy_gas(),
		ascendants(),
		accreted_mass(0),
		host_halo()
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
	bool has_descendant;

	/**
	 * Boolean property indicating if subhalo is a main progenitor of its descendant.
	 */
	bool main_progenitor;

	/**
	 * Boolean property indicating if subhalo is the result of an interpolation in snapshots were descendants were missing. In this case Dhalos puts a subhalo in those snapshots
	 * to ensure continuation of the merger tree.
	 */
	bool IsInterpolated;


	/**
	 * The ID of the descendant of this subhalo.
	 * Valid only if has_descendant is \code{true}
	 */
	id_t descendant_id;

	/**
	 * The ID of the Halo containing the descendant of this subhalo
	 */
	id_t descendant_halo_id;

	/**
	 * The snapshot at which the descendant of this subhalo can be found
	 */
	int descendant_snapshot;

    /**
     * Integer that shows if this subhalo will disappear from the tree in the next snapshot.
     * last_snapshot_identified = 1 if disappears in the next snapshot, =0 otherwise.
     */

    int last_snapshot_identified;

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

	/**
	 * The subhalo type
	 */
	subhalo_type_t subhalo_type;

	/**
	 * Which Halo does this Subhalo belong to
	 */
	id_t haloID;

	/** TODO: Properly document these */
	float Vvir;
	float Mvir;
	xyz<float> L;
	float Vcirc;
	float concentration;
	float lambda;

	/**
	 * This component saves que information of the virial temperature, total halo gas and cooling time history.
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
	 * A list of pointers to the ascendants of this subhalo, sorted by mass in
	 * descending order
	 */
	std::vector<SubhaloPtr> ascendants;

	/**
	 * The accreted baryonic mass onto the subhalo. This information comes from the merger tree.
	 */
	float accreted_mass;

	/**
	 * The halo that holds this subhalo.
	 */
	HaloPtr host_halo;

	// Sort ascendant subhalos by Mvir
	std::vector<SubhaloPtr> ordered_ascendants(){

		if(ascendants.size()==0){
			return std::vector<SubhaloPtr>();
		}
		else if(ascendants.size()>1){
			std::sort(ascendants.begin(), ascendants.end(), [](const SubhaloPtr &lhs, const SubhaloPtr &rhs) {
			return lhs->Mvir > rhs->Mvir;
			});
		}

		return ascendants;

	}

	/// Returns main progenitor subhalo.
	SubhaloPtr main(){
		for (auto &sub: ascendants) {
			if (sub->main_progenitor) {
				return sub;
			}
		}
		return SubhaloPtr();

	}

	/// Returns a pointer to the central galaxy. If no central galaxy is found
	/// in this Subhalo, then an empty pointer is returned.
	GalaxyPtr central_galaxy(){
		for (auto galaxy: galaxies){
			if(galaxy->galaxy_type == Galaxy::CENTRAL){
				return galaxy;
			}
		}
		return GalaxyPtr();
	}

	/// Copies the galaxies from this Subhalo into `target`
	void copy_galaxies_to(SubhaloPtr &target) {
		target->galaxies.insert(target->galaxies.end(), galaxies.begin(), galaxies.end());
	}

	/// Transfers (i.e., moves) the galaxies from this Subhalo into `target`
	void transfer_galaxies_to(SubhaloPtr &target) {

		auto gals_before = target->galaxy_count();
		auto our_gals = galaxies.size();
		LOG(trace) << "Transferring " << our_gals << " galaxies from " << *this << " to " << target << " (currently " << gals_before << " galaxies)";

		copy_galaxies_to(target);
		galaxies.clear();

		auto gals_after = target->galaxy_count();
		assert(gals_before + our_gals == gals_after);
	}

	void remove_galaxies(const std::vector<GalaxyPtr> &to_remove) {
		// TODO: Maybe not most efficiently, but it will do for now
		for(auto &galaxy: to_remove) {
			auto it = std::find(galaxies.begin(), galaxies.end(), galaxy);
			if (it == galaxies.end()) {
				LOG(warning) << "Trying to remove galaxy " << galaxy << " which is not in subhalo " << *this << ", ignoring";
				continue;
			}
			LOG(debug) << "Removing galaxy " << galaxy << " from subhalo " << *this;
			galaxies.erase(it);
		}
	}
	///
	/// Returns the number of galaxies contained in this Halo
	///
	unsigned long galaxy_count() {
		return galaxies.size();
	}

	// Sort galaxies by baryon mass.
	std::vector<GalaxyPtr> ordered_galaxies(){

		if(galaxies.size()==0){
			return std::vector<GalaxyPtr>();
		}
		else if(galaxies.size()>1){
			std::sort(galaxies.begin(), galaxies.end(), [](const GalaxyPtr &lhs, const GalaxyPtr &rhs) {
			return lhs->baryon_mass() > rhs->baryon_mass();
			});
		}

		return galaxies;

	}

	double total_baryon_mass(){

		double mass= 0.0;

		// add halo components.
		mass += hot_halo_gas.mass + cold_halo_gas.mass + ejected_galaxy_gas.mass;

		for (auto &galaxy: galaxies){
			mass += galaxy->baryon_mass() + galaxy->smbh.mass;
		}

		return mass;
	}

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

/**
 * Class to extend the properties a subhalo can have, by allowing it to have more baryon components than the basic subhalo.
 */
class SuperSubhalo : public Subhalo {
public:
	float rscale;
	float sAM;
	Baryon ejected_galaxy_gas; /*Gas that has been ejected by the galaxy but lives inside the halo.*/
	Baryon ejected_halo_gas; /*Gas that has been ejected outside the halo.*/
	Baryon halo_stars; /*Stars that live in the halo and outside the galaxy.*/
};


/**
 * A halo.
 *
 * Halos are the largest gravitationally bound structures in the universe. They
 * must contain at least one subhalo inside.
 */
class Halo : public Identifiable<long>, public Spatial<float> {

public:

	Halo(Halo::id_t halo_id, int snapshot) :
		central_subhalo(),
		satellite_subhalos(),
		mass_fraction_subhalos(-1),
		Vvir(0),
		Mvir(0),
		concentration(0),
		lambda(0),
		cooling_rate(0),
		snapshot(snapshot),
		main_progenitor(false)
	{
		// no-op
		id = halo_id;
	}

	/**
	 * The central subhalo
	 */
	SubhaloPtr central_subhalo;

	/**
	 * The subhalos contained in this halo
	 */
	std::vector<SubhaloPtr> satellite_subhalos;

	///
	/// Returns the total number of subhalos contained in this halo
	///
	unsigned long subhalo_count() {
		unsigned int count = (central_subhalo ? 1 : 0);
		return count + satellite_subhalos.size();
	}

	///
	/// Returns a new vector containing pointers to all subhalos contained in
	/// this halo (i.e., the central and satellite subhalos).
	///
	/// @return A vector with all subhalos
	///
	std::vector<SubhaloPtr> all_subhalos() {

		std::vector<SubhaloPtr> all;

		if (central_subhalo) {
			all.push_back(central_subhalo);
		}
		all.insert(all.end(), satellite_subhalos.begin(), satellite_subhalos.end());

		// If there are more than one subhalo, then return them ordered by mass in decreasing order.
		if(all.size() > 1){
			std::sort(all.begin(), all.end(), [](const SubhaloPtr &lhs, const SubhaloPtr &rhs) {
				return lhs->Mvir > rhs->Mvir;
			});
		}

		return all;
	}

	void remove_subhalo(SubhaloPtr subhalo) {

		if (subhalo == central_subhalo) {
			central_subhalo.reset();
			return;
		}

		auto it = std::find(satellite_subhalos.begin(), satellite_subhalos.end(), subhalo);
		if (it == satellite_subhalos.end()) {
			throw "subhalo not in satellites";
		}
		satellite_subhalos.erase(it);

	}

	/**
	 * The mass contained in the subhalos.
	 * This quantity should be =1 for classic SAMs, but with Rodrigo Canas work
	 * on VELOCIraptor, this quantity could be less than 1.
	 */
	float mass_fraction_subhalos;

	/** TODO: document these */
	float Vvir;
	float Mvir;
	float concentration;
	float lambda;

	float cooling_rate;

	/**
	 * The snapshot at which this halo is found
	 */
	int snapshot;
	bool main_progenitor;

	HaloPtr descendant;
	std::vector<HaloPtr> ascendants;

	/**
	 * The merger tree that holds this halo.
	 */
	MergerTreePtr merger_tree;

	void add_subhalo(const SubhaloPtr &&subhalo) {

		// Assign subhalo to proper member
		if (subhalo->subhalo_type == Subhalo::CENTRAL) {
			central_subhalo = subhalo;
		}
		else {
			satellite_subhalos.emplace_back(subhalo);
		}

		// Add subhalo mass to halo
		Mvir += subhalo->Mvir;
	}

	///
	/// Returns the number of galaxies contained in this Halo
	///
	unsigned long galaxy_count() {
		unsigned long count = 0;
		if (central_subhalo) {
			count = central_subhalo->galaxy_count();
		}
		return std::accumulate(satellite_subhalos.begin(), satellite_subhalos.end(), count,
		[](unsigned long galaxy_count, const SubhaloPtr &subhalo) {
			return galaxy_count + subhalo->galaxy_count();
		});
	}

	// Sort ascendant halos by Mvir
	std::vector<HaloPtr> ordered_ascendants(){

		if(ascendants.size()==0){
			return std::vector<HaloPtr>();
		}
		else if(ascendants.size()>1){
			std::sort(ascendants.begin(), ascendants.end(), [](const HaloPtr &lhs, const HaloPtr &rhs) {
			return lhs->Mvir > rhs->Mvir;
			});
		}
		return ascendants;
	}
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

	/**
	 * All halos contained in this merger tree, indexed by snapshot number
	 */
	std::map<int, std::vector<HaloPtr>> halos;

	void add_halo(const HaloPtr &halo) {
		halos[halo->snapshot].push_back(halo);
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
	 * mastars_burst: total stellar mass formed via bursts.
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
	std::vector<BaryonBase> mstars_burst;
	std::vector<BaryonBase> mhot_halo;
	std::vector<BaryonBase> mcold_halo;
	std::vector<BaryonBase> mejected_halo;
	std::vector<BaryonBase> mBH;
	std::vector<BaryonBase> mHI;
	std::vector<BaryonBase> mH2;
	std::vector<BaryonBase> mDM;

	std::vector<double> SFR_disk;
	std::vector<double> SFR_bulge;

	std::map<int,double> baryon_total_created;
	std::map<int,double> baryon_total_lost;

	std::vector<double> get_masses (const std::vector<BaryonBase> &B){

		std::vector<double> masses(B.size());
		std::transform(B.begin(), B.end(), masses.begin(), [](const BaryonBase &b) {
			return b.mass;
		});

		return masses;
	}

	std::vector<double> get_metals (const std::vector<BaryonBase> &B){

		std::vector<double> masses(B.size());
		std::transform(B.begin(), B.end(), masses.begin(), [](const BaryonBase &b) {
			return b.mass_metals;
		});

		return masses;
	}

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
