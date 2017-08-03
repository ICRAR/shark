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

#include <map>
#include <memory>
#include <vector>

#include "mixins.h"

namespace shark {

// Forward-defines
class Halo;
class MergerTree;

/**
 * The common base for all baryon component types.
 */
class BaryonBase {

public:

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
	 * TODO: add description of this parameter
	 */
	float macc;
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
	 * An enumeration of types of galaxies
	 */
	enum galaxy_type_t {
		CENTRAL = 0,
		SATELLITE,
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
	BlackHole smbh;
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
	std::vector<double> deltat;
	std::vector<double> temp;
	std::vector<double> mass;
	std::vector<double> tcooling;
};


/** This class defines what a subhalo is. In an abstract sense, a subhalo is the quantum units of how dark matter clusters. Subhalos can
 * coexist with other subhalos in the same halo. A subhalo con also host 0, 1 or more galaxies depending on how we allow galaxies to populate
 * subhalos.
 */
class Subhalo : public Identifiable<long>, public Spatial<float> {

public:

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
	 * The ID of the descendant of this subhalo.
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
	 * A pointer to the descendant of this subhalo.
	 * If this pointer is set then descendant_id and descendant_subhalo are
	 * meaningless.
	 */
	std::shared_ptr<Subhalo> descendant;

	/**
	 * The list of galaxies in this subhalo.
	 */
	std::vector<std::shared_ptr<Galaxy>> galaxies;

	/**
	 * The subhalo type
	 */
	subhalo_type_t subhalo_type;

	/** TODO: Properly document these */
	int haloID; /*Which halos does this subhalo belong to*/
	float Vvir;
	float Mvir;
	float L[3];
	float Vcirc;
	float Concentration;

	/**
	 * Hot gas component of the halo and outside the galaxies that is
	 * allowed to cool down and/or fall onto the galaxy.
	 */
	Baryon hot_halo_gas;


	/**
	 * This component saves que information of the virial temperature, total halo gas and cooling time history.
	 */
	CoolingSubhaloTracking cooling_subhalo_tracking;

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
	std::vector<std::shared_ptr<Subhalo>> ascendants;

	/**
	 * The accretion rate onto the subhalo. This information comes from the merger tree
	 */
	std::vector<float> accretion_rate;

	/**
	 * The halo that holds this subhalo.
	 */
	std::shared_ptr<Halo> host_halo;
};

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

	Halo(long halo_id, int snapshot) :
		central_subhalo(),
		satellite_subhalos(),
		mass_fraction_subhalos(-1),
		Vvir(0),
		Mvir(0),
		snapshot(snapshot)
	{
		// no-op
		id = halo_id;
	}

	/**
	 * The central subhalo
	 */
	std::shared_ptr<Subhalo> central_subhalo;

	/**
	 * The subhalos contained in this halo
	 */
	std::vector<std::shared_ptr<Subhalo>> satellite_subhalos;

	std::vector<std::shared_ptr<Subhalo>> all_subhalos() {

		std::vector<std::shared_ptr<Subhalo>> all;

		if (central_subhalo) {
			all.push_back(central_subhalo);
		}
		all.insert(all.end(), satellite_subhalos.begin(), satellite_subhalos.end());

		return all;
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

	/**
	 * The snapshot at which this halo is found
	 */
	int snapshot;

	std::shared_ptr<Halo> descendant;
	std::vector<std::shared_ptr<Halo>> ascendants;

	/**
	 * The merger tree that holds this halo.
	 */
	std::shared_ptr<MergerTree> merger_tree;

	void add_subhalo(const std::shared_ptr<Subhalo> &subhalo) {

		// Assign subhalo to proper member
		if (subhalo->subhalo_type == Subhalo::CENTRAL) {
			central_subhalo = subhalo;
		}
		else {
			satellite_subhalos.push_back(subhalo);
		}

		// Add subhalo mass to halo
		Mvir += subhalo->Mvir;
	}

};

/**
 * A merger tree.
 *
 * A merger tree contains halos, which are indexed by snapshot,
 * and an ID to identify it.
 */
class MergerTree {
public:

	/**
	 * All halos contained in this merger tree, indexed by snapshot number
	 */
	std::map<int, std::vector<std::shared_ptr<Halo>>> halos;

	void add_halo(const std::shared_ptr<Halo> &halo) {
		halos[halo->snapshot].push_back(halo);
	}
};

}  // namespace shark

#endif // SHARK_COMPONENTS_H_
