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
 * Baryon-related classes and functionality
 */

#ifndef INCLUDE_BARYON_H_
#define INCLUDE_BARYON_H_

namespace shark {

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

	BaryonBase &operator-=(const BaryonBase &b) {
		mass -= b.mass;
		mass_metals -= b.mass_metals;
		return *this;
	}

	friend BaryonBase &operator+(BaryonBase &lhs, const BaryonBase &rhs)
	{
		lhs += rhs;
		return lhs;
	}

	friend BaryonBase &operator-(BaryonBase &lhs, const BaryonBase &rhs)
	{
		lhs -= rhs;
		return lhs;
	}

	void restore_baryon()
	{
		mass = 0;
		mass_metals = 0;
	}

};

/**
 * A baryon that rotates, and therefore has angular momentum.
 */
class RotatingBaryonBase : public BaryonBase {

public:

	/**
	 * Specific angular momentum
	 */
	float sAM = 0;

	/**
	 * @return The angular momentum of this baryon
	 */
	float angular_momentum() const
	{
		return mass * sAM;
	}

	void restore_baryon()
	{
		BaryonBase::restore_baryon();
		sAM = 0;
	}

	friend RotatingBaryonBase &operator+(RotatingBaryonBase &lhs, const RotatingBaryonBase &rhs)
	{
		lhs += rhs;
		return lhs;
	}

};

/**
 * A common baryon component.
 * Note that black holes are not baryon components as they use their own class.
 */
class Baryon : public RotatingBaryonBase {
public:

	/**
	 * A scale radius
	 */
	float rscale = 0;

	friend Baryon &operator+(Baryon &lhs, const Baryon &rhs)
	{
		lhs += rhs;
		return lhs;
	}

	friend Baryon &operator-(Baryon &lhs, const BaryonBase &rhs)
	{
		lhs -= rhs;
		return lhs;
	}

	void restore_baryon()
	{
		RotatingBaryonBase::restore_baryon();
		rscale = 0;
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

	/** macc_sb: accretion rate onto the black hole during galaxy starbursts. */
	float macc_sb = 0;

	/** massembly: mass that has been brought from BH-BH mergers. */
	float massembly = 0;

	/** spin: black hole spin. */
	float spin = 0;

};

}  // namespace shark


#endif /* INCLUDE_BARYON_H_ */
