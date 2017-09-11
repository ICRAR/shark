//
// Useful mix-in classes
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

#ifndef SHARK_MIXINS_H_
#define SHARK_MIXINS_H_

#include "cmath"

namespace shark {

/**
 * A class with three components representing values for the X, Y and Z
 * coordinates.
 *
 * @tparam T The type used to store the underlying values
 */
template <typename T>
class xyz {

public:

	/**
	 * The value in the X coordinate
	 */
	T x;

	/**
	 * The value in the Y coordinate
	 */
	T y;

	/**
	 * The value in the Z coordinate
	 */
	T z;

	T norm(){
		T v2= std::pow(x,2)+std::pow(y,2)+std::pow(z,2);
		return std::pow(v2,0.5);
	}

};

/**
 * An object that has a position and a velocity.
 *
 * @tparam T The type used to store the position and velocity
 */
template <typename T>
class Spatial {

public:
	xyz<T> position;
	xyz<T> velocity;
};


/**
 * Base class for all objects that require unique identification.
 *
 * @param T the type used to represent the ID
 */
template <typename T>
class Identifiable {
public:

	typedef T id_t;

	/**
	 * The ID of this object
	 */
	id_t id;

};

}  // namespace shark

#endif // SHARK_MIXINS_H_
