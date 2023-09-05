//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2019
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

/// @file Common algorithms for shark structural components

#ifndef INCLUDE_COMPONENTS_ALGORITHMS_H_
#define INCLUDE_COMPONENTS_ALGORITHMS_H_

#include <algorithm>

#include "components/traits.h"

namespace shark {

/**
 * Gets the ID of component x
 *
 * @param x a shark component
 * @return The ID of the shark component
 */
template <typename T>
typename std::enable_if<!is_component_pointer<T>::value, id_type<T>>::type
get_id(const T &x)
{
	return x.id;
}

template <typename T>
typename std::enable_if<is_component_pointer<T>::value, id_type<T>>::type
get_id(const T &x)
{
	return x->id;
}

/**
 * Template function that compares two component objects (of the same type) by
 * their IDs.
 * @param x The first component
 * @param y The second component
 * @return true if the ID of the first component is less than the ID of the
 * second, false otherwise
 */
template <typename T>
enable_if_component_t<T, bool>
less_by_id(const T &x, const T &y)
{
	return get_id(x) < get_id(y);
}

/**
 * Sorts the given collection of identifiable components by their IDs in
 * ascending order.
 * @param identifiables The objects to sort by ID
 */
template <typename Range>
enable_if_component_t<typename Range::value_type>
sort_by_id(Range &identifiables)
{
	using value_type = typename Range::value_type;
	std::sort(identifiables.begin(), identifiables.end(), less_by_id<value_type>);
}

}  // namespace shark

#endif /* INCLUDE_COMPONENTS_ALGORITHMS_H_ */