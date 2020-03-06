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

/**
 * @file Traits for shark component types
 */

#ifndef INCLUDE_COMPONENTS_TRAITS_H_
#define INCLUDE_COMPONENTS_TRAITS_H_

#include <type_traits>

#include "components.h"

namespace shark {

/**
 * Defines a static @pre value boolean member if @pre T is one of the shark
 * component types, including the pointer types.
 *
 * @tparam T The type to inspect
 */
template <typename T>
struct is_component : public std::false_type {};
template <>
struct is_component<Galaxy> : public std::true_type {};
template <>
struct is_component<Subhalo> : public std::true_type {};
template <>
struct is_component<Halo> : public std::true_type {};
template <>
struct is_component<MergerTree> : public std::true_type {};
template <>
struct is_component<GalaxyPtr> : public std::true_type {};
template <>
struct is_component<SubhaloPtr> : public std::true_type {};
template <>
struct is_component<HaloPtr> : public std::true_type {};
template <>
struct is_component<MergerTreePtr> : public std::true_type {};

/**
 * Defines a static @pre value boolean member if @pre T is one of the shark
 * pointer types.
 *
 * @tparam T the type to inspect
 */
template <typename T>
struct is_component_pointer : public std::false_type {};
template <>
struct is_component_pointer<GalaxyPtr> : public std::true_type {};
template <>
struct is_component_pointer<SubhaloPtr> : public std::true_type {};
template <>
struct is_component_pointer<HaloPtr> : public std::true_type {};
template <>
struct is_component_pointer<MergerTreePtr> : public std::true_type {};

/**
 * Defines a static @pre value boolean member if @pre T is one of the shark
 * smart pointer types.
 *
 * @tparam T the type to inspect
 */
template <typename T>
struct is_component_smart_pointer : public std::false_type {};
template <>
struct is_component_smart_pointer<SubhaloPtr> : public std::true_type {};
template <>
struct is_component_smart_pointer<HaloPtr> : public std::true_type {};
template <>
struct is_component_smart_pointer<MergerTreePtr> : public std::true_type {};

/// Like std's in C++14, but in the shark namespace and for C++11
template <bool B, class T=void>
using enable_if_t = typename std::enable_if<B, T>::type;

/// Like enable_if_t but specific for shark component checking
template <class ComponentT, class T=void>
using enable_if_component_t = enable_if_t<is_component<ComponentT>::value, T>;

namespace detail {

	template <typename T, bool smart_pointer=false>
	struct id_type_trait {
		using type = typename T::id_t;
	};

	template <typename T>
	struct id_type_trait<T, true> {
		using type = typename T::element_type::id_t;
	};
}

/**
 * The ID type used by the shark type @pre T, which can be any of the shark
 * component types
 */
template <typename T, typename Void=enable_if_component_t<T>>
using id_type = typename detail::id_type_trait<T, is_component_smart_pointer<T>::value>::type;

}  // namespace shark

#endif /* INCLUDE_COMPONENTS_TRAITS_H_ */