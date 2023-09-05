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
 * Merger tree-related classes and functionality
 */

#ifndef INCLUDE_MERGER_TREE_H_
#define INCLUDE_MERGER_TREE_H_

#include <algorithm>
#include <iosfwd>
#include <map>
#include <memory>
#include <numeric>
#include <vector>

#include "halo.h"
#include "mixins.h"
#include "ranges.h"

namespace shark {

// Forward definitions
class Halo;

using HaloPtr = std::shared_ptr<Halo>;

namespace detail {

class is_root {
public:
	bool operator()(const HaloPtr &halo) const
	{
		return halo->ascendants.empty();
	}
};

}  // namespace detail

/**
 * A merger tree.
 *
 * A merger tree contains halos, which are indexed by snapshot,
 * and an ID to identify it.
 */
class MergerTree : public Identifiable<merger_tree_id_t> {
public:

	using Identifiable::Identifiable;

	using root_subrange = range_filter<std::vector<HaloPtr>, detail::is_root>;
	using snapshot_subrange = range<std::vector<HaloPtr>::iterator>;

	void add_halo(const HaloPtr &halo) {
		if (last_snapshot < halo->snapshot) {
			last_snapshot = halo->snapshot;
		}
		halos.push_back(halo);
	}

	snapshot_subrange halos_at(int snapshot)
	{
		return {
			std::lower_bound(halos.begin(), halos.end(), snapshot, by_snapshot{}),
			std::upper_bound(halos.begin(), halos.end(), snapshot, by_snapshot{})
		};
	}

	snapshot_subrange halos_at_last_snapshot()
	{
		return halos_at(last_snapshot);
	}

	void consolidate()
	{
		std::sort(halos.begin(), halos.end(), by_snapshot{});
	}

	/**
	 * Get all the roots of this merger tree -- that is, all Halos
	 * that don't have an ascendant.
	 */
	root_subrange roots()
	{
		return {halos};
	}

	/**
	 * Calculate and return the number of galaxies contained in this merger tree
	 * @return The number of galaxies contained in this merger tree
	 */
	std::size_t galaxy_count()
	{
		return accumulate(halos.begin(), halos.end(), size_t(0), [](const size_t n, const HaloPtr &halo) {
			return n + halo->galaxy_count();
		});
	}

	/// All halos contained in this merger tree
	std::vector<HaloPtr> halos;
	/// Last snapshot included by this MergerTree
	int last_snapshot = -1;
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

#endif /* INCLUDE_MERGER_TREE_H_ */