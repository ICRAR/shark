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

#include <cstdint>
#include <iosfwd>
#include <map>
#include <memory>
#include <vector>

#include "halo.h"
#include "mixins.h"

namespace shark {

// Forward definitions
class Halo;

using HaloPtr = std::shared_ptr<Halo>;


/**
 * A merger tree.
 *
 * A merger tree contains halos, which are indexed by snapshot,
 * and an ID to identify it.
 */
class MergerTree : public Identifiable<std::int32_t> {
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
		auto it = halos.find(snapshot);
		if (it == halos.end()) {
			return NONE;
		}
		return it->second;
	}

	std::vector<HaloPtr> &halos_at_last_snapshot()
	{
		return halos.rbegin()->second;
	}

	/**
	 * Get all the roots of this merger tree -- that is, all Halos
	 * that don't have an ascendant.
	 */
	std::vector<HaloPtr> roots()
	{
		std::vector<HaloPtr> roots;
		for (auto &snapshot_and_halos: halos) {
			for (auto &halo: snapshot_and_halos.second) {
				if (halo->ascendants.empty()) {
					roots.push_back(halo);
				}
			}
		}
		return roots;
	}

private:
	static std::vector<HaloPtr> NONE;
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