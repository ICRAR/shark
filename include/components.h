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
 * Forward definition of commonly-used structural component classes
 */

#ifndef SHARK_COMPONENTS_H_
#define SHARK_COMPONENTS_H_

#include <cstdint>
#include <memory>

namespace shark {

using galaxy_id_t = int;
using subhalo_id_t = std::int64_t;
using halo_id_t = std::int64_t;
using merger_tree_id_t = std::int32_t;

class Galaxy;
class Subhalo;
class Halo;
class MergerTree;
class TotalBaryon;

using GalaxyPtr = Galaxy *;
using ConstGalaxyPtr = const Galaxy *;
using SubhaloPtr = std::shared_ptr<Subhalo>;
using HaloPtr = std::shared_ptr<Halo>;
using MergerTreePtr = std::shared_ptr<MergerTree>;
using galaxies_size_type = std::size_t;

} // namespace shark

#endif // SHARK_COMPONENTS_H_
