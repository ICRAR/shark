//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2018
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
 * Shark forward class declarations
 */

#ifndef SHARKFWD_H_
#define SHARKFWD_H_

#include <memory>

namespace shark {

// Components and their pointer types
class Baryon;
class Galaxy;
class Subhalo;
class Halo;
class MergerTree;

typedef std::shared_ptr<Galaxy> GalaxyPtr;
typedef std::shared_ptr<Subhalo> SubhaloPtr;
typedef std::shared_ptr<Halo> HaloPtr;
typedef std::shared_ptr<MergerTree> MergerTreePtr;

// Physical processes, their parameters, and their pointer types
class AGNFeedback;
class AGNFeedbackParameters;
typedef std::shared_ptr<AGNFeedback> AGNFeedbackPtr;

class Cosmolog;
class CosmologicalParameters;
typedef std::shared_ptr<Cosmology> CosmologyPtr;

class DarkMatterHalos;
class DarkMatterHaloParameters;
typedef std::shared_ptr<DarkMatterHalos> DarkMatterHalosPtr;

class DiskInstability;
class DiskInstabilityParameters;

class GalaxyMerger;
class GalaxyMergerParameters;

class Reionisation;
class ReionisationParameters;
typedef std::shared_ptr<Reionisation> ReionisationPtr;

class Reincorporation;
class ReincorporationParameters;
typedef std::shared_ptr<Reincorporation> ReincorporationPtr;

// High-level physical models
class GalaxyCreator;
template <int N>
class PhysicalModel;
class BasicPhysicalModel;

// I/O classes
class GalaxyWriter;
class SURFSReader;

// Others
class Options;
class Simulation;
class SimulationParameters;
class TreeBuilder;

} // namespace shark

#endif /* SHARKFWD_H_ */