//
// Halo evolution routines
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

#ifndef SHARK_EVOLVE_HALOS_H_
#define SHARK_EVOLVE_HALOS_H_

#include "components.h"
#include "physical_model.h"

namespace shark {

void populate_halos(BasicPhysicalModel &physical_model, std::shared_ptr<Halo> halo, int snapshot, double z, double delta_t);

}  // namespace shark

#endif // SHARK_EVOLVE_HALOS_H_
