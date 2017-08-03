//
// Merger tree builder classes
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

#ifndef SHARK_TREE_BUILDER_H_
#define SHARK_TREE_BUILDER_H_

#include <memory>
#include <vector>

#include "tree_builder.h"
#include "components.h"
#include "execution.h"

namespace shark {

class TreeBuilder {

public:
	TreeBuilder(ExecutionParameters exec_params);
	std::vector<std::shared_ptr<MergerTree>> build_trees(std::vector<std::shared_ptr<Halo>> halos);

private:
	ExecutionParameters exec_params;

	void link(const std::shared_ptr<Subhalo> &subhalo, const std::shared_ptr<Subhalo> &d_subhalo,
	          const std::shared_ptr<Halo> &halo, const std::shared_ptr<Halo> &d_halo);
};

}  // namespace shark

#endif // SHARK_TREE_BUILDER_H_