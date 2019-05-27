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

#include <algorithm>

#include "total_baryon.h"

namespace shark {

std::vector<double> TotalBaryon::get_masses (const std::vector<BaryonBase> &B) const
{
	std::vector<double> masses(B.size());
	std::transform(B.begin(), B.end(), masses.begin(), [](const BaryonBase &b) {
		return b.mass;
	});
	return masses;
}

std::vector<double> TotalBaryon::get_metals (const std::vector<BaryonBase> &B) const
{
	std::vector<double> masses(B.size());
	std::transform(B.begin(), B.end(), masses.begin(), [](const BaryonBase &b) {
		return b.mass_metals;
	});
	return masses;
}

}  // namespace shark