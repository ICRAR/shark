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
 * Interpolation2D class implementation
 */

#include <algorithm>
#include <sstream>

#include "exceptions.h"
#include "gsl_utils.h"
#include "interpolator.h"

namespace shark {

Interpolator::Interpolator(std::vector<double> xvals, std::vector<double> yvals,
		std::vector<double> zvals, InterpolatorType type) :
	type(to_gsl(type)),
	x(std::move(xvals)), y(std::move(yvals)), z(std::move(zvals))
{
	if (x.size() * y.size() != z.size()) {
		std::ostringstream os;
		os << "Grid size (" << x.size() << "x" << y.size() << " = " << x.size() * y.size();
		os << ") does not correspond with values size (" << z.size() << ")";
		throw invalid_argument(os.str());
	}
	_init_gsl_objects();
}

void Interpolator::_init_gsl_objects()
{
	interp2d.reset(gsl_interp2d_alloc(type, x.size(), y.size()));
	gsl_invoke(gsl_interp2d_init, interp2d.get(), x.data(), y.data(), z.data(), x.size(), y.size());
	xacc.reset(gsl_interp_accel_alloc());
	yacc.reset(gsl_interp_accel_alloc());
}

Interpolator::Interpolator(const Interpolator &other) :
	type(other.type),
	x(other.x), y(other.y), z(other.z)
{
	_init_gsl_objects();
}

double Interpolator::get(double x, double y) const
{
	// Truncate x and y to min/max values in ranges
	x = std::min(std::max(x, this->x.front()), this->x.back());
	y = std::min(std::max(y, this->y.front()), this->y.back());
	return gsl_interp2d_eval_extrap(interp2d.get(),
			this->x.data(), this->y.data(), this->z.data(),
			x, y, xacc.get(), yacc.get());
}

const gsl_interp2d_type *Interpolator::to_gsl(InterpolatorType type) const
{
	if (type == BILINEAR) {
		return gsl_interp2d_bilinear;
	}
	else if (type == BICUBIC) {
		return gsl_interp2d_bicubic;
	}

	std::ostringstream os;
	os << "Invalid interpolation type: " << type;
	throw invalid_argument(os.str());
}

}  // namespace shark
