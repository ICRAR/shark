//
// Interpolation2D class implementation
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

#include <algorithm>
#include <sstream>

#include "exceptions.h"
#include "interpolator.h"

namespace shark {

Interpolator::Interpolator(std::vector<double> xvals, std::vector<double> yvals,
		std::vector<double> zvals, InterpolatorType type) :
	spline2d(nullptr),
	xacc(nullptr),
	yacc(nullptr),
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
	spline2d.reset(gsl_spline2d_alloc(type, x.size(), y.size()));
	gsl_spline2d_init(spline2d.get(), x.data(), y.data(), z.data(), x.size(), y.size());
	xacc.reset(gsl_interp_accel_alloc());
	yacc.reset(gsl_interp_accel_alloc());
}

Interpolator::Interpolator(const Interpolator &other) :
	spline2d(nullptr),
	xacc(nullptr),
	yacc(nullptr),
	type(other.type),
	x(other.x), y(other.y), z(other.z)
{
	_init_gsl_objects();
}

Interpolator::Interpolator(Interpolator &&other) :
	spline2d(nullptr),
	xacc(nullptr),
	yacc(nullptr),
	type(other.type),
	x(std::move(other.x)), y(std::move(other.y)), z(std::move(other.z))
{
	std::swap(other.spline2d, spline2d);
	std::swap(other.xacc, xacc);
	std::swap(other.yacc, yacc);
}

Interpolator::~Interpolator()
{
	if (xacc) {
		gsl_interp_accel_free(xacc.release());
	}
	if (yacc) {
		gsl_interp_accel_free(yacc.release());
	}
	if (spline2d) {
		gsl_spline2d_free(spline2d.release());
	}
}

double Interpolator::get(double x, double y) const
{
	return gsl_spline2d_eval(spline2d.get(), x, y, xacc.get(), yacc.get());
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