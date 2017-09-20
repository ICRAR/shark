//
// An 2D interpolation object
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

#include <memory>
#include <vector>

#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#ifndef SHARK_INTERPOLATION_H_
#define SHARK_INTERPOLATION_H_

namespace shark {

class Interpolator {

public:

	enum InterpolatorType {
		BILINEAR = 0,
		BICUBIC,
	};

	Interpolator(std::vector<double> xvals, std::vector<double> yvals,
	             std::vector<double> zvals, InterpolatorType type = BILINEAR);

	// Copy/move constructors, destructor
	Interpolator(const Interpolator &other);
	Interpolator(Interpolator &&other);
	~Interpolator();

	double get(double x, double y) const;

private:
	std::unique_ptr<gsl_interp2d> interp2d;
	std::unique_ptr<gsl_interp_accel> xacc;
	std::unique_ptr<gsl_interp_accel> yacc;
	const gsl_interp2d_type *type;
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;

	void _init_gsl_objects();
	const gsl_interp2d_type *to_gsl(InterpolatorType type) const;
};

}  // namespace shark

#endif // SHARK_INTERPOLATION_H_