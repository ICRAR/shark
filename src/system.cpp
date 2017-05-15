//
// System classes implementation
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

#include "system.h"
#include "parameters.h"

namespace shark {

static
int basic_system_evaluator(double t, const double y[], double f[], void *data) {

	/*
	 * f[0]: stellar mass of galaxy.
	 * f[1]: cold gas mass of galaxy.
	 * f[2]: hot gas mass;
	 * f[3]: metals locked in the stellar mass of galaxies.
	 * f[4]: metals locked in the cold gas mass of galaxies.
	 * f[5]: metals locked in the hot gas mass.
	 */

	BasicSystem *system = dynamic_cast<BasicSystem *>(data);

	double tau = 2.0; /*star formation timescale assumed to be 2Gyr*/
	double R = system->recycling_parameters.recycle; /*recycling fraction of newly formed stars*/
	double yield = system->recycling_parameters.yield; /*yield of newly formed stars*/
	double mcoolrate = 5e8; /*cooling rate in units of Msun/Gyr*/
	double beta = system->stellar_feedback_outflow_rate(y); /*mass loading parameter*/
	double SFR = y[1] * starformation_parameters.nu_sf; /*star formation rate assumed to be cold gas mass divided by time*/
	double zcold = y[4] / y[1]; /*cold gas metallicity*/
	double zhot = y[5] / y[2]; /*hot gas metallicity*/

	f[0] = SFR * (1-R);
	f[1] = mcoolrate - (1 - R + beta) * SFR;
	f[2] = beta * SFR - mcoolrate;
	f[3] = (1 - R) * zcold * SFR;
	f[4] = mcoolrate * zhot + SFR * (yield - (1 + beta - R) * zcold);
	f[5] = beta * zcold * SFR - mcoolrate * zhot;

	return 0;
}

BasicSystem::BasicSystem(double t0, double delta_t, double ode_solver_precision) :
	System(t0, delta_t, ode_solver_precision, basic_system_evaluator)
{
	// no-op
}

static
int basic_system_with_satellites_evaluator(double t, const double y[], double f[], void *) {

	/*
	 * f[0]: stellar mass of galaxy.
	 * f[1]: cold gas mass of galaxy.
	 * f[2]: hot gas mass;
	 * f[3]: metals locked in the stellar mass of galaxies.
	 * f[4]: metals locked in the cold gas mass of galaxies.
	 * f[5]: metals locked in the hot gas mass.
	 */

	double tau = 2.0; /*star formation timescale assumed to be 2Gyr*/
	double R = 0.5; /*recycling fraction of newly formed stars*/
	double yield = 0.029; /*yield of newly formed stars*/
	double mcoolrate = 5e8; /*cooling rate in units of Msun/Gyr*/
	double beta = 2; /*mass loading parameter*/
	double SFR = y[1] / tau; /*star formation rate assumed to be cold gas mass divided by time*/
	double zcold = y[4] / y[1]; /*cold gas metallicity*/
	double zhot = y[5] / y[2]; /*hot gas metallicity*/

	f[0] = SFR * (1-R);
	f[1] = mcoolrate - (1 - R + beta) * SFR;
	f[2] = beta * SFR - mcoolrate;
	f[3] = (1 - R) * zcold * SFR;
	f[4] = mcoolrate * zhot + SFR * (yield - (1 + beta - R) * zcold);
	f[5] = beta * zcold * SFR - mcoolrate * zhot;
	f[6] = 0;

	return 0;
}

BasicSystemWithSatellites::BasicSystemWithSatellites(double t0, double delta_t, double ode_solver_precision) :
	System(t0, delta_t, ode_solver_precision, basic_system_with_satellites_evaluator)
{
	// no-op
}

}  // namespace shark
