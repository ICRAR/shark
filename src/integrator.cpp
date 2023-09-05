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
 * Integrator class implementation
 */

#include <string>

#include "gsl_utils.h"
#include "integrator.h"
#include "logging.h"


namespace shark {

Integrator::Integrator(size_t max_intervals) :
	max_intervals(max_intervals),
	num_intervals(0)
{
	init_gsl_objects();
}

Integrator::Integrator(const Integrator &other) :
	max_intervals(other.max_intervals),
	num_intervals(other.num_intervals)
{
	init_gsl_objects();
}

void Integrator::init_gsl_objects()
{
	workspace.reset(gsl_integration_workspace_alloc(max_intervals));
}

struct integration_data_t {
	Integrator::func_t func;
	void *params;
	int status;
};

extern "C" double gsl_integration_target(double x, void * params)
{
	// Catch C++ exceptions, set our internal status flag to GSL_EBADFUNC,
	// and report NaN as a result hoping that integration will quickly fail
	auto integration_data = static_cast<integration_data_t *>(params);
	try {
		return integration_data->func(x, integration_data->params);
	} catch (std::exception &e) {
		std::string msg = "Error while evaluating integration function: ";
		msg += e.what();
		LOG(error) << msg;
		::gsl_error(msg.c_str(), __FILE__, __LINE__, GSL_EBADFUNC);
		integration_data->status = GSL_EBADFUNC;
		return GSL_NAN;
	} catch (...) {
		auto msg = "Unexpected error while evaluating integration function";
		LOG(error) << msg << ", aborting integration";
		::gsl_error(msg, __FILE__, __LINE__, GSL_EBADFUNC);
		integration_data->status = GSL_EBADFUNC;
		return GSL_NAN;
	}
}

double Integrator::integrate(func_t f, void *params, double from, double to, double epsabs, double epsrel)
{
	// User-provided C++ functions for integration with GSL (parameter `f` here)
	// might raise exceptions, which do not necessarily play well with the GSL
	// C-based stack (i.e., exceptions might not propagate correctly).
	// To overcome this we first pack both the user-provided function and
	// parameters into our own integration_data_t object, and then point GSL
	// to integrate the gsl_integration_target function instead. The latter
	// receives the integration_data_t object as parameter, unpacks the
	// user-provided function and parameters, and executes that in a try/catch
	// block, returning NaN in case of an error.
	//
	// In GSL the user-provided functions also don't have a clear way of
	// signaling errors, so even though we return NaN on failure there is still
	// a change that the integrator won't see that as a failure. We thus include
	// our own status flag in the integration_data_t object, which is set to
	// GSL_EBADFUNC in the catch block of gsl_integration_target. We then check
	// this flag for any errors that might not have been caught by the GSL
	// integration routine.
	integration_data_t integration_data;
	integration_data.func = f;
	integration_data.params = params;
	integration_data.status = GSL_SUCCESS;
	gsl_function F;
	F.function = gsl_integration_target;
	F.params = &integration_data;

	double result, abserr;
	gsl_invoke(gsl_integration_qag, &F, from, to, epsabs, epsrel, max_intervals,
	           GSL_INTEG_GAUSS15, workspace.get(), &result, &abserr);
	if (integration_data.status != GSL_SUCCESS) {
		throw to_gsl_error(integration_data.status);
	}

	num_intervals += workspace->size;
	return result;
}

std::size_t Integrator::get_num_intervals()
{
	return num_intervals;
}

void Integrator::reset_num_intervals()
{
	num_intervals = 0;
}

}  // namespace shark
