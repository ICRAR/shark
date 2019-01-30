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
 */

#include <cmath>
#include <fstream>
#include <map>
#include <tuple>

#include "logging.h"
#include "reionisation.h"
#include "utils.h"

namespace shark {

ReionisationParameters::ReionisationParameters(const Options &options)
{
	options.load("reionisation.vcut", vcut, true);
	options.load("reionisation.zcut", zcut, true);
	options.load("reionisation.alpha_v", alpha_v);
	options.load("reionisation.model", model, true);

}

template <>
ReionisationParameters::ReionisationModel
Options::get<ReionisationParameters::ReionisationModel>(const std::string &name, const std::string &value) const {
	auto lvalue = lower(value);
	if (lvalue == "lacey16") {
		return ReionisationParameters::LACEY16;
	}
	else if (lvalue == "sobacchi13") {
		return ReionisationParameters::SOBACCHI13;
	}

	std::ostringstream os;
	os << name << " option value invalid: " << value << ". Supported values are Lacey16 and Sobacchi13";
	throw invalid_option(os.str());
}

Reionisation::Reionisation(const ReionisationParameters &parameters) :
	parameters(parameters)
{
	// no-op
}

Reionisation::~Reionisation() = default;

bool Lacey16Reionisation::reionised_halo(double v, double z) const
{
	auto &params = get_reionisation_params();
	return (v < params.vcut && z < params.zcut);
}

bool Sobacchi13Reionisation::reionised_halo(double v, double z) const
{
	using std::pow;
	auto &params = get_reionisation_params();
	double vthresh = params.vcut * pow(1.0 + z, params.alpha_v) * pow((1.0 - pow((1.0 + z) /(1.0 + params.zcut), 2.0)), 0.833);
	return v < vthresh;
}

}