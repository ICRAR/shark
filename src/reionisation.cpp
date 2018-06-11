/*
 * reionisation.cpp
 *
 *  Created on: 14Jun.,2017
 *      Author: clagos
 */


#include <cmath>
#include <fstream>
#include <map>
#include <tuple>

#include "reionisation.h"
#include "logging.h"

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
	if ( value == "Lacey16" ) {
		return ReionisationParameters::LACEY16;
	}
	else if ( value == "Sobacchi13" ) {
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

Reionisation::~Reionisation()
{
}

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