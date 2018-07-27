//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2018
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

#ifndef INCLUDE_GAS_COOLING_H_
#define INCLUDE_GAS_COOLING_H_

#include <map>
#include <memory>
#include <string>
#include <vector>


#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "agn_feedback.h"
#include "components.h"
#include "dark_matter_halos.h"
#include "environment.h"
#include "interpolator.h"
#include "options.h"
#include "reincorporation.h"
#include "reionisation.h"
#include "star_formation.h"

namespace shark {
/**
 * An element of the cooling tables
 */
class CoolingTable {

public:

	void add_metallicity_measurements(double zmetal, const std::map<double, double> &records);

	std::vector<double> get_temperatures();
	std::vector<double> get_metallicities();
	std::vector<double> get_lambda();

private:
	std::map<double, std::map<double, double>> _table;
};

class GasCoolingParameters {

public:

	typedef std::map<double, std::string> tables_idx;

	enum LambdaCoolingModel {
		CLOUDY = 0,
		SUTHERLAND
	};

	enum CoolingModel {
		CROTON06 = 0,
		BENSON10
	};

	GasCoolingParameters(const Options &options);

	double rcore = 0;
	double pre_enrich_z = 1e-7;
	double tau_cooling = 1;

	LambdaCoolingModel lambdamodel = CLOUDY;
	CoolingModel model = CROTON06;

	//cooling tables
	CoolingTable cooling_table {}; //these should be an array of parameters.


private:
	tables_idx find_tables(const std::string &cooling_tables_dir);
	void load_tables(const std::string &cooling_tables_dir, const tables_idx &metallicity_tables);
};


class GasCooling {

public:
	GasCooling(GasCoolingParameters parameters,
			StarFormationParameters params_sf,
			const ReionisationPtr &reionisation,
			const CosmologyPtr &cosmology,
			const AGNFeedbackPtr &agnfeedback,
			const DarkMatterHalosPtr &darkmatterhalos,
			const ReincorporationPtr &reincorporation,
			const EnvironmentPtr &environment);

	double cooling_rate(Subhalo &subhalo, Galaxy &galaxy, double z, double deltat);
	double cooling_time(double Tvir, double logl, double nh_density);
	double mean_density(double mhot, double rvir);
	double cooling_radius(double mhot, double rvir, double tcharac, double logl, double Tvir);
	double density_shell(double mhot, double rvir, double r);
	double cooling_luminosity(double logl, double rcool, double rvir, double mhot);
	double disk_size_cooling(Subhalo &subhalo);

private:

	GasCoolingParameters parameters;
	StarFormationParameters params_sf;
	ReionisationPtr reionisation;
	CosmologyPtr cosmology;
	AGNFeedbackPtr agnfeedback;
	DarkMatterHalosPtr darkmatterhalos;
	ReincorporationPtr reincorporation;
	EnvironmentPtr environment;
	Interpolator cooling_lambda_interpolator;

};

}  // namespace shark

#endif /* INCLUDE_GAS_COOLING_H_ */
