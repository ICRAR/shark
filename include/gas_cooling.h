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

#include "agn_feedback.h"
#include "components.h"
#include "dark_matter_halos.h"
#include "environment.h"
#include "execution.h"
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

	using tables_idx = std::map<double, std::string>;

	enum LambdaCoolingModel {
		CLOUDY = 0,
		SUTHERLAND
	};

	enum CoolingModel {
		CROTON06 = 0,
		BENSON10
	};

	explicit GasCoolingParameters(const Options &options);

	double pre_enrich_z = 1e-7;
	double tau_cooling = 1;
	bool limit_fbar = false;
        double rcore = 0.01;

	LambdaCoolingModel lambdamodel = CLOUDY;
	CoolingModel model = CROTON06;

	//cooling tables
	CoolingTable cooling_table; //these should be an array of parameters.


private:
	tables_idx find_tables(const std::string &cooling_tables_dir);
	void load_tables(const std::string &cooling_tables_dir, const tables_idx &metallicity_tables);
};


class GasCooling {

public:
	GasCooling(GasCoolingParameters parameters,
			StarFormationParameters params_sf,
			ExecutionParameters exec_params,
			ReionisationPtr reionisation,
			CosmologyPtr cosmology,
			AGNFeedbackPtr agnfeedback,
			DarkMatterHalosPtr darkmatterhalos,
			ReincorporationPtr reincorporation,
			EnvironmentPtr environment);

	double cooling_rate(Subhalo &subhalo, Galaxy &galaxy, double z, double deltat);
	double cooling_time(double Tvir, double logl, double nh_density);
	double mean_density(double mhot, double rvir);
	double cooling_radius(double mhot, double rvir, double tcharac, double logl, double Tvir);
	double density_shell(double mhot, double rvir, double r);
	double cooling_luminosity(double logl, double rcool, double rvir, double mhot);
	bool quasi_hydrostatic_halo(double mhot, double lambda, double nh_density,
			double mass, double Tvir, double redshift);

private:

	GasCoolingParameters parameters;
	StarFormationParameters params_sf;
	ExecutionParameters exec_params;
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
