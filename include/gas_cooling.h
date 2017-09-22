/*
 * gas_cooling.h
 *
 *  Created on: 17May,2017
 *      Author: clagos
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
#include "interpolator.h"
#include "options.h"
#include "reincorporation.h"
#include "reionisation.h"

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

	double rcore;
	double pre_enrich_z;

	LambdaCoolingModel lambdamodel;
	CoolingModel model;
	//cooling tables
	CoolingTable cooling_table; //these should be an array of parameters.


private:
	tables_idx find_tables(const std::string &cooling_tables_dir);
	void load_tables(const std::string &cooling_tables_dir, const tables_idx &metallicity_tables);
};


class GasCooling {

public:
	GasCooling(GasCoolingParameters parameters, ReionisationParameters reio_parameters, std::shared_ptr<Cosmology> cosmology, std::shared_ptr<AGNFeedback> agnfeedback, std::shared_ptr<DarkMatterHalos> darkmatterhalos, std::shared_ptr<Reincorporation> reincorporation);

	double cooling_rate(Subhalo &subhalo, double z, double deltat);
	double cooling_time(double Tvir, double logl, double nh_density);
	double mean_density(double mhot, double rvir);
	double cooling_radius(double mhot, double rvir, double tcharac, double logl, double Tvir);
	double density_shell(double mhot, double rvir, double r);
	double cooling_luminosity(double logl, double rcool, double rvir, double mhot);
	double disk_size_cooling(Subhalo &subhalo);

private:

	ReionisationParameters reio_parameters;
	GasCoolingParameters parameters;
	std::shared_ptr<Cosmology> cosmology;
	std::shared_ptr<AGNFeedback> agnfeedback;
	std::shared_ptr<DarkMatterHalos> darkmatterhalos;
	std::shared_ptr<Reincorporation> reincorporation;
	Interpolator cooling_lambda_interpolator;

};

}  // namespace shark



#endif /* INCLUDE_GAS_COOLING_H_ */
