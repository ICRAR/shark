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

#include "options.h"
#include "components.h"

namespace shark {
/**
 * An element of the cooling tables
 */
struct CoolingTable {
	std::vector<double> log10lam;
	std::vector<double> log10temp;
	std::vector<double> zmetal;
};

class GasCoolingParameters : public Options {

public:

	typedef std::map<double, std::string> tables_idx;

	enum LambdaCoolingModel {
		CLOUDY = 0,
		SUTHERLAND
	};

	enum CoolingModel {
		CROTON06 = 0,
		GALFORM
	};

	GasCoolingParameters(const std::string &filename);

	double rcore;
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
	GasCooling(GasCoolingParameters parameters);

	double cooling_rate(std::shared_ptr<Subhalo> &subhalo, double deltat);

private:
	GasCoolingParameters parameters;
	std::shared_ptr<gsl_interp2d> interp;

};

}  // namespace shark



#endif /* INCLUDE_GAS_COOLING_H_ */
