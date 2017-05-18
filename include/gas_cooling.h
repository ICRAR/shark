/*
 * gas_cooling.h
 *
 *  Created on: 17May,2017
 *      Author: clagos
 */

#ifndef INCLUDE_GAS_COOLING_H_
#define INCLUDE_GAS_COOLING_H_

#include <memory>
#include <string>
#include <vector>


#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>

#include "options.h"

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
};


class GasCooling {

public:
	GasCooling(GasCoolingParameters parameters);

	double cooling_rate(double mhot, double vvir, double mvir, double zhot);

private:
	GasCoolingParameters parameters;
	std::shared_ptr<gsl_interp2d> interp;

};

}  // namespace shark



#endif /* INCLUDE_GAS_COOLING_H_ */
