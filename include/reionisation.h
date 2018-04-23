/*
 * reionisation.h
 *
 *  Created on: 14Jun.,2017
 *      Author: clagos
 */

#ifndef SHARK_REIONISATION_H_
#define SHARK_REIONISATION_H_

#include <vector>
#include <string>

#include "options.h"

namespace shark {

class ReionisationParameters {

public:
	ReionisationParameters(const Options &options);

	enum ReionisationModel {
		LACEY16 = 0,
		SOBACCHI13
	};

	double zcut;
	double vcut;
	double alpha_v;

	ReionisationModel model;


};


class Reionisation{

public:
	Reionisation(ReionisationParameters parameters);

	bool reionised_halo (double v, double z);

private:

	ReionisationParameters parameters;

};

}//end namespace shark


#endif //SHARK_REIONISATION_H_
