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

	double zcut;
	double vcut;

};
}

#endif //SHARK_REIONISATION_H_
