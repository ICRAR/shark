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

ReionisationParameters::ReionisationParameters(const Options &options) :
	zcut(0),
	vcut(0)
{
	options.load("reionisation.vcut",vcut);
	options.load("reionisation.zcut",zcut);
}

}
