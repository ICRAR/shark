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

ReionisationParameters::ReionisationParameters(const std::string &filename) :
	Options(filename),
	zcut(0),
	vcut(0)
	{
	load("reionisation.vcut",vcut);
	load("reionisation.zcut",zcut);
}

}
