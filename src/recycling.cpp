/*
 * recycling.cpp
 *
 *  Created on: 14Sep.,2017
 *      Author: clagos
 */

#include <cmath>

#include "recycling.h"

namespace shark {

RecyclingParameters::RecyclingParameters(const Options &options) :
	yield(0),
	recycle(0),
	zsun(0)
{
	options.load("recycling.yield", yield, true);
	options.load("recycling.recycle", recycle, true);
	options.load("recycling.zsun", zsun, true);

}

}


