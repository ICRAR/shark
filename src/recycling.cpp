/*
 * recycling.cpp
 *
 *  Created on: 14Sep.,2017
 *      Author: clagos
 */

#include "recycling.h"

namespace shark {

RecyclingParameters::RecyclingParameters(const Options &options)
{
	options.load("recycling.yield", yield, true);
	options.load("recycling.recycle", recycle, true);
	options.load("recycling.zsun", zsun, true);
}

}