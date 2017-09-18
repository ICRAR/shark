/*
 * write_output.cpp

 *
 *  Created on: 18Sep.,2017
 *      Author: clagos
 */

#include "exceptions.h"
#include "logging.h"
#include "write_output.h"

using namespace std;

namespace shark {

WriteOutput::WriteOutput(ExecutionParameters exec_params):
	exec_params(exec_params)
	{
		//no-opt
	}

void WriteOutput::write_galaxies(int snapshot, std::vector<HaloPtr> halos){

	for (auto &halo: halos){

	}


}

}// namespace shark

