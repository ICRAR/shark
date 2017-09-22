/*
 * write_output.h
 *
 *  Created on: 18Sep.,2017
 *      Author: clagos
 */

#ifndef INCLUDE_WRITE_OUTPUT_H_
#define INCLUDE_WRITE_OUTPUT_H_

#include <map>
#include <memory>

#include "cosmology.h"
#include "execution.h"
#include "simulation.h"
#include "star_formation.h"

namespace shark {

class WriteOutput{

public:

	WriteOutput(ExecutionParameters exec_params, CosmologicalParameters cosmo_params,  SimulationParameters sim_params, StarFormation starformation);

	void write_galaxies(int snapshot, const std::vector<HaloPtr> &halos);

private:

	ExecutionParameters exec_params;
	SimulationParameters sim_params;
	CosmologicalParameters cosmo_params;
	StarFormation starformation;

};

}


#endif /* INCLUDE_WRITE_OUTPUT_H_ */
