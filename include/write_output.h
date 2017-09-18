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

#include "execution.h"
namespace shark {

class WriteOutput{

public:

	WriteOutput(ExecutionParameters exec_params);

private:

	ExecutionParameters exec_params;

	void write_galaxies(int snapshot, std::vector<HaloPtr> halos);
};

}


#endif /* INCLUDE_WRITE_OUTPUT_H_ */
