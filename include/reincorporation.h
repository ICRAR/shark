/*
 * reincorporation.h
 *
 *  Created on: 22Sep.,2017
 *      Author: clagos
 */

#ifndef INCLUDE_REINCORPORATION_H_
#define INCLUDE_REINCORPORATION_H_

#include <memory>
#include <string>
#include <vector>

#include "components.h"
#include "dark_matter_halos.h"

namespace shark {

class ReincorporationParameters{

public:
	ReincorporationParameters(const Options &options);

	double alpha_reheat;
	double mhalo_norm;
	double halo_mass_power;

};

class Reincorporation{

public:
	Reincorporation(ReincorporationParameters parameters, std::shared_ptr<DarkMatterHalos> darkmatterhalo);

	double reincorporated_mass_rate (HaloPtr halo, double z);

private:

	ReincorporationParameters parameters;
	std::shared_ptr<DarkMatterHalos> darkmatterhalo;

};


}//end namespace shark

#endif /* INCLUDE_REINCORPORATION_H_ */
