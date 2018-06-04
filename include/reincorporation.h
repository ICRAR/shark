/*
 * reincorporation.h
 *
 *  Created on: 22Sep.,2017
 *      Author: clagos
 */

#ifndef INCLUDE_REINCORPORATION_H_
#define INCLUDE_REINCORPORATION_H_

#include <memory>
#include <utility>

#include "components.h"
#include "dark_matter_halos.h"
#include "options.h"

namespace shark {

class ReincorporationParameters{

public:
	ReincorporationParameters(const Options &options);

	double alpha_reheat = 0;
	double mhalo_norm = 0;
	double halo_mass_power = 0;

};

class Reincorporation{

public:
	Reincorporation(const ReincorporationParameters &parameters, const DarkMatterHalosPtr &darkmatterhalo);

	double reincorporated_mass (HaloPtr halo, double z, double delta_t);

private:

	ReincorporationParameters parameters;
	DarkMatterHalosPtr darkmatterhalo;
};

typedef std::shared_ptr<Reincorporation> ReincorporationPtr;

template <typename ...Ts>
ReincorporationPtr make_reincorporation(Ts&&...ts)
{
	return std::make_shared<Reincorporation>(std::forward<Ts>(ts)...);
}

}//end namespace shark

#endif /* INCLUDE_REINCORPORATION_H_ */
