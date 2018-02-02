/*
 * reincorporation.cpp
 *
 *  Created on: 22Sep.,2017
 *      Author: clagos
 */

#include <cmath>
#include <memory>

#include "components.h"
#include "dark_matter_halos.h"
#include "reincorporation.h"

namespace shark {

ReincorporationParameters::ReincorporationParameters(const Options &options):
		alpha_reheat(0),
		mhalo_norm(0),
		halo_mass_power(0)
{
	options.load("reincorporation.alpha_reheat",alpha_reheat);
	options.load("reincorporation.mhalo_norm",mhalo_norm);
	options.load("reincorporation.halo_mass_power",halo_mass_power);
}

Reincorporation::Reincorporation(ReincorporationParameters parameters, std::shared_ptr<DarkMatterHalos> darkmatterhalo):
	parameters(parameters),
	darkmatterhalo(darkmatterhalo)
{
	//no-opt
}

double Reincorporation::reincorporated_mass(HaloPtr halo, double z, double delta_t){

	double mvir = halo->Mvir;
	double tdyn = darkmatterhalo->halo_dynamical_time(halo, z);
	double mgas = halo->central_subhalo->ejected_galaxy_gas.mass;

	double treinc = parameters.alpha_reheat * std::pow( (mvir / parameters.mhalo_norm), parameters.halo_mass_power);
	//double reincor_rate = mgas * parameters.alpha_reheat/tdyn * std::pow( (mvir / parameters.mhalo_norm), parameters.halo_mass_power);

	if(treinc == 0){
		return 0;
	}

	return mgas / treinc * delta_t;

}

}
