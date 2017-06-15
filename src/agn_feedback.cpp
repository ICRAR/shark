/*
 * agn_feedback.cpp
 *
 *  Created on: 15Jun.,2017
 *      Author: clagos
 */

#include <cmath>
#include <memory>

#include "numerical_constants.h"
#include "agn_feedback.h"

using namespace std;

namespace shark {

AGNFeedbackParameters::AGNFeedbackParameters(const std::string &filename) :
	Options(filename),
	mseed(0),
	alpha_cool(0),
	epsilon_smbh(0),
	accretion_eff_cooling(0),
	accretion_eff_bursts(0),
	mhalo_seed(0)
{
	load("agn_feedback.mseed",mseed);
	load("agn_feedback.alpha_cool",alpha_cool);
	load("agn_feedback.epsilon_smbh",epsilon_smbh);
	load("agn_feedback.accretion_eff_cooling",accretion_eff_cooling);
	load("agn_feedback.accretion_eff_bursts",accretion_eff_bursts);
	load("agn_feedback.mhalo_seed",mhalo_seed);

}


AGNFeedback::AGNFeedback(AGNFeedbackParameters parameters, std::shared_ptr<Cosmology> cosmology) :
	parameters(parameters),
	cosmology(cosmology)
{
	// no-op
}

void AGNFeedback::plant_seed_smbh(Subhalo &subhalo){

	if(subhalo.Mvir > parameters.mhalo_seed){
		for(shared_ptr<Galaxy> &galaxies: subhalo.galaxies) {
			if(galaxies->galaxy_type == Galaxy::CENTRAL && galaxies->smbh.mass == 0){
				galaxies->smbh.mass = 0;
			}
		}
	}
}

double AGNFeedback::eddington_luminosity(double mbh){
	// Numerical constants appearing in the expression for the Eddington luminosity (4 pi c G M_sun M_H/sigma_T)

    //Eddington luminosity in units of $10^{40} erg/s$
    return constants::Eddngtn_Lmnsty_Scale_Factor*mbh/constants::ERG2J;

}

}
