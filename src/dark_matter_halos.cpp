/*
 * dark_matter_halos.cpp
 *
 *  Created on: 4Aug.,2017
 *      Author: clagos
 */

#include <cmath>
#include <fstream>
#include <map>
#include <stdexcept>
#include <tuple>

#include "cosmology.h"
#include "dark_matter_halos.h"
#include "logging.h"
#include "numerical_constants.h"


namespace shark {

DarkMatterHaloParameters::DarkMatterHaloParameters(const std::string &filename) :
	Options(filename),
	haloprofile(NFW)
{
	load("dark_matter_halo.halo_profile", haloprofile);
}

namespace detail {

template <>
DarkMatterHaloParameters::DarkMatterProfile Helper<DarkMatterHaloParameters::DarkMatterProfile>::get(const std::string &name, const std::string &value) {
	if ( value == "nfw" ) {
		return DarkMatterHaloParameters::NFW;
	}
	else if ( value == "einasto" ) {
		return DarkMatterHaloParameters::EINASTO;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value << ". Supported values are nfw and einasto";
	throw invalid_option(os.str());
}
}

DarkMatterHalos::DarkMatterHalos(DarkMatterHaloParameters parameters, std::shared_ptr<Cosmology> cosmology) :
	parameters(parameters),
	cosmology(cosmology)
	{
	// no-op
}

double DarkMatterHalos::grav_potential_halo (double r, double c){

	if(r <=0 ){
		return 0;
	}
	else{

		if(parameters.haloprofile == parameters.NFW){

			double x = r/c;

			double g_c = std::log(1+1/c)-1/(1+c);

			double rhos = 1/(constants::PI4 * std::pow(c,3) * g_c);

			double f_x = 0;

			if(r < 1){
				f_x = std::log(1+x)/x-1/(1+1/c);
			}
			else{
				f_x = (std::log(1+1/c)-1/(1+c))/x;
			}
			return -1 * constants::PI4 * rhos * std::pow(c,2) * f_x;

		}
		else if(parameters.haloprofile == parameters.EINASTO){

			//TODO: implement Einasto profile.
			return 0;
		}

	}
}

double DarkMatterHalos::energy_circular (double r, double c){

    return 0.5 * enclosed_mass(r,c)/r+grav_potential_halo(r, c);

}

double DarkMatterHalos::enclosed_mass(double r, double c){

	if(parameters.haloprofile == parameters.NFW){

		double fa = std::log(1+1/c)-1.0/(1.0+c);
		double r_c = r/c;
		if(r < constants::tolerance * c){
			return  fa* (std::pow(r_c,2)) * (0.5+r_c * (-0.6666667+r_c * (0.75-0.8*r_c) ) );
		}
		else{
			 return fa*(std::log(1.0+r/c)-r/(r+c));
		}

	}
	if(parameters.haloprofile == parameters.EINASTO){

		//TODO: implement Einasto profile.
		return 0;
	}
}

double DarkMatterHalos::halo_dynamical_time (std::shared_ptr<Halo> &halo){

	return constants::MPCKM2GYR * halo_virial_radius(halo) / halo->Vvir / cosmology->parameters.Hubble_h;
}

double DarkMatterHalos::halo_virial_radius(std::shared_ptr<Halo> &halo){

	/**
	 * Function to calculate the halo virial radius. Returns virial radius in Mpc/h.
	 */
	return constants::G * halo->Mvir / std::pow(halo->Vvir,2);
}

} // namespace shark

