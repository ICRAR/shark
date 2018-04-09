/*
 * dark_matter_halos.cpp
 *
 *  Created on: 4Aug.,2017
 *      Author: clagos
 */

#include <cmath>
#include <fstream>
#include <map>
#include <random>
#include <stdexcept>
#include <tuple>

#include "cosmology.h"
#include "dark_matter_halos.h"
#include "logging.h"
#include "numerical_constants.h"


namespace shark {

DarkMatterHaloParameters::DarkMatterHaloParameters(const Options &options) :
	haloprofile(NFW),
	random_lambda()
{
	int lambda;
	options.load("dark_matter_halo.halo_profile", haloprofile);
	options.load("dark_matter_halo.lambda_random", lambda);

	if(lambda == 1){
		random_lambda = true;
	}


}

template <>
DarkMatterHaloParameters::DarkMatterProfile
Options::get<DarkMatterHaloParameters::DarkMatterProfile>(const std::string &name, const std::string &value) const {
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

DarkMatterHalos::DarkMatterHalos(DarkMatterHaloParameters &params, std::shared_ptr<Cosmology> cosmology, SimulationParameters &sim_params) :
	params(params),
	cosmology(cosmology),
	sim_params(sim_params),
	generator(),
	distribution(-3.5,-0.69)
	{
	// no-op
}

double DarkMatterHalos::energy_circular (double r, double c){

    return 0.5 * enclosed_mass(r,c)/r+grav_potential_halo(r, c);

}

double DarkMatterHalos::halo_virial_velocity (double mvir, double redshift){

	double hparam = cosmology->hubble_parameter(redshift);

	double vvir = std::cbrt(10.0 *constants::G * mvir * hparam);

	return vvir;
}

double DarkMatterHalos::halo_dynamical_time (HaloPtr &halo, double z){

	auto subhalo_central = halo->central_subhalo;

	double r = halo_virial_radius(*subhalo_central);

	return constants::MPCKM2GYR * cosmology->comoving_to_physical_size(r, z) / subhalo_central->Vvir;
}

double DarkMatterHalos::halo_virial_radius(Subhalo &subhalo){

	/**
	 * Function to calculate the halo virial radius. Returns virial radius in Mpc/h.
	 */
	return constants::G * subhalo.Mvir / std::pow(subhalo.Vvir,2);
}

double DarkMatterHalos::halo_lambda (xyz<float> L, double mvir, double rvir, double z){

	//Spin parameter calculated from j=sqrt(2) * lambda *G^2/3 M^2/3 / (10*H)^1/3.
	//Weak redshift evolution applied according to Elahi et al. (2018, arXiv:1712.01988) in the case of random distribution.

	if(params.random_lambda){
		return distribution(generator) * std::pow(1+z,0.4);
	}

	double E = constants::G * std::pow(mvir,2.0) / (2.0 * rvir);

	double lambda = L.norm() * std::sqrt(E) / constants::G * std::pow(mvir, -2.5);

	if(lambda > 1){
		lambda = 1;
	}

	return lambda;
}

double DarkMatterHalos::disk_size_theory (Subhalo &subhalo, double z){

	//Calculation comes from assuming rdisk = 2/sqrt(2) *lambda *Rvir;
	double Rvir = halo_virial_radius(subhalo);

	double lambda = halo_lambda(subhalo.L, subhalo.Mvir, Rvir, z);

	// Limit value of lambda to 1.
	if(lambda > 1) lambda =1;

	double rdisk = 3/constants::SQRT2 * lambda * Rvir;

	//Numerical factor comes from 1/5 * 1.67. The 1/5 comes from scaling the size to a scale length, and the 1.67 comes from
	//assuming an exponential disk and scaling the scale length to a half mass radius.

	return 0.334 * rdisk;

	//return 0.03 * Rvir;
}

double NFWDarkMatterHalos::grav_potential_halo(double r, double c) const
{

	if (r <=0) {
		return 0;
	}

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

double NFWDarkMatterHalos::enclosed_mass(double r, double c) const
{
	double fa = std::log(1+1/c)-1.0/(1.0+c);
	double r_c = r/c;
	if(r < constants::tolerance * c){
		return  fa* (std::pow(r_c,2)) * (0.5+r_c * (-0.6666667+r_c * (0.75-0.8*r_c) ) );
	}
	return fa*(std::log(1.0+r/c)-r/(r+c));
}

double EinastoDarkMatterHalos::grav_potential_halo(double r, double c) const
{
	//TODO: implement Einasto profile.
	return 0;
}

double EinastoDarkMatterHalos::enclosed_mass(double r, double c) const
{
	//TODO: implement Einasto profile.
	return 0;
}

void DarkMatterHalos::galaxy_velocity(Subhalo &subhalo, Galaxy &galaxy){

	double rvir = halo_virial_radius(subhalo);

	//disk properties.
	double rdisk = galaxy.disk_gas.rscale;
	double mdisk = galaxy.disk_mass();
	double cd = 0;
	if(mdisk > 0 and rdisk > 0){
		cd = rvir / (rdisk / constants::RDISK_HALF_SCALE);
	}

	//bulge properties.
	double rbulge = galaxy.bulge_stars.rscale;
	double mbulge = galaxy.bulge_mass();
	double cb = 0.0;

	if(rbulge > 0){
		cb = rvir/(rbulge / constants::RDISK_HALF_SCALE);
	}

	//halo properties.
	double ch = subhalo.concentration;
	double mvir = subhalo.host_halo->Mvir;

	double xd = rdisk / rvir;

	//Rotational velocity at the half-mass radius of the disk.
	double vd = std::sqrt(v2disk(xd, mdisk, cd, rvir));
	double vb = std::sqrt(v2bulge(xd, mbulge, cb, rvir));
	double vh = std::sqrt(v2halo(xd, mvir, ch, rvir) );

	double v2tot_d = v2halo(xd, mvir, ch, rvir) + v2disk(xd, mdisk, cd, rvir) + v2bulge(xd, mbulge, cb, rvir);

	galaxy.disk_gas.sAM = rdisk * std::sqrt(v2tot_d);

	if(rbulge > 0){

		double xb = rbulge / rvir;
		double v2tot_b = v2halo(xb, mvir, ch, rvir) + v2disk(xb, mdisk, cd, rvir) + v2bulge(xb, mbulge, cb, rvir);

		galaxy.bulge_stars.sAM = rbulge * std::sqrt(v2tot_b);
		galaxy.bulge_gas.sAM   = galaxy.bulge_stars.sAM;
	}

}

double DarkMatterHalos::v2halo (double x, double m, double c, double r){

	double cpx = 1 + c * x;
	double cx = c * x;
	double cp = 1 + c;

	double v = constants::G * m / r * (std::log(cpx) - cx/cpx) / (x * (std::log(cp) - c / cp));

	return v;
}


double DarkMatterHalos::v2disk (double x, double m, double c, double r){

	double cx = c * x;

	double nom = c + 4.8 * c * std::exp(-0.35 * cx - 3.5 /cx);
	double denom = cx + std::pow(cx,-2.0) + 2.0 * std::pow(cx,-0.5);

	double v = constants::G * m / r * nom/denom;

	return v;
}

double DarkMatterHalos::v2bulge (double x, double m, double c, double r){

	if(m <= 0){
		return 0;
	}

	double cx = c * x;

	double nom = std::pow(cx,2.0) * c;
	double denom = std::pow( 1 + std::pow(cx,2.0), 1.5);

	double v = constants::G * m / r * nom/denom;

	return v;

}

double DarkMatterHalos::nfw_concentration(double mvir, double z){

	return 12.3/(1.0+z) * std::pow(mvir/1.3e13,-0.13);

}


} // namespace shark

