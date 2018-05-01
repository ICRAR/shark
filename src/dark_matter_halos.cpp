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

DarkMatterHaloParameters::DarkMatterHaloParameters(const Options &options)
{
	options.load("dark_matter_halo.halo_profile", haloprofile);
	options.load("dark_matter_halo.size_model", sizemodel);
	options.load("dark_matter_halo.lambda_random", random_lambda);

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


template <>
DarkMatterHaloParameters::SizeModel
Options::get<DarkMatterHaloParameters::SizeModel>(const std::string &name, const std::string &value) const {
	if ( value == "Mo98" ) {
		return DarkMatterHaloParameters::MO98;
	}
	else if ( value == "Cole00" ) {
		return DarkMatterHaloParameters::COLE00;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value << ". Supported values are Mo98 and Cole00";
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

float DarkMatterHalos::halo_lambda (float lambda, double z){

	//Spin parameter either read from the DM files or assumed a random distribution.
	//In the case of the random distribution, a weak redshift evolution is
	//applied according to Elahi et al. (2018, arXiv:1712.01988).

	if(params.random_lambda){
		return distribution(generator) * std::pow(1+z,0.4);
	}
	else{
		//take the value read from the DM merger trees, but limit it to a maximum of 1.
		if(lambda > 1){
			lambda = 1;
		}
		return lambda;
	}

}

double DarkMatterHalos::disk_size_theory (Subhalo &subhalo, double z){

	if(params.sizemodel == DarkMatterHaloParameters::MO98){
		//Calculation comes from assuming rdisk = 2/sqrt(2) *lambda *Rvir;
		double Rvir = halo_virial_radius(subhalo);

		double lambda = subhalo.lambda;

		double rdisk = 3/constants::SQRT2 * lambda * Rvir;

		//Numerical factor comes from 1/5 * 1.67. The 1/5 comes from scaling the size to a scale length, and the 1.67 comes from
		//assuming an exponential disk and scaling the scale length to a half mass radius.
		return 0.334 * rdisk;

	}
	else if (params.sizemodel == DarkMatterHaloParameters::COLE00){
		//TODO
	}
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

void DarkMatterHalos::cooling_gas_sAM(Subhalo &subhalo, double z){

	if(params.random_lambda){
		double H0 = 10.0* cosmology->hubble_parameter(z);
		subhalo.cold_halo_gas.sAM = constants::SQRT2 * std::pow(constants::G,0.66) *
								    subhalo.lambda * std::pow(subhalo.Mvir,0.66) / std::pow(H0,0.33);
	}
	else{
		subhalo.cold_halo_gas.sAM = subhalo.hot_halo_gas.sAM;
	}
	/*if(params.sizemodel == DarkMatterHaloParameters::MO98){
		//Assumes that cooled gas brings a specific angular momentum that is equivalent to
		//the disk specific angular momentum that is set by the disk_gas.rscale and disk_gas.sAM.
		subhalo.cold_halo_gas.sAM = subhalo.hot_halo_gas.sAM; //Mo98_j; //galaxy.disk_gas.sAM;
	}
	else if (params.sizemodel == DarkMatterHaloParameters::COLE00){
		//TODO
	}*/

}

void DarkMatterHalos::disk_sAM(Subhalo &subhalo, Galaxy &galaxy){

	double rvir = halo_virial_radius(subhalo);

	//disk properties. Use composite size of disk.
	double rdisk = galaxy.disk_size();
	double mdisk = galaxy.disk_mass();
	double cd = 0;

	if(rdisk > 0){
		cd = rvir / (rdisk / constants::RDISK_HALF_SCALE);
	}

	//bulge properties. Use composite size of bulge.
	double rbulge = galaxy.bulge_size();
	double mbulge = galaxy.bulge_mass();
	double cb = 0.0;

	if(rbulge > 0){
		cb = rvir/(rbulge);
	}

	//halo properties.
	double ch = subhalo.concentration;
	double mvir = subhalo.host_halo->Mvir;

	double xd = rdisk / rvir;

	//Rotational velocity at the half-mass radius of the disk.
	double v2tot_d = v2halo(xd, mvir, ch, rvir) + v2disk(xd, mdisk, cd, rvir) + v2bulge(xd, mbulge, cb, rvir);

	galaxy.disk_gas.sAM = 2.0 * rdisk / constants::RDISK_HALF_SCALE * std::sqrt(v2tot_d);

}

void DarkMatterHalos::bulge_sAM(Subhalo &subhalo, Galaxy &galaxy){

	double rvir = halo_virial_radius(subhalo);

	//disk properties. Use composite size of disk.
	double rdisk = galaxy.disk_size();
	double mdisk = galaxy.disk_mass();
	double cd = 0;

	if(rdisk > 0){
		cd = rvir / (rdisk / constants::RDISK_HALF_SCALE);
	}

	//bulge properties. Use composite size of bulge.
	double rbulge = galaxy.bulge_size();
	double mbulge = galaxy.bulge_mass();
	double cb = 0.0;

	if(rbulge > 0){
		cb = rvir/(rbulge);
	}

	//halo properties.
	double ch = subhalo.concentration;
	double mvir = subhalo.host_halo->Mvir;

	if(rbulge > 0){
		double xb = rbulge / rvir;
		double v2tot_b = v2halo(xb, mvir, ch, rvir) + v2disk(xb, mdisk, cd, rvir) + v2bulge(xb, mbulge, cb, rvir);
		double vbulge = std::sqrt(v2tot_b);

		if(vbulge > 1e3){
			double vb = std::sqrt(constants::G * mbulge / rbulge);
			double medd=0;
		}

		galaxy.bulge_gas.sAM = 2.0 * rbulge / constants::RDISK_HALF_SCALE * vbulge;
	}

}

void DarkMatterHalos::transfer_bulge_am(SubhaloPtr &subhalo, GalaxyPtr &galaxy, double z){

	//modify AM based on mass weighting. How to do this should depend on size model.
	if(params.sizemodel == DarkMatterHaloParameters::MO98){
  	   	galaxy->disk_gas.rscale = disk_size_theory(*subhalo, z);

  	   	//define disk angular momentum.
  	   	disk_sAM(*subhalo, *galaxy);
	}
	else if (params.sizemodel == DarkMatterHaloParameters::COLE00){
		//TODO
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

