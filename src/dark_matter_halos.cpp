//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2017
// Copyright by UWA (in the framework of the ICRAR)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//

/**
 * @file
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
#include "nfw_distribution.h"
#include "numerical_constants.h"
#include "utils.h"


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
	auto lvalue = lower(value);
	if (lvalue == "nfw") {
		return DarkMatterHaloParameters::NFW;
	}
	else if (lvalue == "einasto") {
		return DarkMatterHaloParameters::EINASTO;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value << ". Supported values are nfw and einasto";
	throw invalid_option(os.str());
}


template <>
DarkMatterHaloParameters::SizeModel
Options::get<DarkMatterHaloParameters::SizeModel>(const std::string &name, const std::string &value) const {
	auto lvalue = lower(value);
	if (lvalue == "mo98") {
		return DarkMatterHaloParameters::MO98;
	}
	else if (lvalue == "cole00") {
		return DarkMatterHaloParameters::COLE00;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value << ". Supported values are Mo98 and Cole00";
	throw invalid_option(os.str());
}

DarkMatterHalos::DarkMatterHalos(const DarkMatterHaloParameters &params, const CosmologyPtr &cosmology, SimulationParameters &sim_params) :
	params(params),
	cosmology(cosmology),
	sim_params(sim_params),
	generator(),
	distribution(-3.5,-0.69),
	flat_distribution(0,1)
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

double DarkMatterHalos::halo_dynamical_time (HaloPtr &halo, double z)
{
	return subhalo_dynamical_time(*halo->central_subhalo, z);
}

double DarkMatterHalos::subhalo_dynamical_time (Subhalo &subhalo, double z){

	double r = halo_virial_radius(subhalo);

	return constants::MPCKM2GYR * cosmology->comoving_to_physical_size(r, z) / subhalo.Vvir;
}

double DarkMatterHalos::halo_virial_radius(Subhalo &subhalo){

	/**
	 * Function to calculate the halo virial radius. Returns virial radius in Mpc/h.
	 */
	return constants::G * subhalo.Mvir / std::pow(subhalo.Vvir,2);
}

float DarkMatterHalos::halo_lambda (float lambda, double z){

	//Spin parameter either read from the DM files or assumed a random distribution.

	if(params.random_lambda){
		return distribution(generator);
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
		return 0;
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

	if (std::isnan(galaxy.disk_gas.sAM) or std::isnan(galaxy.disk_gas.rscale)) {
		throw invalid_argument("rgas or sAM are NaN, cannot continue at disk_sAM in dark_matter_halos");
	}

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

	// From Duffy et al. (2008). Full sample from z=0-2 for Virial masses.

	return 7.85 * std::pow(1.0+z, -0.71) * std::pow(mvir/2.0e12,-0.081);

}

/// Specialization of lambert_w0 implemented using GSL
template <>
struct lambert_w0<double>
{
	double operator()(double x) {
		return gsl_sf_lambert_W0(x);
	}
};

void DarkMatterHalos::generate_random_orbits(xyz<float> &pos, xyz<float> &v, xyz<float> &L, double total_am, const HaloPtr &halo){

	double c = halo->concentration;

	auto   subhalo = halo->central_subhalo;
	double rvir = halo_virial_radius(*subhalo);

	auto   pos_halo = subhalo->position;

	// Assign positions based on an NFW halo of concentration c.
	nfw_distribution<double> r(c);
	double rproj = r(generator);
	double theta = std::acos(flat_distribution(generator)*2.0 - 1); //flat between -1 and 1.
	double phi   = flat_distribution(generator)*constants::PI2; //flat between 0 and 2PI.

	pos.x = rproj * std::cos(theta) * std::cos(phi) * rvir + pos_halo.x;
	pos.y = rproj * std::cos(theta) * std::sin(phi) * rvir + pos_halo.y;
	pos.z = rproj * std::sin(theta) * rvir + pos_halo.z;

	// Assign velocities using NFW velocity dispersion and assuming isotropy.
	// f_c equation from Manera et al. (2013; eq. 23).
	double f_c = c * (0.5 * c / (c + 1) - std::log(1 + c) / (1 + c))/ std::pow(std::log(1 + c) - c / (1 + c), 2.0);

	if(f_c < 0){
		//Negative values happen if c <<~ 2.6, and thus we set the minimum value to f_c evaluated in c = 2.6.
		f_c = 0.04411218227;
	}

	// 1D velocity dispersion.
	double sigma = std::pow(0.333  * constants::G * halo->Mvir / rvir * f_c, 0.5);

	std::normal_distribution<double> normal_distribution(0,sigma);

	v.x = normal_distribution(generator) + halo->velocity.x;
	v.y = normal_distribution(generator) + halo->velocity.y;
	v.z = normal_distribution(generator) + halo->velocity.z;

	// Assign angular momentum based on random angles,
	// drawn random angles again
	theta = std::acos(flat_distribution(generator)*2.0 - 1); //flat between -1 and 1.
	phi   = flat_distribution(generator)*constants::PI2; //flat between 0 and 2PI.
	L.x   = total_am *  std::sin(theta) * std::cos(phi);
	L.y   = total_am *  std::sin(theta) * std::sin(phi);
	L.z   = total_am *  std::cos(theta);


}





} // namespace shark

