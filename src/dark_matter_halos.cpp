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

#include <gsl/gsl_sf_lambert.h>

#include "cosmology.h"
#include "dark_matter_halos.h"
#include "galaxy.h"
#include "halo.h"
#include "logging.h"
#include "nfw_distribution.h"
#include "numerical_constants.h"
#include "subhalo.h"
#include "utils.h"


namespace shark {

DarkMatterHaloParameters::DarkMatterHaloParameters(const Options &options)
{
	options.load("dark_matter_halo.halo_profile", haloprofile);
	options.load("dark_matter_halo.size_model", sizemodel);
	options.load("dark_matter_halo.lambda_random", random_lambda);
	options.load("dark_matter_halo.concentration_model", concentrationmodel);
	options.load("dark_matter_halo.use_converged_lambda_catalog", use_converged_lambda_catalog);
	options.load("dark_matter_halo.min_part_convergence", min_part_convergence);

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

template <>
DarkMatterHaloParameters::ConcentrationModel
Options::get<DarkMatterHaloParameters::ConcentrationModel>(const std::string &name, const std::string &value) const {
	auto lvalue = lower(value);
	if (lvalue == "duffy08") {
		return DarkMatterHaloParameters::DUFFY08;
	}
	else if (lvalue == "dutton14") {
		return DarkMatterHaloParameters::DUTTON14;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value << ". Supported values are Duffy08 and Dutton14";
	throw invalid_option(os.str());
}


DarkMatterHalos::DarkMatterHalos(
	const DarkMatterHaloParameters &params,
	CosmologyPtr cosmology,
	SimulationParameters &sim_params,
	ExecutionParameters exec_params) :
	params(params),
	cosmology(std::move(cosmology)),
	sim_params(sim_params),
	exec_params(std::move(exec_params))
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
	double v = halo_virial_velocity(halo->Mvir, z);
	double r = constants::G * halo->Mvir / std::pow(v,2);
	double t = constants::MPCKM2GYR * cosmology->comoving_to_physical_size(r, z) / v;
	return t;
}


double DarkMatterHalos::subhalo_dynamical_time (Subhalo &subhalo, double z){

	double r = 0;
	double v = 0;

	if(subhalo.subhalo_type == Subhalo::CENTRAL){
		r = halo_virial_radius(subhalo.host_halo, z);
		v = subhalo.Vvir;
	}
	else{
		r = subhalo.rvir_infall;
		v = halo_virial_velocity(subhalo.Mvir_infall, subhalo.infall_t);
	}

	return constants::MPCKM2GYR * cosmology->comoving_to_physical_size(r, z) / v;
}

double DarkMatterHalos::halo_virial_radius(const HaloPtr &halo, double z){

	/**
	 * Function to calculate the halo virial radius. Returns virial radius in physical Mpc/h.
	 */
        double v = halo_virial_velocity(halo->Mvir, z);
	return constants::G * halo->Mvir / std::pow(v,2);
}


float DarkMatterHalos::halo_lambda (const Subhalo &subhalo, float m, double z, double npart){

	//Spin parameter either read from the DM files or assumed a random distribution.
	double H0 = cosmology->hubble_parameter(z);
	double lambda = subhalo.L.norm() / m * 1.5234153 / std::pow(constants::G * m, 0.666) * std::pow(H0,0.33);

	if(lambda > 1){
			lambda = 1;
	}

	// Prime the generator with a known seed to allow for reproducible runs
	// using a very weak dependence on Mhalo for the spin distribution, following Kim et al. (2015): arxiv:1508.06037
	double lambda_cen_mhalo = 0.00895651600584195 * std::log10(m)  - 0.07580254755439589;
	std::default_random_engine generator(exec_params.get_seed(subhalo));
	std::lognormal_distribution<double> distribution(std::log(lambda_cen_mhalo), std::abs(std::log(0.5)));
	auto lambda_random = distribution(generator);

	// Avoid zero values. In that case assume small lambda value.
	if(lambda_random == 0){
		lambda_random = 1e-3;
	}
	else if (lambda_random > 1){
		 lambda_random = 1;
	}

	if(params.random_lambda && !params.use_converged_lambda_catalog){
		return lambda_random;
	}
	else if (params.random_lambda && params.use_converged_lambda_catalog && npart < params.min_part_convergence){
		return lambda_random;
	} 
	else {
		//take the value read from the DM merger trees, that has been limited to a maximum of 1.
		return lambda;
	}
 
}

double DarkMatterHalos::disk_size_theory (Subhalo &subhalo, double z){

	if(params.sizemodel == DarkMatterHaloParameters::MO98){
		//Calculation comes from assuming rdisk = 2/sqrt(2) *lambda *Rvir;
		double Rvir = halo_virial_radius(subhalo.host_halo, z);

		double lambda = subhalo.lambda;

		double rdisk = 3/constants::SQRT2 * lambda * Rvir;

		//Numerical factor comes from 1/5 * 1.67. The 1/5 comes from scaling the size to a scale length, and the 1.67 comes from
		//assuming an exponential disk and scaling the scale length to a half mass radius.
		return 0.334 * rdisk;

	}
	throw std::runtime_error("Only the MO98 disk size model is implemented");
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
	// r is normalized by the virial radius.

	double nom = 1.0 / (1.0 + c * r) - 1.0 + std::log(1.0 + c*r);
	double denom = 1.0 / (1.0 + c) - 1.0 + std::log(1.0 + c);

	double frac = nom/denom;

	if(frac > 1){
		frac =1;
	}
	return frac ;

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

float DarkMatterHalos::enclosed_total_mass(const Subhalo &subhalo, double z, float r){

	ConstGalaxyPtr galaxy;
	if(subhalo.subhalo_type == Subhalo::SATELLITE){
		galaxy = subhalo.type1_galaxy();
	}
	else{
		galaxy = subhalo.central_galaxy();
	}
	double mgal = 0;

	auto rvir = halo_virial_radius(subhalo.host_halo, z);
	auto rnorm = r/rvir;

	//calculate enclosed DM mass
	auto mdm = subhalo.Mvir * enclosed_mass(rvir, subhalo.concentration);

	//calculate enclosed hot gas mass (only relevant for isothermal sphere)
	auto mhot = subhalo.hot_halo_gas.mass * std::pow(rnorm,2);

	//calculate enclosed galaxy mass
	if(galaxy){
		mgal = galaxy->enclosed_mass_disk(r) + galaxy->enclosed_bulge_mass(r);
	}

	return mdm + mhot + mgal;
}

void DarkMatterHalos::disk_sAM(Subhalo &subhalo, Galaxy &galaxy, double z){

	double rvir = halo_virial_radius(subhalo.host_halo, z);

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

	if (std::isnan(galaxy.disk_gas.sAM) || std::isnan(galaxy.disk_gas.rscale)) {
		throw invalid_argument("rgas or sAM are NaN, cannot continue at disk_sAM in dark_matter_halos");
	}

}

void DarkMatterHalos::bulge_sAM(Subhalo &subhalo, Galaxy &galaxy, double z){

	double rvir = halo_virial_radius(subhalo.host_halo, z);

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

void DarkMatterHalos::transfer_bulge_am(SubhaloPtr &subhalo, Galaxy &galaxy, double z){

	//modify AM based on mass weighting. How to do this should depend on size model.
	if(params.sizemodel == DarkMatterHaloParameters::MO98){
		galaxy.disk_gas.rscale = disk_size_theory(*subhalo, z);

		// define disk angular momentum.
		disk_sAM(*subhalo, galaxy, z);
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

	if(params.concentrationmodel == DarkMatterHaloParameters::DUFFY08){ 
		// From Duffy et al. (2008). Full sample from z=0-2 for Virial masses.
		return 7.85 * std::pow(1.0+z, -0.71) * std::pow(mvir/2.0e12,-0.081);
	}
	else if(params.concentrationmodel == DarkMatterHaloParameters::DUTTON14){
		//From Dutton & Maccio (2014) for virial masses.
		double b = -0.097 + 0.024 * z;
		double a = 0.537 + (1.025 - 0.537) * std::exp(-0.718 * std::pow(z,1.08));
		return std::pow(10.0, a + b * std::log10(mvir/1e12));
	}
	throw std::runtime_error("Only the Duffy08 and Dutton14 concentration models are implemented");
}

/// Specialization of lambert_w0 implemented using GSL
template <>
struct lambert_w0<double>
{
	double operator()(double x) {
		return gsl_sf_lambert_W0(x);
	}
};

xyz<float> DarkMatterHalos::random_point_in_sphere(float r, std::default_random_engine &generator)
{
	// We distribute cos_theta flatly instead of theta itself to end up with a
	// more uniform distribution of points in the sphere
	std::uniform_real_distribution<float> flat_distribution(0, 1);
	float cos_theta = flat_distribution(generator) * 2.0 - 1; //flat between -1 and 1.
	float theta = std::acos(cos_theta);
	float sin_theta = std::sin(theta);
	float phi = flat_distribution(generator) * constants::PI2; //flat between 0 and 2PI.
	return {
		sin_theta * std::cos(phi) * r,
		sin_theta * std::sin(phi) * r,
		cos_theta * r
	};
}

void DarkMatterHalos::generate_random_orbits(xyz<float> &pos, xyz<float> &v, xyz<float> &L, double total_am, const HaloPtr &halo, const Galaxy &galaxy)
{

	// Prime the generator with a known seed to allow for reproducible runs
	std::default_random_engine generator(exec_params.get_seed(galaxy));

	double c = halo->concentration;

	double rvir = constants::G * halo->Mvir / std::pow(halo->Vvir,2);

	// Assign positions based on an NFW halo of concentration c.
	nfw_distribution<double> r(c);
	double rproj = r(generator);
	pos = halo->position + random_point_in_sphere(rvir * rproj, generator);

	// Assign velocities using NFW velocity dispersion at the radius in which the galaxy is and assuming isotropy.
	double sigma = std::sqrt(0.333 * constants::G * halo->Mvir * enclosed_mass(rproj, c) / (rvir * rproj));

	// Add the velocity bias of Dieman et al. (2005) found between subhalos and DM particles.
	sigma = sigma * 1.12 * std::pow(rproj, -0.1);

	std::normal_distribution<double> normal_distribution(0, sigma);
	xyz<double> delta_v {normal_distribution(generator), normal_distribution(generator), normal_distribution(generator)};

	//delta_v and velocity are in physical km/s.
	v = halo->velocity + delta_v;

	// Assign angular momentum based on random angles,
	L = random_point_in_sphere(total_am, generator);

}

} // namespace shark

