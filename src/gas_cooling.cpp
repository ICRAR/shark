/*
 * gas_cooling.cpp
 *
 *  Created on: 17May,2017
 *      Author: clagos
 */

#include <cmath>
#include <fstream>
#include <map>
#include <tuple>

#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "cosmology.h"
#include "gas_cooling.h"
#include "logging.h"
#include "numerical_constants.h"
#include "components.h"

using namespace std;

namespace shark {

GasCoolingParameters::GasCoolingParameters(const std::string &filename) :
	Options(filename),
	rcore(0),
	model(CROTON06),
	lambdamodel(CLOUDY),
	cooling_table()

{
	string cooling_tables_dir;
	load("gas_cooling.model", model);
	load("gas_cooling.rcore", rcore);
	load("gas_cooling.lambdamodel", lambdamodel);
	load("gas_cooling.cooling_tables_dir", cooling_tables_dir, true);

    //read cooling tables and load values in
	string prefix;
	if (lambdamodel == CLOUDY) {
		prefix = "C08.00_";
	}
	else if (lambdamodel == SUTHERLAND) {
		prefix = "S93_";
	}

	string tables = cooling_tables_dir + "/" + prefix + "tables.txt";

	// Collect metallicity tables information
	LOG(debug) << "Reading metallicity table index" << tables;
	string line;
	map<double, string> metallicity_tables;
	ifstream f = open_file(tables);
	while ( getline(f, line) ) {

		trim(line);
		if (is_skipable(line)) {
			continue;
		}

		double metallicity;
		string table_fname;
		istringstream iss(line);
		iss >> metallicity >> table_fname;
		trim(table_fname);

		metallicity_tables[metallicity] = table_fname;
	}
	f.close();

	// Populate cooling_table
	for(auto &kv: metallicity_tables) {

		double metallicity = std::get<0>(kv);
		string fname = std::get<1>(kv);
		fname = cooling_tables_dir + "/" + fname;

		LOG(debug) << "Reading table " << fname << " for metallicity " << metallicity;

		ifstream f = open_file(fname);
		while ( getline(f, line) ) {

			trim(line);
			if (is_skipable(line)) {
				continue;
			}

			double t, ne, nh, nt, logl;
			istringstream iss(line);
			iss >> t >> ne >> nh >> nt >> logl;

			cooling_table.log10lam.push_back(logl);
			cooling_table.log10temp.push_back(t);
			cooling_table.zmetal.push_back(metallicity);

		}
		f.close();
	}
}

namespace detail {

template <>
GasCoolingParameters::LambdaCoolingModel Helper<GasCoolingParameters::LambdaCoolingModel>::get(const std::string &name, const std::string &value) {
	if ( value == "cloudy" ) {
		return GasCoolingParameters::CLOUDY;
	}
	else if ( value == "sutherland" ) {
		return GasCoolingParameters::SUTHERLAND;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value;
	throw invalid_option(os.str());
}


template <>
GasCoolingParameters::CoolingModel Helper<GasCoolingParameters::CoolingModel>::get(const std::string &name, const std::string &value) {
	if ( value == "Croton06" ) {
		return GasCoolingParameters::CROTON06;
	}
	else if ( value == "Galform" ) {
		return GasCoolingParameters::GALFORM;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value;
	throw invalid_option(os.str());
}

} // namespace detail

GasCooling::GasCooling(GasCoolingParameters parameters) :
	parameters(parameters),
	interp(nullptr)
{
	interp.reset(gsl_interp2d_alloc(gsl_interp2d_bilinear, parameters.cooling_table.log10temp.size(), parameters.cooling_table.zmetal.size()));
}

double GasCooling::cooling_rate(std::shared_ptr<Subhalo> &subhalo, double deltat) {


    gsl_interp_accel *xacc = gsl_interp_accel_alloc();
    gsl_interp_accel *yacc = gsl_interp_accel_alloc();

    double coolingrate;

    //TODO: see if here it should be
    double mhot = subhalo->hot_halo_gas.mass;
    double mzhot = subhalo->hot_halo_gas.mass_metals;
    double vvir = subhalo->Mvir;
    double mvir = subhalo->Vvir;

    double zhot = (mzhot/mhot);

    double Tvir = 35.9*std::pow(vvir,2); //in K.
    double lgTvir = log10(Tvir); //in K.

    double Rvir = constants::G*mvir/std::pow(vvir,2); //in Mpc.

    /**
     * This corresponds to the very simple model of Croton06, in which the cooling time is assumed to be equal to the
     * dynamical timescale of the halo. With that assumption, and using an isothermal halo, the calculation of rcool and mcool
     * is trivial.
     */
    if(parameters.model == GasCoolingParameters::CROTON06)
    {

    	double tcoolGyr = Rvir/vvir*constants::KMS2MPCGYR; //in Gyr.

    	double tcool = tcoolGyr*constants::GYR2S; //in seconds.

    	double rho_shell = mhot*constants::MSOLAR_g/constants::PI4/(Rvir*constants::MPC2CM); //in cgs.

    	double logl = gsl_interp2d_eval_extrap(interp.get(), parameters.cooling_table.log10temp.data(), parameters.cooling_table.zmetal.data(), parameters.cooling_table.log10lam.data(), lgTvir, zhot, xacc, yacc); //in cgs

    	double denominator_temp = 1.5*constants::M_Atomic_g*constants::mu_Primordial*constants::k_Boltzmann_erg*Tvir; //in cgs

    	double r_cool = pow(rho_shell*tcool*pow(logl,10)/denominator_temp,0.5)/constants::MPC2CM; //in Mpc.

    	if(r_cool < Rvir){
    		//cooling radius smaller than virial radius
    		coolingrate = 0.5*(r_cool/Rvir)*(mhot/tcoolGyr); //in Msun/Gyr.
    	}
    	else {
    		//cooling radius larger than virial radius
    		coolingrate = mhot/tcoolGyr;
    	}
    }


    /**
     * This corresponds to the GALFORM-like model. Specifically here we are implementing the Benson & Bower (2010) model
     * for cooling, which calculated a time available for cooling depending on the integral of the the temperature
     * mass and cooling time history.
     */

    if(parameters.model == GasCoolingParameters::GALFORM){

    }

    subhalo->cold_halo_gas.mass = coolingrate*deltat;
    subhalo->cold_halo_gas.mass_metals = subhalo->cold_halo_gas.mass/mhot;

	return coolingrate;
}


}  // namespace shark
