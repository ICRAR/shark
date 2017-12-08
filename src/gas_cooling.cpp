/*
 * gas_cooling.cpp
 *
 *  Created on: 17May,2017
 *      Author: clagos
 */

#include <cmath>
#include <fstream>
#include <iterator>
#include <map>
#include <numeric>
#include <tuple>

#include "logging.h"
#include "components.h"
#include "cosmology.h"
#include "gas_cooling.h"
#include "numerical_constants.h"
#include "reincorporation.h"

using namespace std;

namespace shark {

void CoolingTable::add_metallicity_measurements(double zmetal, const std::map<double, double> &records)
{
	// Check that the keys (i.e., the list of temperatures) on the incoming map
	// are the same than those from other measurements.
	// In other words, make sure that for all metallicities we have measurements
	// at the same temperatures, which gives us a nicely populated grid
	// over which we can interpolate later
	if (!_table.empty()) {

		auto first_metallicity = std::begin(_table)->first;
		auto &first_measurement = std::begin(_table)->second;
		vector<double> first_keys = get_keys(first_measurement);
		vector<double> these_keys = get_keys(first_measurement);

		if (first_keys != these_keys) {
			std::ostringstream os;
			os << "Temperature measurements for metallicity " << first_metallicity;
			os << " are different from those measured for metallicity " << zmetal;
			throw invalid_data(os.str());
		}
	}

	_table[zmetal] = records;
}

std::vector<double> CoolingTable::get_metallicities()
{
	return get_keys(_table);
}

std::vector<double> CoolingTable::get_temperatures()
{
	return get_keys(std::begin(_table)->second);
}

std::vector<double> CoolingTable::get_lambda()
{
	std::vector<double> lambda_values;
	for(auto &metallicity_measurement: _table) {
		for(auto &t_record: metallicity_measurement.second) {
			lambda_values.push_back(t_record.second);
		}
	}
	return lambda_values;
}

GasCoolingParameters::GasCoolingParameters(const Options &options) :
	rcore(0),
	pre_enrich_z(1e-4),
	lambdamodel(CLOUDY),
	model(CROTON06),
	cooling_table()
{
	string cooling_tables_dir;
	options.load("gas_cooling.model", model, true);
	options.load("gas_cooling.rcore", rcore);
	options.load("gas_cooling.lambdamodel", lambdamodel, true);
	options.load("gas_cooling.cooling_tables_dir", cooling_tables_dir, true);
	options.load("gas_cooling.pre_enrich_z", pre_enrich_z);

	tables_idx metallicity_tables = find_tables(cooling_tables_dir);
	load_tables(cooling_tables_dir, metallicity_tables);

}

GasCoolingParameters::tables_idx GasCoolingParameters::find_tables(
	const string &cooling_tables_dir)
{
	//read cooling tables and load values in
	string prefix;
	if (lambdamodel == CLOUDY) {
		prefix = "C08.00_";
	}
	else if (lambdamodel == SUTHERLAND) {
		prefix = "S93_";
	}
	else {
		throw invalid_argument("Cooling model is not valid");
	}

	string tables = cooling_tables_dir + "/" + prefix + "tables.txt";

	// Collect metallicity tables information
	LOG(debug) << "Reading metallicity table index" << tables;
	string line;
	map<double, string> metallicity_tables;
	ifstream f = open_file(tables);
	while ( getline(f, line) ) {

		trim(line);
		if (empty_or_comment(line)) {
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

	return metallicity_tables;
}

void GasCoolingParameters::load_tables(
	const string &cooling_tables_dir,
	const tables_idx &metallicity_tables)
{

	// Populate cooling_table
	for(auto &kv: metallicity_tables) {

		double metallicity = std::get<0>(kv);
		const string fname = cooling_tables_dir + "/" + std::get<1>(kv);

		LOG(debug) << "Reading table " << fname << " for metallicity " << metallicity;

		ifstream f = open_file(fname);
		string line;

		std::map<double, double> measurements;
		while ( getline(f, line) ) {

			trim(line);
			if (empty_or_comment(line)) {
				continue;
			}

			double t, ne, nh, nt, logl;
			istringstream iss(line);
			iss >> t >> ne >> nh >> nt >> logl;

			measurements[t] = logl;
		}
		f.close();

		cooling_table.add_metallicity_measurements(metallicity, measurements);
	}

}

template <>
GasCoolingParameters::LambdaCoolingModel
Options::get<GasCoolingParameters::LambdaCoolingModel>(const std::string &name, const std::string &value) const {
	if ( value == "cloudy" ) {
		return GasCoolingParameters::CLOUDY;
	}
	else if ( value == "sutherland" ) {
		return GasCoolingParameters::SUTHERLAND;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value << ". Supported values are cloudy and sutherland";
	throw invalid_option(os.str());
}


template <>
GasCoolingParameters::CoolingModel
Options::get<GasCoolingParameters::CoolingModel>(const std::string &name, const std::string &value) const {
	if ( value == "Croton06" ) {
		return GasCoolingParameters::CROTON06;
	}
	else if ( value == "Benson10" ) {
		return GasCoolingParameters::BENSON10;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value << ". Supported values are Croton06 and Galform";
	throw invalid_option(os.str());
}

GasCooling::GasCooling(GasCoolingParameters parameters,
		ReionisationParameters reio_parameters,
		std::shared_ptr<Cosmology> cosmology,
		std::shared_ptr<AGNFeedback> agnfeedback,
		std::shared_ptr<DarkMatterHalos> darkmatterhalos,
		std::shared_ptr<Reincorporation> reincorporation) :
	reio_parameters(reio_parameters),
	parameters(parameters),
	cosmology(cosmology),
	agnfeedback(agnfeedback),
	darkmatterhalos(darkmatterhalos),
	reincorporation(reincorporation),
	cooling_lambda_interpolator(parameters.cooling_table.get_temperatures(), parameters.cooling_table.get_metallicities(), parameters.cooling_table.get_lambda())
{
	//no-opt
}

double GasCooling::cooling_rate(Subhalo &subhalo, Galaxy &galaxy, double z, double deltat) {

	using namespace constants;

    double coolingrate = 0;

    //Define host halo
    auto halo = subhalo.host_halo;

    /**
     * For now assume that gas can cool only in central subhalos and to central galaxies.
     */

    if(subhalo.subhalo_type == Subhalo::SATELLITE){
    	//Assume satellite subhalos have 0 cooling rate, and give any hot mass to the central subhalo.
    	//TODO: add gradual ram pressure stripping.

    	if(subhalo.hot_halo_gas.mass > 0){
    		subhalo.host_halo->central_subhalo->hot_halo_gas.mass += subhalo.hot_halo_gas.mass;
    		subhalo.host_halo->central_subhalo->hot_halo_gas.mass_metals += subhalo.hot_halo_gas.mass_metals;

    		subhalo.hot_halo_gas.mass = 0;
    		subhalo.hot_halo_gas.mass_metals = 0;
    	}

    	if(subhalo.ejected_galaxy_gas.mass > 0){
    		subhalo.host_halo->central_subhalo->ejected_galaxy_gas.mass += subhalo.ejected_galaxy_gas.mass;
    		subhalo.host_halo->central_subhalo->ejected_galaxy_gas.mass_metals += subhalo.ejected_galaxy_gas.mass_metals;

        	subhalo.host_halo->central_subhalo->hot_halo_gas.mass += subhalo.hot_halo_gas.mass;
        	subhalo.host_halo->central_subhalo->hot_halo_gas.mass_metals += subhalo.hot_halo_gas.mass_metals;
    	}

    	subhalo.cold_halo_gas.mass = 0;
    	subhalo.cold_halo_gas.mass_metals = 0;

    	return 0;
    }

    if(galaxy.galaxy_type != Galaxy::CENTRAL){
    	return 0;
    }

   	/**
   	 * Estimate disk size and specific angular momentum.
   	 */
   	auto central_galaxy = subhalo.central_galaxy();
   	central_galaxy->disk_gas.rscale = darkmatterhalos->disk_size_theory(subhalo);
   	central_galaxy->disk_stars.rscale = central_galaxy->disk_gas.rscale;

   	darkmatterhalos->galaxy_velocity(subhalo);
   	//TODO: remove this part and calculate rscale and sAM of the stellar disk properly.
   	central_galaxy->disk_stars.sAM = central_galaxy->disk_gas.sAM;

    /**
     * Plant black hole seed if necessary.
     */
    agnfeedback->plant_seed_smbh(*halo);

    // Calculate Eddington luminosity of BH in central galaxy.

    double Ledd = agnfeedback->eddington_luminosity(central_galaxy->smbh.mass);

    /**
     * Test for subhalos that are affected by reionisation
     */
    if(subhalo.Vvir < reio_parameters.vcut && z < reio_parameters.zcut){
    	return 0;
    }


    // Calculate reincorporated mass and metals.
	double mreinc_mass = reincorporation->reincorporated_mass_rate(halo) * deltat;

   	if(mreinc_mass > subhalo.ejected_galaxy_gas.mass){
   		mreinc_mass = subhalo.ejected_galaxy_gas.mass;
   	}

   	double mreinc_mass_metals = 0.0;
   	// If reincorporation mass is >0 then modify gas budget.
   	if(mreinc_mass > 0){
   		mreinc_mass_metals = mreinc_mass/subhalo.ejected_galaxy_gas.mass * subhalo.ejected_galaxy_gas.mass_metals;
      	subhalo.ejected_galaxy_gas.mass -= mreinc_mass;
      	subhalo.ejected_galaxy_gas.mass_metals -= mreinc_mass_metals;
      	subhalo.hot_halo_gas.mass += mreinc_mass;
      	subhalo.hot_halo_gas.mass_metals += mreinc_mass_metals;
   	}

  	// Avoid negative values.
  	if(subhalo.ejected_galaxy_gas.mass < constants::tolerance){
   		subhalo.ejected_galaxy_gas.mass = 0;
   		subhalo.ejected_galaxy_gas.mass_metals =0;
   	}
   	if(subhalo.ejected_galaxy_gas.mass_metals < 0){
   		subhalo.ejected_galaxy_gas.mass_metals = 0;
   	}
    	// Add up accreted mass and metals.
   	subhalo.hot_halo_gas.mass += subhalo.accreted_mass;
   	subhalo.hot_halo_gas.mass_metals += subhalo.accreted_mass * parameters.pre_enrich_z;

   	/**
   	 * We need to convert masses and velocities to physical units before proceeding with calculation.
   	 */
   	double mhot = cosmology->comoving_to_physical_mass(subhalo.hot_halo_gas.mass+subhalo.cold_halo_gas.mass);
   	double mzhot = cosmology->comoving_to_physical_mass(subhalo.hot_halo_gas.mass_metals+subhalo.cold_halo_gas.mass_metals);

   	double vvir = halo->Vvir;

 	double zhot = 0;
 	if(mhot > 0){
 		zhot = (mzhot/mhot);
 	}

 	if(mhot < 0 or mhot >1e17 or std::isnan(mhot)){
		std::ostringstream os;
		os << halo << " has hot halo gas mass not well defined";
		throw invalid_data(os.str());
 	}

   	double Tvir = 35.9*std::pow(vvir,2.0); //in K.
   	double lgTvir = log10(Tvir); //in K.
	double Rvir = darkmatterhalos->halo_virial_radius(halo)/cosmology->parameters.Hubble_h;//Mpc

   	/**
   	 * Calculates the cooling Lambda function for the metallicity and temperature of this halo.
   	 */
   	double logl = cooling_lambda_interpolator.get(lgTvir, zhot); //in cgs

	/**
	 * Calculate mean density for notional cooling profile.
	 */
	double nh_density  = mean_density(mhot, Rvir); //in units of cm^-3.

   	double tcool = 0;
   	double tcharac = 0;

   	/**
   	 * This corresponds to the very simple model of Croton06, in which the cooling time is assumed to be equal to the
   	 * dynamical timescale of the halo. With that assumption, and using an isothermal halo, the calculation of rcool and mcool
   	 * is trivial.
   	 */
   	if(parameters.model == GasCoolingParameters::CROTON06)
   	{
		tcool = darkmatterhalos->halo_dynamical_time(halo); //in Gyr.
		tcharac = tcool*constants::GYR2S; //in seconds.
	}

  	/**
   	 * This corresponds to the GALFORM-like model. Specifically here we are implementing the Benson & Bower (2010) model
   	 * for cooling, which calculated a time available for cooling depending on the integral of the the temperature
   	 * mass and cooling time history.
   	 */

   	if(parameters.model == GasCoolingParameters::BENSON10){
   		tcool = cooling_time(Tvir, logl,nh_density); //cooling time at notional density in Gyr.

   		/**
		 * Push back the cooling properties at this timestep.
		 */
   		subhalo.cooling_subhalo_tracking.tcooling.push_back(tcool);
   		subhalo.cooling_subhalo_tracking.temp.push_back(Tvir);
   		subhalo.cooling_subhalo_tracking.deltat.push_back(deltat);
   		//In the case of mass we convert back to comoving units.
   		subhalo.cooling_subhalo_tracking.mass.push_back(cosmology->physical_to_comoving_mass(mhot));

   		double integral = 0.0;//will save integral(T*M/tcool)

   		for(unsigned i=0;i<subhalo.cooling_subhalo_tracking.deltat.size();++i){
   			integral += subhalo.cooling_subhalo_tracking.temp[i]*subhalo.cooling_subhalo_tracking.mass[i]/subhalo.cooling_subhalo_tracking.tcooling[i]*subhalo.cooling_subhalo_tracking.deltat[i];
   		}
   		tcharac = integral/(Tvir*mhot/tcool); //available time for cooling in Gyr.
   	}

   	//I STILL NEED TO ADD A LIMIT TO THE TOTAL RADIATED ENERGY TO THE TOTAL THERMAL ENERGY OF THE HALO. SEE EQ. 18 AND 19 IN BENSON ET AL. (2010).

   	double r_cool = cooling_radius(mhot, Rvir, tcharac, logl, Tvir); //in Mpc.

   	if(r_cool < Rvir){
   		//cooling radius smaller than virial radius
   		coolingrate = 0.5*(r_cool/Rvir)*(mhot/tcool); //in Msun/Gyr.
   		r_cool = Rvir;
   	}
   	else {
   		//cooling radius larger than virial radius
   		coolingrate = mhot/tcool;
   	}

   	if(agnfeedback->parameters.model == AGNFeedbackParameters::GALFORM){
   		if(agnfeedback->parameters.alpha_cool > 0){
   			double tdyn_rcool = r_cool/vvir * constants::MPCKM2GYR; //Dynamical timescale at cooling radius.
    		/**
    		 * Calculate mean density at rcool.
    		 */
    		double nh_rcool  = density_shell(mhot, Rvir, r_cool); //in units of cm^-3.
    		double tcool_rcool = cooling_time(Tvir, logl, nh_rcool); //cooling time at rcool in Gyr.
    		double timescale_ratio = tcool_rcool/tdyn_rcool;
    		/**
       		 * Now we compare the cooling timescale with the dynamical timescale at rcool to
       		 * find out if galaxy is in hot halo mode. Only galaxies in the hot halo regime are
       		 * eligible for (radio) AGN fedback.
       		 */
       		if(timescale_ratio > 1.0/agnfeedback->parameters.alpha_cool){
       			//Halo is eligible for AGN feedback.
       			//Calculate cooling luminosity.
       			double Lcool = cooling_luminosity(logl, r_cool, Rvir, mhot);

       			if(Lcool < agnfeedback->parameters.f_edd * Ledd && Lcool > 0){
       				central_galaxy->smbh.macc_hh = agnfeedback->accretion_rate_hothalo_smbh(Lcool, central_galaxy->smbh.mass);
       				//now convert mass accretion rate to comoving units.
       				central_galaxy->smbh.macc_hh = cosmology->physical_to_comoving_mass(central_galaxy->smbh.macc_hh);
  		    		// set cooling rate to 0.
   					coolingrate = 0;
   				}
   				else{
   					central_galaxy->smbh.macc_hh = 0;
  				}
   			}//end if of hot halo mode.
		}// end if of AGN feedback model
	}// end if of GALFORM AGN feedback model.

    else if(agnfeedback->parameters.model == AGNFeedbackParameters::LGALAXIES){
    	//a pseudo cooling luminosity k*T/lambda(T,Z)
    	double Lpseudo_cool = constants::k_Boltzmann_erg * Tvir / std::pow(10.0,logl) / 1e40;
   		central_galaxy->smbh.macc_hh = agnfeedback->accretion_rate_hothalo_smbh(Lpseudo_cool, central_galaxy->smbh.mass);

		//now convert mass accretion rate to comoving units.
		central_galaxy->smbh.macc_hh = cosmology->physical_to_comoving_mass(central_galaxy->smbh.macc_hh);

    	//Mass heating rate from AGN in units of Msun/Gyr.
    	double mheatrate = agnfeedback->agn_bolometric_luminosity(central_galaxy->smbh.macc_hh) * 1e40 / (0.5*std::pow(vvir*KM2CM,2.0)) * MACCRETION_cgs_simu;

    	//modify cooling rate according to heating rate.
    	if(mheatrate < coolingrate){
    		coolingrate = coolingrate - mheatrate;
    	}
    	else{
    		/*CHECK: If mheatrate is > than cooling rate, then one needs to truncate the BH accretion rate to match that heating rate?*/
    		coolingrate = 0;
    	}

   	}
   	else{//In the case none of the two models above are included.
   		central_galaxy->smbh.macc_hh = 0;
   	}

   	/**
   	 * Now, we modify the SMBH mass and metals, and the hot gas mass accordingly.
   	 */
   	if(central_galaxy->smbh.macc_hh > 0){
   		//Now calculate new BH mass and metals due to gas accretion from hot halo.
		central_galaxy->smbh.mass += central_galaxy->smbh.macc_hh * deltat ;
		if(mhot > 0){
			central_galaxy->smbh.mass_metals += central_galaxy->smbh.macc_hh * deltat /mhot*mzhot;
		}

    	/**
    	 * Update hot halo gas properties as a response of how much the SMBH is accreting;
    	 */
    	subhalo.hot_halo_gas.mass -= central_galaxy->smbh.mass;
    	subhalo.hot_halo_gas.mass_metals -= central_galaxy->smbh.mass_metals;
   	}

   	if(coolingrate > 0){//perform calculations below ONLY if cooling rate >0.
   		/**
   		 * Convert cooling rate back to comoving units. This conversion is only necessary because mass was previously converted to physical.
   		 */
   		coolingrate = cosmology->physical_to_comoving_mass(coolingrate);
   		/**
   		 * Save properties of cooling gas in the halo gas component that tracks the cold gas.
   		 */
   		subhalo.cold_halo_gas.mass = coolingrate*deltat;
   		subhalo.cold_halo_gas.mass_metals = coolingrate*deltat/mhot*mzhot;//fraction of mass in the cold gas is the same as in metals.

   		/**
   		 * Update hot halo gas properties as a response of how much cooling there is in this timestep;
   		 */
   		subhalo.hot_halo_gas.mass -=subhalo.cold_halo_gas.mass;
   		subhalo.hot_halo_gas.mass_metals -= subhalo.cold_halo_gas.mass_metals;
   	}
   	else {
   		//avoid negative numbers.
   		coolingrate = 0;
   	}

 	if(subhalo.cold_halo_gas.mass < 0 or subhalo.cold_halo_gas.mass >1e17 or std::isnan(subhalo.cold_halo_gas.mass)){
		std::ostringstream os;
		os << halo << " has cold halo gas mass not well defined";
		throw invalid_data(os.str());
 	}

  	// Avoid negative values for the hot gas mass.
  	if(subhalo.hot_halo_gas.mass < constants::tolerance){
   		subhalo.hot_halo_gas.mass = 0;
   		subhalo.hot_halo_gas.mass_metals = 0;
   	}
  	if(subhalo.hot_halo_gas.mass_metals < constants::tolerance){
  		subhalo.hot_halo_gas.mass_metals = 0;
  	}

   	return coolingrate;

}

double GasCooling::cooling_time(double Tvir, double logl, double nh_density){

	using namespace constants;

	return 3.0*k_Boltzmann_erg*Tvir/(2.0*std::pow(10.0,logl)*nh_density)/GYR2S; //cooling time at notional density in Gyr.;
}

double GasCooling::mean_density(double mhot, double rvir){

	using namespace constants;

	return mhot*MSOLAR_g/(SPI * std::pow(rvir * MPC2CM,3.0)) / (M_Atomic_g * mu_Primordial); //in units of cm^-3.
}

double GasCooling::cooling_radius(double mhot, double rvir, double tcharac, double logl, double Tvir){

	using namespace constants;

	double pseudo_density = mhot*MSOLAR_g/(PI4*rvir*MPC2CM); //in units of gr/cm.

	double denominator_temp = 1.5 * k_Boltzmann_erg * Tvir * (M_Atomic_g*mu_Primordial) / std::pow(10.0, logl); //in cgs

	return std::pow((pseudo_density/denominator_temp*tcharac),0.5)/MPC2CM; //in Mpc.
}

double GasCooling::density_shell(double mhot, double rvir, double r) {

	using namespace constants;

	/**
	 * rho_shell as defined by an isothermal profile.
	 * Any other hot gas profile should modify rho_shell.
	 */
	return mhot*MSOLAR_g /PI4 /(rvir*MPC2CM) / std::pow(r*MPC2CM,2.0) / (M_Atomic_g*mu_Primordial); //in cgs.

}

double GasCooling::cooling_luminosity(double logl, double rcool, double rvir, double mhot){

	/**
	 *  This function calculates the cooling luminosity for a given cooling function
	 *  and a notional gas density profile. Units are returned in 10^40 erg/s.
	 *
	 */
	using namespace constants;

	if(rcool < rvir){

		/**
		 * For an isothermal profile, we define mass enclosed between rcool and rvir in csg.
		 */
		double mass_enclosed = mhot /PI4 /rvir * (rvir-rcool) * MSOLAR_g;

		//Define cooling luminosity in $10^{40} erg/s$.
		double Lcool = PI4 * std::pow(10.0,logl) * mass_enclosed / std::pow(10.0,40.0);

		return Lcool;
	}
	else{
		return 0;
	}
}

double GasCooling::disk_size_cooling(Subhalo &subhalo){

	//Do the basics first, and then more complicated stuff.
	return 0;
}


}  // namespace shark
