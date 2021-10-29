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
#include <iterator>
#include <map>
#include <numeric>
#include <tuple>

#include "cosmology.h"
#include "data.h"
#include "galaxy.h"
#include "gas_cooling.h"
#include "halo.h"
#include "logging.h"
#include "numerical_constants.h"
#include "reincorporation.h"
#include "subhalo.h"
#include "utils.h"

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
		std::vector<double> first_keys = get_keys(first_measurement);
		std::vector<double> these_keys = get_keys(first_measurement);

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

GasCoolingParameters::GasCoolingParameters(const Options &options)
{
	options.load("gas_cooling.model", model, true);
	options.load("gas_cooling.lambdamodel", lambdamodel, true);
	options.load("gas_cooling.pre_enrich_z", pre_enrich_z);
	options.load("gas_cooling.tau_cooling", tau_cooling);

	auto cooling_tables_dir = get_static_data_filepath("cooling");
	tables_idx metallicity_tables = find_tables(cooling_tables_dir);
	load_tables(cooling_tables_dir, metallicity_tables);

}

GasCoolingParameters::tables_idx GasCoolingParameters::find_tables(
	const std::string &cooling_tables_dir)
{
	//read cooling tables and load values in
	std::string prefix;
	if (lambdamodel == CLOUDY) {
		prefix = "C08.00_";
	}
	else if (lambdamodel == SUTHERLAND) {
		prefix = "S93_";
	}
	else {
		throw invalid_argument("Cooling model is not valid");
	}

	std::string tables = cooling_tables_dir + "/" + prefix + "tables.txt";

	// Collect metallicity tables information
	LOG(debug) << "Reading metallicity table index" << tables;
	std::string line;
	std::map<double, std::string> metallicity_tables;
	std::ifstream f = open_file(tables);
	while ( std::getline(f, line) ) {

		trim(line);
		if (empty_or_comment(line)) {
			continue;
		}

		double metallicity;
		std::string table_fname;
		std::istringstream iss(line);
		iss >> metallicity >> table_fname;
		trim(table_fname);

		metallicity_tables[metallicity] = table_fname;
	}
	f.close();

	return metallicity_tables;
}

void GasCoolingParameters::load_tables(
	const std::string &cooling_tables_dir,
	const tables_idx &metallicity_tables)
{

	// Populate cooling_table
	for(auto &kv: metallicity_tables) {

		double metallicity = std::get<0>(kv);
		const std::string fname = cooling_tables_dir + "/" + std::get<1>(kv);

		LOG(debug) << "Reading table " << fname << " for metallicity " << metallicity;

		std::ifstream f = open_file(fname);
		std::string line;

		std::map<double, double> measurements;
		while ( std::getline(f, line) ) {

			trim(line);
			if (empty_or_comment(line)) {
				continue;
			}

			double t, ne, nh, nt, logl;
			std::istringstream iss(line);
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
	auto lvalue = lower(value);
	if (lvalue == "cloudy") {
		return GasCoolingParameters::CLOUDY;
	}
	else if (lvalue == "sutherland") {
		return GasCoolingParameters::SUTHERLAND;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value << ". Supported values are cloudy and sutherland";
	throw invalid_option(os.str());
}


template <>
GasCoolingParameters::CoolingModel
Options::get<GasCoolingParameters::CoolingModel>(const std::string &name, const std::string &value) const {
	auto lvalue = lower(value);
	if (lvalue == "croton06") {
		return GasCoolingParameters::CROTON06;
	}
	else if (lvalue == "benson10") {
		return GasCoolingParameters::BENSON10;
	}
	std::ostringstream os;
	os << name << " option value invalid: " << value << ". Supported values are croton06 and benson10";
	throw invalid_option(os.str());
}

GasCooling::GasCooling(GasCoolingParameters parameters,
		StarFormationParameters params_sf,
		ExecutionParameters exec_params,
		ReionisationPtr reionisation,
		CosmologyPtr cosmology,
		AGNFeedbackPtr agnfeedback,
		DarkMatterHalosPtr darkmatterhalos,
		ReincorporationPtr reincorporation,
		EnvironmentPtr environment) :
	parameters(parameters),
	params_sf(params_sf),
	exec_params(exec_params),
	reionisation(std::move(reionisation)),
	cosmology(std::move(cosmology)),
	agnfeedback(std::move(agnfeedback)),
	darkmatterhalos(std::move(darkmatterhalos)),
	reincorporation(std::move(reincorporation)),
	environment(std::move(environment)),
	cooling_lambda_interpolator(parameters.cooling_table.get_temperatures(), parameters.cooling_table.get_metallicities(), parameters.cooling_table.get_lambda())
{
	//no-opt
}

double GasCooling::cooling_rate(Subhalo &subhalo, Galaxy &galaxy, double z, double deltat) {

	using namespace constants;

	double coolingrate = 0;

	//Define host halo
	auto halo = subhalo.host_halo;

	halo->cooling_rate = 0;

	if(subhalo.subhalo_type == Subhalo::SATELLITE){
		//Compute how much hot gas there is in this satellite_subhalo based on the environmental processes applied to it.
		environment->process_satellite_subhalo_environment(subhalo, *subhalo.host_halo->central_subhalo);
	}

	// If galaxy is type 2, then they don't have a hot halo.
	if ( galaxy.galaxy_type == Galaxy::TYPE2) {
		return 0;
	}

	// Define main galaxy, which would accrete the cooled gas if any.
	Galaxy *central_galaxy = &galaxy;

	// If subhalo does not have a main galaxy (which could happen in satellite subhalos), or subhalo does not have a hot halo, return 0.
	if(subhalo.hot_halo_gas.mass <= 0){
		return 0;
	}

	// Check if user is running the code to ignore galaxy formation in late, massive forming halos, and if so,
	// check whether this is one of the halos we have to ignore.
	if(exec_params.ignore_late_massive_halos && halo->ignore_gal_formation){
		return 0;
	}


	/**
	 * Calculate reincorporated mass and metals.
	 */
	double mreinc_mass = reincorporation->reincorporated_mass(*halo, subhalo, z, deltat);

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
		subhalo.ejected_galaxy_gas.mass_metals = 0;
	}
	if(subhalo.ejected_galaxy_gas.mass_metals < 0){
		subhalo.ejected_galaxy_gas.mass_metals = 0;
	}


	/**
	 * Include accreted gas only if the subhalo is central.
	 */
	if (subhalo.subhalo_type == Subhalo::CENTRAL) {
		// Calculate maximum accreted mass allowed by the universal baryon fraction.
		float max_allowed_baryon_accreted = halo->Mvir * cosmology->universal_baryon_fraction() - halo->total_baryon_mass();

		if (max_allowed_baryon_accreted > 0) {
			// Add up accreted mass and metals.
			auto accreted_mass = std::min(subhalo.accreted_mass, max_allowed_baryon_accreted);
			subhalo.hot_halo_gas.mass += accreted_mass;
			subhalo.hot_halo_gas.mass_metals += accreted_mass * parameters.pre_enrich_z;
		}
	}

	/**
	* Plant black hole seed if necessary.
	*/
	agnfeedback->plant_seed_smbh(*halo);

	// Calculate Eddington luminosity of BH in galaxy.

	double Ledd = agnfeedback->eddington_luminosity(central_galaxy->smbh.mass);

	/**
	* Test for subhalos that are affected by reionisation
	*/
	auto reionised_halo = reionisation->reionised_halo(subhalo.Vvir, z);

	if(reionised_halo){
		return 0;
	}


	//Assume hot halo has the same specific angular momentum of DM halo.
	subhalo.hot_halo_gas.sAM = subhalo.L.norm() / subhalo.Mvir;

	/**
	* We need to convert masses and velocities to physical units before proceeding with calculation.
	*/
	double mhot = cosmology->comoving_to_physical_mass(subhalo.hot_halo_gas.mass + subhalo.cold_halo_gas.mass);
	double mhot_ejec = cosmology->comoving_to_physical_mass(subhalo.ejected_galaxy_gas.mass);
	double mzhot = cosmology->comoving_to_physical_mass(subhalo.hot_halo_gas.mass_metals + subhalo.cold_halo_gas.mass_metals);

	double vvir = subhalo.Vvir;

	double zhot = 0;
	if(mhot > 0){
		zhot = (mzhot/mhot);
	}

	// Check for undefined cases.
	if(mhot < 0 || mhot >1e17 || std::isnan(mhot)){
		std::ostringstream os;
		os << halo << " has hot halo gas mass not well defined";
		throw invalid_data(os.str());
	}

	// If galaxy is central then redefine vmax according to subhalo Vmax.
	if(galaxy.galaxy_type == Galaxy::CENTRAL) {
		galaxy.vmax  = subhalo.Vcirc;
	}

	double Tvir   = 97.48*std::pow(vvir,2.0); //in K.
	double lgTvir = log10(Tvir); //in K.
	double Rvir   = cosmology->comoving_to_physical_size(darkmatterhalos->halo_virial_radius(halo, z), z);//physical Mpc

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
		tcool = parameters.tau_cooling * darkmatterhalos->halo_dynamical_time(halo, z); 
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
		tcharac = integral/(Tvir*mhot/tcool) *constants::GYR2S; //available time for cooling in seconds.
	}

	double r_cool = cooling_radius(mhot, Rvir, tcharac, logl, Tvir); //in physical Mpc.

	if(r_cool < Rvir){
		//cooling radius smaller than virial radius
		coolingrate = 0.5*(r_cool/Rvir)*(mhot/tcool); //in Msun/Gyr.
	}
	else {
		//cooling radius larger than virial radius, and set cooling radius to virial radius.
		coolingrate = mhot/tcool;
				r_cool = Rvir;
	}

	if(agnfeedback->parameters.model == AGNFeedbackParameters::BOWER06){
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
				//Calculate cooling luminosity using total hot gas mass as in GALFORM.
				double Lcool = cooling_luminosity(logl, r_cool, Rvir, mhot + mhot_ejec);
	
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
	}// end if of BOWER06 AGN feedback model.

	else if(agnfeedback->parameters.model == AGNFeedbackParameters::BRAVO19 || agnfeedback->parameters.model == AGNFeedbackParameters::CROTON16){
		//a pseudo cooling luminosity k*T/lambda(T,Z)
		double Lpseudo_cool = constants::k_Boltzmann_erg * Tvir / std::pow(10.0,logl) / 1e40;
		central_galaxy->smbh.macc_hh = agnfeedback->accretion_rate_hothalo_smbh(Lpseudo_cool, central_galaxy->smbh.mass);

		//now convert mass accretion rate to comoving units.
		central_galaxy->smbh.macc_hh = cosmology->physical_to_comoving_mass(central_galaxy->smbh.macc_hh);

		//Mass heating rate from AGN in units of Msun/Gyr.
		double mheatrate = 0;
		if(agnfeedback->parameters.model == AGNFeedbackParameters::BRAVO19){
			// decide whether this halo is in a quasi-hydrostatic regime or not.
			bool hothalo = quasi_hydrostatic_halo(mhot, std::pow(10.0,logl), nh_density, halo->Mvir, Tvir, Rvir, z);

			// radio mode feedback only applies in situations where there is a hot halo
			if(hothalo){
				mheatrate = agnfeedback->agn_mechanical_luminosity(central_galaxy->smbh.macc_sb + central_galaxy->smbh.macc_hh , central_galaxy->smbh.mass) * agnfeedback->parameters.kappa_radio
						* 1e40 / (0.5 * std::pow(vvir * KM2CM,2.0)) * MACCRETION_cgs_simu;
			}
		}
		else if(agnfeedback->parameters.model == AGNFeedbackParameters::CROTON16){
			mheatrate = agnfeedback->agn_bolometric_luminosity(central_galaxy->smbh.macc_hh,central_galaxy->smbh.mass) * 1e40 /
					(0.5 * std::pow(vvir * KM2CM,2.0)) * MACCRETION_cgs_simu;
		}

		// Calculate heating radius
		double rheat = mheatrate/coolingrate * r_cool;

		double r_ratio = rheat/r_cool;
		if(agnfeedback->parameters.model == AGNFeedbackParameters::CROTON16){
			if(subhalo.cooling_subhalo_tracking.rheat < rheat){
				subhalo.cooling_subhalo_tracking.rheat = rheat;
			}

			r_ratio = subhalo.cooling_subhalo_tracking.rheat/r_cool;
		}

		if(r_ratio > agnfeedback->parameters.alpha_cool){
			r_ratio = 1;
			//Redefine mheatrate and macc_h accordingly.
			mheatrate = r_ratio * coolingrate;
			central_galaxy->smbh.macc_hh = agnfeedback->accretion_rate_hothalo_smbh_limit(mheatrate,vvir);
		}

		//modify cooling rate according to heating rate.
		coolingrate = (1 - r_ratio) * coolingrate;
		if(coolingrate < 0){
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

		double delta_mass_bh = central_galaxy->smbh.macc_hh * deltat;
		double delta_metals_bh = 0;
		if(mhot > 0) {
			delta_metals_bh = delta_mass_bh/mhot * mzhot;
		}

		central_galaxy->smbh.mass += delta_mass_bh;
		central_galaxy->smbh.mass_metals += delta_metals_bh;

		/**
		* Update hot halo gas properties as a response of how much the SMBH is accreting;
		*/
		subhalo.hot_halo_gas.mass -= delta_mass_bh;
		subhalo.hot_halo_gas.mass_metals -= delta_metals_bh;
	}

	if(coolingrate > 0){//perform calculations below ONLY if cooling rate >0.
		/**
		* Convert cooling rate back to comoving units. This conversion is only necessary because mass was previously converted to physical.
		*/
		coolingrate = cosmology->physical_to_comoving_mass(coolingrate);

		double mcooled = coolingrate * deltat;

		// Limit cooled mass to the amount available in the halo and adjust cooling rate accordingly.
		if(mcooled > subhalo.hot_halo_gas.mass){
			mcooled = subhalo.hot_halo_gas.mass;
			coolingrate = mcooled / deltat;
		}

		/**
		* Save properties of cooling gas in the halo gas component that tracks the cold gas.
		*/
		subhalo.cold_halo_gas.mass += mcooled;
		subhalo.cold_halo_gas.mass_metals += mcooled * mzhot/mhot;//fraction of mass in the cold gas is the same as in metals.

		/**
		* Update hot halo gas properties as a response of how much cooling there is in this timestep;
		*/
		subhalo.hot_halo_gas.mass -= mcooled;
		subhalo.hot_halo_gas.mass_metals -= mcooled * mzhot/mhot;
	}
	else {
		//avoid negative numbers.
		coolingrate = 0;
	}

	// check for undefined values.
	if(subhalo.cold_halo_gas.mass < 0 || subhalo.cold_halo_gas.mass >1e17 || std::isnan(subhalo.cold_halo_gas.mass)){
		std::ostringstream os;
		os << halo << " has cold halo gas mass not well defined";
		throw invalid_data(os.str());
	}

	// check for undefined values.
	if(subhalo.hot_halo_gas.mass < 0 || subhalo.hot_halo_gas.mass >1e17 || std::isnan(subhalo.hot_halo_gas.mass)){
		std::ostringstream os;
		os << halo << " has hot halo gas mass not well defined";
		throw invalid_data(os.str());
	}

	// Avoid negative values for the hot gas mass.
	if(subhalo.hot_halo_gas.mass < constants::tolerance){
		subhalo.hot_halo_gas.restore_baryon();
	}
	if(subhalo.hot_halo_gas.mass_metals < constants::tolerance){
		subhalo.hot_halo_gas.mass_metals = 0;
	}

	// Save net cooling rate.
	halo->cooling_rate = coolingrate;

	if(coolingrate > 0){
		// define cooled gas angular momentum.
		darkmatterhalos->cooling_gas_sAM(subhalo, z);
	}
	else{
		subhalo.cold_halo_gas.sAM = 0;
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
		* For an isothermal profile, we define a small core radius.
		*/
		double rcore = 0.01 * rvir;

		double r1 = rvir/rcore;
		double r2 = rcool/rcore;

		//Define cooling luminosity in $10^{40} erg/s$.
		double func1 = std::atan(r1) - r1/(std::pow(r1,2.0) + 1);
		double func2 = std::atan(r2) - r2/(std::pow(r2,2.0) + 1);

		double ave_pseudo_density = std::pow(mhot, 2.0) / std::pow(rcore,3.0); //in Msun^2/Mpc^3.
		double factor_geometry = (func1 - func2)/ std::pow(r1 - std::atan(r1), 2.0);

		double Lcool = lcool_conversion_factor / (8.0*PI) * std::pow(10.0,logl) * ave_pseudo_density  * factor_geometry;

		return Lcool;
	}
	else{
		return 0;
	}
}

bool GasCooling::quasi_hydrostatic_halo(double mhot, double lambda, double nh_density, double mass, double Tvir, double rvir, double redshift){
		/**
		 *  This function uses the model of Correa et al. (2018) to determine if a hot halo has formed or not.
		 **/

		using namespace constants;

		// cooling rate in cgs.
		double gamma_cool = mhot * MSOLAR_g * lambda * nh_density / (M_Atomic_g * mu_Primordial * 0.1);

		double omega_term = std::sqrt(cosmology->parameters.OmegaM * std::pow(redshift + 1.0, 3.0) + cosmology->parameters.OmegaL);

		// growth rate of halo in Msun/Gyr.
		double mdot = 71.6 * GIGA * (cosmology->comoving_to_physical_mass(mass) / 1e12) * (cosmology->parameters.Hubble_h/0.7)  * (-0.24 + 0.75 * (redshift + 1.0)) * omega_term;

		// heating rate in cgs.
		double gamma_heat = 1.5 * k_Boltzmann_erg * Tvir / (M_Atomic_g * mu_Primordial) * cosmology->universal_baryon_fraction() * mdot / MACCRETION_cgs_simu * 0.8333;

		double ratio = gamma_cool/gamma_heat;
		if(ratio <  agnfeedback->parameters.hot_halo_threshold){
			return true;
		}
		else{
			return false;
		}

}


}  // namespace shark
