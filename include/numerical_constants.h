//
// A collection of numerical constants used throughout shark
//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2017
// Copyright by UWA (in the framework of the ICRAR)
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307  USA
//

#ifndef SHARK_NUMERICAL_CONSTANTS_H_
#define SHARK_NUMERICAL_CONSTANTS_H_

namespace shark {

/**
 * A collection of numerical constants used throughout SHArk
 */
namespace constants {


	// Templated constexpr pow function to define constant powers with positive
	// integral exponents in our constants
	// Defined inside the constants namespace to avoid polluting the parennt shark
	// namespace with a pow function
	template <int __exp>
	constexpr double pow (double base) {
		return base * pow<__exp - 1>(base);
	};

	template <>
	constexpr double pow<1>(double base) {
		return base;
	};

	constexpr float SQRT2=1.4142135624,LN2=0.6931471806,LN10=2.3025850930,ISO_FAC=2.5/LN10,log4=0.602059991;

	/** Euler's constant */
	constexpr float Eulers_Constant = 0.5772156649015328606;

	/* PI and friends */
	constexpr float PI = 3.1415926536;
	constexpr float logPI = 0.497149873;
	constexpr float PIO2 = PI / 2.0;
	constexpr float PIO4 = PI / 4.0;
	constexpr float PI4 = 4.0 * PI;
	constexpr float PI2 = 2.0 * PI;
	constexpr float logPI4 = logPI + log4;
	constexpr float PISQ = pow<2>(PI);
	constexpr float SQRTPI = 1.7724538509;
	constexpr float SQRT2PI = 2.5066282746;
	constexpr float SQRT2OPI = 0.7978845608;
	constexpr double DSQRTPI = 1.77245385090551602729816748334114518279754945612239;
	constexpr double DPI = 3.14159265358979323846264338327950288419716939937510;
	constexpr double SPI = DPI * 4.0 / 3.0;

	/*Standard multipliers.*/
	constexpr float MILLI=1.0e-3, logMILLI=-3.0;
	constexpr float DECA =1.0e+1, logDECA =+1.0;
	constexpr float HECTO=1.0e+2, logHECTO=+2.0;
	constexpr float KILO =1.0e+3, logKILO =+3.0;
	constexpr float MEGA =1.0e+6, logMEGA =+6.0;
	constexpr float GIGA =1.0e+9, logGIGA =+9.0;

	/* Small numbers (for tolerances usually).*/
	constexpr float EPS3=1.0e-3;
	constexpr float EPS4=1.0e-4;
	constexpr float EPS6=1.0e-6;

	/*Unit conversions.*/
	constexpr float J2ERG=1.0e7, logJ2ERG=7.0; /*Number of ergs in a Joule.*/
	constexpr float GYR2S=3.15576e16, logGYR2S=16.499103967; /*The number of seconds in a Gyr (based on the Julian year
	                                                of exactly 365.25 days - Allen's Astrophysical Quantities, page 15)*/
	constexpr float GYR2YR=1.0e9; /*The number of years in a Gyr.*/
	constexpr float YR2S=GYR2S/GYR2YR; /*The number of seconds in a year.*/
	constexpr float Ang2M=1.0e-10, logAng2M=-10.0; /*The number of metres in an Angstrom.*/
	constexpr float MPC2M=3.0856775807e22, sqrtMPC2M=1.75660968365e11; /*The number of metres in a Mpc (Particle Data Book 2002,
	                                                           page 6).*/
	constexpr float KPC2M=MPC2M/KILO; /*The number of metres in a kpc (Particle Data Book 2002, page 6).*/
	constexpr float PC2M=3.0856775807e16, logPC2M=16.489350545; /*The number of metres in a pc (Particle Data Book 2002, page 6).*/
	constexpr float MPC2CM=MPC2M*HECTO; /*The number of centimetres in a Mpc (Particle Data Book 2002, page 6).*/
	constexpr float PC2CM=PC2M*HECTO; /*The number of centimetres in a pc (Particle Data Book 2002, page 6).*/
	constexpr float KM2M=1.0e3; /* The number of metres in a kilometre.*/
	constexpr float BARN2M2=1.0e-28; /* Convert barns to m^2.*/
	constexpr float KMS2MPCGYR=KM2M*GYR2S/MPC2M; /*Convert velocity in km/s to Mpc/Gyr.*/
	constexpr float M2CM=100.0, logM2CM=2.0; /*Number of cm in a m.*/
	constexpr float KM2CM = M2CM*KM2M; /*Number of cm in a km.*/
	constexpr float M32CM3=pow<3>(M2CM); /*Number of cm^3 in a m^3.*/
	constexpr float CM32M3=1.0/M32CM3; /*Number of m^3 in a cm^3.*/
	constexpr float eV2ERG=1.60217733e-12, logeV2ERG=-11.7952894176; /*Convert electron volts to ergs.*/
	constexpr float eV2J=1.60217733e-19; /*Convert electron volts to Joules.*/
	constexpr float ERG2J=1.0e-7; /*Number of Joules in an erg.*/
	constexpr float MPCKM2GYR=MPC2M/GYR2S/KM2M; /*Convert units of Mpc/(km/s) to Gyr.*/
	constexpr float RAD2AS=180.0*60.0*60.0/PI, logRAD2AS=5.81157500587-logPI; /*Convert radians to arcseconds.*/
	constexpr float RAD2AM=180.0*60.0/PI; /*Convert radians to arcminutes.*/
	constexpr float RAD2DEGREES=180.0/PI;
	constexpr float Right_Angle_Degrees=90.0; /*Number of degrees in a right angle.*/
	constexpr float Semi_Circle_Degrees=180.0; /*Number of degrees in a semi-circle.*/
	constexpr float MPC2ASat10pc=(MEGA/DECA)*RAD2AS, logMPC2ASat10pc=logMEGA-logDECA+logRAD2AS; /*Converts Mpc to arcseconds for
	                                                                                    objects placed at a distance of 10pc.*/

	/*Properties of the Sun.*/
	/*The mass of the Sun in kg (Allen's Astrophysical Quantities, page 12).*/
	constexpr float MSOLAR=1.9891e30, MSOLAR_1030kg=MSOLAR/1.0e30, sqrtMSOLAR=1.4103545653e15, MSOLAR_g=MSOLAR*KILO, logMSOLAR = 30.298656617391146;
	constexpr float LSUN=3.845e26, logLSUN=26.5848963; /*Bolometric luminosity of the Sun in W (Allen's Astrophysical Quantities, page 340).*/
	constexpr float LSUN_ERG=LSUN*J2ERG; /*Bolometric luminosity of the Sun in erg/s (Allen's Astrophysical Quantities, page 340).*/
	constexpr float LOG_LSUN_ERG=logLSUN+logJ2ERG; /*Bolometric luminosity of the Sun in erg/s (Allen's Astrophysical Quantities, page 340).*/
	constexpr float LSUN_BC=3.826e26, logLSUN_BC=26.58274497; /*Bolometric luminosity of the Sun as used by Bruzual & Charlot [W]*/
	constexpr float LSUN40_BC=3.826e-7; /*Bolometric luminosity of the Sun as used by Bruzual & Charlot [10^40 ergs/s].*/
	constexpr float X_Hydrogen_Solar=0.707; /*Mass fraction of hydrogen in Solar plasma (Allen's Atrophysical Quantities, page 28).*/
	constexpr float Y_Helium_Solar=0.274; /*Mass fraction of helium in Solar plasma (Allen's Atrophysical Quantities, page 28).*/
	constexpr float Z_Metals_Solar=0.0189; /*Mass fraction of metals in Solar plasma (Allen's Atrophysical Quantities, page 28).*/
	constexpr float MACCRETION_cgs_simu = 1/MSOLAR_g*GYR2S; /*Conversion of accretion rate from gr/s to Msun/Gyr.*/

	/* Physical constants.*/
	constexpr float G_SI=6.67259e-11, sqrtG_SI=8.16859228998e-6; /*The gravitational constant in units of m^3/kg/s^2 (Allen's Astrophysical Quantities, page 8).*/
	constexpr float G=G_SI*MSOLAR/MPC2M/pow<2>(KM2M), sqrtG=sqrtG_SI*sqrtMSOLAR/sqrtMPC2M/KM2M; /*The gravitational constant in units of (km/s)^2 Mpc/Msun.*/
	constexpr float G_cgs = 6.67259e-8; /*Gravitational constant in units of cm/gr/s^2.*/
	constexpr float G_MPCGYR=G_SI*MSOLAR*(GYR2S/MPC2M)/KM2M/MPC2M; /*The gravitational constant in units of km/s Mpc^2 Msun^-1 Gyr^-1*/
	constexpr float G_MPCGYR2=G_SI*MSOLAR*pow<2>(GYR2S/MPC2M)/MPC2M; /*The gravitational constant in units of Mpc^3 Msun^-1 Gyr^-2*/
	constexpr float G_GYRKMS3=G_SI*MSOLAR/GYR2S/pow<3>(KM2M); /*The gravitational constant in units of Gyr Msol^-1 km^3 s^-3*/
	constexpr double G_KMSPC=G*1e6; /*The gravitational constant in units of (km/s)^2 pc/Msun. Needs to be* in double precision because it is used in the sn_dynamical_feedback routines.*/
	constexpr float GPI=G*PI*MPC2M/GYR2S/KM2M; /*G*pi in units of Gyr (km/s)^3 Msun^-1*/
	constexpr float k_Boltzmann=1.3806503e-23; /*Boltzmann's constant in J/K (Particle Data Book 2002, page 5).*/
	constexpr float k_Boltzmann_erg=k_Boltzmann*J2ERG, logk_Boltzmann=-22.859916308, sqrtk_Boltzmann=3.715710295489e-12;
	constexpr float K2eV=k_Boltzmann/eV2J; /*Convert temperature to eV.*/
	constexpr float c_light=2.99792458e8, logc_light=8.47682070; /*Speed of light in m/s (Allen's Astrophysical Quantities, page 8).*/
	constexpr float c_light_Angstroms=c_light/Ang2M, logc_light_Angstroms=logc_light-logAng2M; /*Speed of light in Angstroms/s (Allen's Astrophysical Quantities, page 8)*/
	constexpr float c_light_cm=c_light*M2CM; /*Speed of light in cm/s (Allen's Astrophysical Quantities, page 8).*/
	constexpr float h_Planck=6.6260755e-34; /*Planck's constant in J s (Allen's Astrophysical Quantities, page 8).*/
	constexpr float h_Planck_eV=h_Planck/eV2J; /*Planck's constant in eV s (Allen's Astrophysical Quantities, page 8).*/
	constexpr float h_Planck_erg=h_Planck*J2ERG; /*Planck's constant in ergs s (Allen's Astrophysical Quantities, page 8).*/
	constexpr float sigma_Thomson=6.65245854e-29; /*Thomson cross section for zero energy photons in m^2 (Particle Data Book 2002, page 5).*/
	constexpr float Q_Electron=1.602176487e-19; /*Charge of the electron in Coulombs (NIST: http://physics.nist.gov/cgi-bin/cuu/Value?e).*/
	constexpr float M_Electron=9.10938188e-31; /*Mass of the electron in kg (Particle Data Book 2002, page 4).*/
	constexpr float M_Atomic=1.66053873e-27, sqrtM_Atomic=4.07497083425e-14; /*Mass of unit atomic weight in kg (12C=12 scale);*/
	constexpr float M_Atomic_g=M_Atomic*KILO; /*Mass of unit atomic weight in g (12C=12 scale).*/
	constexpr float a_Radiation=8.0*pow<5>(PI)*k_Boltzmann*pow<3>(k_Boltzmann/c_light/h_Planck)/15.0; /*Radiation constant (J m^-3 K^-4).*/
	constexpr float Atomic_Mass_Hydrogen=1.00794, sqrtAtomic_Mass_Hydrogen=1.00396215; /*Mass of hydrogen in units of M_Atomic (Particle Data Book 2002, page 283).*/
	constexpr float Atomic_Mass_Helium=4.002602; /*Mass of helium in units of M_Atomic (Particle Data Book 2002, page 283).*/

	constexpr double Pressure_Conv =  PIO2 * G * MSOLAR_g / pow<3>(MPC2CM) *  pow<2>(HECTO * KILO)/ k_Boltzmann_erg;

	/*Electromagnetism.*/
	constexpr float Permeability_of_Free_Space=4.0e-7*PI; /*Permeability of free space in units of N A^-2 (definition).*/
	constexpr float Permittivity_of_Free_Space=1.0/Permeability_of_Free_Space/pow<2>(c_light); /*Permittivity of free space in units of Coulombs^2/N/m^2 (definition).*/

	/*Cosmological constants.*/
	constexpr float H0100=100.0; /*The Hubble constant in units of h km/s/Mpc.*/
	constexpr float H0100PGYR=H0100*KM2M*GYR2S/MPC2M; /*The Hubble constant in units of h/Gyr.*/
	constexpr float RHOCRIT=3.0*pow<2>(H0100)/8.0/PI/G, sqrtRHOCRIT=0.61237243570*H0100/SQRTPI/sqrtG; /*Critical density of the Universe (3*H0^2/8*PI*G) in h^2 Msun/Mpc^3.*/

	constexpr float X_Hydrogen=0.778; /*Mass fraction of hydrogen in primordial plasma.*/
	constexpr float Y_Helium  =0.222; /*Mass fraction of helium in primordial plasma.*/
	constexpr float Z_Metals  =5.36e-10; /*Mass fraction of metals in primordial plasma.*/
	constexpr float mu_Primordial=1.0/(2.0*X_Hydrogen/Atomic_Mass_Hydrogen+3.0*Y_Helium/Atomic_Mass_Helium); /*Mean atomic weight for fully ionized plasma of primordial composition.*/

	constexpr float M8CRIT=RHOCRIT*4.0*PI*pow<3>(8.0)/3.0; /*Mass in a sphere of 8Mpc/h radius in a critical density Universe.*/
	constexpr float KHORIZON=H0100*KM2M/c_light; /*Defined as H_0/c in h Mpc^-1.*/

	/*Dark matter halo parameters.*/
	constexpr float v_Halo_to_T_Virial=0.5e6*M_Atomic/k_Boltzmann; /*Conversion factor from virial velocity (km/s) to virial temperature (K).*/
	constexpr float Delta200=200.0, sqrtDelta200=14.142135623; /*Mean density contrast of a halo at radius r_200.*/
	constexpr float DeltaEdS=18.0*PI*PI; /*Density contrast at virial radius of a halo in an Einstein-de Sitter universe.*/
	constexpr float fx_peak_NFW=0.2162165956; /*The peak value of the function [ln(1+x)-x/(1+x)]/x as appears in the rotation curve of the NFW halo.*/

	/*Galaxy structural constants.*/
	constexpr float kstrc_1=0.6251543028; /*Constant relating V_c(rdisk)^2 to GMdisk/rdisk in the disk plane.*/
	constexpr float kstrc_2=1.191648617; /*Constant relating disk angular momentum to r*V_c derived assuming a flat rotation curve.*/
	constexpr float RDISK_HALF_SCALE=1.678346990; /*The half-mass radius of an exponential disk in units of the disk scale-length.*/
	constexpr float EM_RDISK_HALF_SCALE_O2=0.432067481; /*exp(-RDISK_HALF_SCALE/2).*/
	constexpr float SURF_DENS_NORM_HALF=2.97114e-2; /*exp(-RDISK_HALF_SCALE)/2/PI.*/

	/*gas content.*/
	constexpr float corr_factor_He=1.35; /*correction factor to account for helium when only hydrogen is given: H=Mcold/corr_factor_He.*/
	constexpr float PressureConst=4.33e-12; /*Constant converting the gravity and Boltzmann constants from the MKS system to the units necessary to obtain the pressure in units of K*cm^-3 (cgs).*/
	constexpr double Eddngtn_Lmnsty_Scale_Factor=4.0*PI*c_light*G_SI*MSOLAR*1.0e-20*M_Atomic*Atomic_Mass_Hydrogen/(sigma_Thomson*1.0e20), Eddngtn_Mdot_Constant = 0.1;
	constexpr double sigma_gas_mw = 2.5 * std::pow(MEGA,2.0);
	constexpr double lcool_conversion_factor = 1.0e-6 * std::pow(1.0e-20 * MSOLAR / MPC2M / Atomic_Mass_Hydrogen / M_Atomic, 2.0) / MPC2M; // Convert M_sun^2 kg^-2 ergs s^-1 cm^3 Mpc^-3 to 10^40 erg s^-1;

	/*define a tolerance for evaluating negative mass*/
	constexpr float tolerance = 1e-10;

	/*define a tolerance of gas mass for evaluating starbursts*/
	constexpr float tolerance_mass = 1e3;

	/* define a maximum cooling luminosity in units of 10^40erg/s.*/
	constexpr float MAXLUM = 1e30;

};

}  // namespace shark

#endif // SHARK_NUMERICAL_CONSTANTS_H_
