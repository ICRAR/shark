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
 *
 * A collection of numerical constants used throughout shark
 */

#ifndef SHARK_NUMERICAL_CONSTANTS_H_
#define SHARK_NUMERICAL_CONSTANTS_H_

namespace shark {

/**
 * A collection of numerical constants used throughout shark
 */
namespace constants {

	// Templated constexpr pow function to define constant powers with positive
	// integral exponents in our constants
	// Defined inside the constants namespace to avoid polluting the parennt shark
	// namespace with a pow function
	template <int __exp>
	constexpr double pow (double base) {
		return base * pow<__exp - 1>(base);
	}

	template <>
	constexpr double pow<1>(double base) {
		return base;
	}

	/** @name Numerical constants */
	// @{
	/** Square root of 2 */
	constexpr double SQRT2 = 1.4142135623730951;

	/** Square root of 3 */
	constexpr double SQRT3 = 1.7320508075688772;

	/** pi */
	constexpr double PI = 3.14159265358979323846;

	/** pi over 2 */
	constexpr double PIO2 = PI / 2.0;

	constexpr double SPI = 4.0 / 3.0 * PI;

	/** pi times 2 */
	constexpr double PI2 = 2.0 * PI;

	/** pi times 4 */
	constexpr double PI4 = 4.0 * PI;
	// @}

	/** @name Standard multipliers */
	// @{
	constexpr double HECTO = 1.0e+2;
	constexpr double KILO = 1.0e+3;
	constexpr double MEGA = 1.0e+6;
	constexpr double GIGA = 1.0e+9;
	constexpr double EPS3 = 1.0e-3;
	constexpr double EPS6 = 1.0e-6;
	// @}

	/** @name Unit conversions */
	// @{
	/** Number of ergs in a Joule */
	constexpr double J2ERG = 1.0e7;

	/** Number of seconds in a Gyr (based on the Julian year of exactly 365.25 days - Allen's Astrophysical Quantities, page 15) */
	constexpr double GYR2S = 3.15576e16;

	/** Number of metres in a Mpc (Particle Data Book 2002, page 6) */
	constexpr double MPC2M = 3.0856775807e22;

	/** Number of centimetres in a Mpc (Particle Data Book 2002, page 6) */
	constexpr double MPC2CM = MPC2M * HECTO;

	/** Number of centimetres in a Mpc (Particle Data Book 2002, page 6) cubed*/
	constexpr double MPC2CM_cube = pow<3>(MPC2CM);

	/** Number of cm in a km */
	constexpr double KM2CM = 1e5;

	/** Number of Joules in an erg */
	constexpr double ERG2J = 1.0e-7;

	/** Convert units of Mpc/(km/s) to Gyr */
	constexpr double MPCKM2GYR = MPC2M / GYR2S / KILO;
	// @}


	/** @name Properties of the Sun */
	// @{
	/** The mass of the Sun in kg */
	constexpr double MSOLAR = 1.9891e30;

	/**Solar luminosity in 10^40 ergs/s.**/

	constexpr double LSOLAR = 3.828e-7;

	/** The mass of the Sun in g */
	constexpr double MSOLAR_g = MSOLAR * KILO;

	/** Conversion of accretion rate from gr/s to Msun/Gyr */
	constexpr double MACCRETION_cgs_simu = 1 / MSOLAR_g * GYR2S;
	// @}

	/** @name Physical constants */
	// @{
	/** The gravitational constant in units of m^3/kg/s^2 */
	constexpr double G_SI = 6.67259e-11;

	/** The gravitational constant in units of (km/s)^2 Mpc/Msun */
	constexpr double G = G_SI * MSOLAR / MPC2M / pow<2>(KILO);

	/** Gravitational constant in units of cm/gr/s^2 */
	constexpr double G_cgs = 6.67259e-8;

	/** Boltzmann's constant in J/K (Particle Data Book 2002, page 5) */
	constexpr double k_Boltzmann = 1.3806503e-23;

	constexpr double k_Boltzmann_erg = k_Boltzmann * J2ERG;

	/** Speed of light in m/s (Allen's Astrophysical Quantities, page 8) */
	constexpr double c_light = 2.99792458e8;

	/** Speed of light in m/s (Allen's Astrophysical Quantities, page 8) */
	constexpr double c_light_km = 2.99792458e8 / KILO;

	/** Speed of light in cm/s (Allen's Astrophysical Quantities, page 8) */
	constexpr double c_light_cm = c_light * 100;

	/** Thomson cross section for zero energy photons in m^2 (Particle Data Book 2002, page 5) */
	constexpr double sigma_Thomson = 6.65245854e-29;

	/** Mass of unit atomic weight in kg (12C=12 scale) */
	constexpr double M_Atomic = 1.66053873e-27;

	/** Mass of unit atomic weight in g (12C=12 scale) */
	constexpr double M_Atomic_g = M_Atomic * KILO;

	/** Mass of hydrogen in units of M_Atomic (Particle Data Book 2002, page 283) */
	constexpr double Atomic_Mass_Hydrogen = 1.00794;

	/** Mass of helium in units of M_Atomic (Particle Data Book 2002, page 283) */
	constexpr double Atomic_Mass_Helium = 4.002602;

	constexpr double Pressure_Conv =  PIO2 * G * MSOLAR_g / pow<3>(MPC2CM) *  pow<2>(HECTO * KILO)/ k_Boltzmann_erg;
	// @}

	/** @name Cosmological constants */
	// @{
	/** The Hubble constant in units of h km/s/Mpc */
	constexpr double H0100 = 100.0;

	/** The Hubble constant in units of h/Gyr */
	constexpr double H0100PGYR = H0100 * KILO * GYR2S / MPC2M;

	/** Mass fraction of hydrogen in primordial plasma */
	constexpr double X_Hydrogen = 0.778;

	/** Mass fraction of helium in primordial plasma */
	constexpr double Y_Helium = 0.222;

	/** Mean atomic weight for fully ionized plasma of primordial composition */
	constexpr double mu_Primordial = 1.0 / (2.0 * X_Hydrogen / Atomic_Mass_Hydrogen + 3.0 * Y_Helium / Atomic_Mass_Helium);
	// @}

	/** @name Galaxy structural constants */
	// @{
	/** The half-mass radius of an exponential disk in units of the disk scale-length */
	constexpr double RDISK_HALF_SCALE = 1.678346990;
	// @}

	/** @name Gas content constants */
	// @{
	constexpr double Eddngtn_Lmnsty_Scale_Factor = 4.0 * PI * c_light * G_SI * MSOLAR * 1.0e-20 * M_Atomic * Atomic_Mass_Hydrogen / (sigma_Thomson * 1.0e20);
	constexpr double sigma_gas_mw = 2.5 * pow<2>(MEGA);
	// @}

	/** @name Miscellaneous */
	// @{
	/** Convert M_sun^2 kg^-2 ergs s^-1 cm^3 Mpc^-3 to 10^40 erg s^-1 */
	constexpr double lcool_conversion_factor = 1.0e-6 * pow<2>(1.0e-20 * MSOLAR / MPC2M / Atomic_Mass_Hydrogen / M_Atomic) / MPC2M;

	/** Maximum cooling luminosity in units of 10^40erg/s */
	constexpr double MAXLUM = 1e30;

	/** Define a tolerance for evaluating negative mass */
	constexpr double tolerance = 1e-10;

	/** Conversion between j/v and the half-mass radius from EAGLE */
	constexpr double EAGLEJconv = 0.835;
	//0.67714285714;
	// @}

} // namespace constants

}  // namespace shark

#endif // SHARK_NUMERICAL_CONSTANTS_H_
