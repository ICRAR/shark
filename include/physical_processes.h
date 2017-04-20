//
// Description of physical processes used in SHArk
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

#ifndef SHARK_PHYSICAL_PROCESSES_H_
#define SHARK_PHYSICAL_PROCESSES_H_


namespace shark {

/**
 * This routine generates a set of ODEs based on the models
 * that we are using in SHARk.
 * There is a set of basic ODEs based on a minimum number
 * of physical processes that are included in SHArk.
 * If more sophisticated models are included this is the place
 * where the equations need to be extended.
 * The basic set of equations hold for a central galaxy without AGN feedback.
 * If either of the latter cases apply, these ODEs are modified accordingly.
 */

class Basic_Set_ODEs_Describing_Physics {
	/*The number of ODEs should be equal to the number of baryonic components in a subhalo.*/
	void star_formation();
	void recyclying();
	void stellar_feedback();
	void halo_cooling();
	void ejected_gas_reincorporation();
	void halo_gas_accretion();
};

class Modify_ODEs_Satellites: public Basic_Set_ODEs_Describing_Physics {
	void halo_cooling() override;
	void ejected_gas_reincorporation() override;
	void halo_gass_accretion();
};

class Modify_ODEs_AGN: public Basic_Set_ODEs_Describing_Physics {
	void black_hole_growth();
	void halo_cooling() override;
};

class Modify_ODEs_HaloEjection: public Basic_Set_ODEs_Describing_Physics{
	void ejected_gas_reincorporation() override;
	void outside_halo_ejection();
};

}  // namespace shark

#endif // SHARK_PHYSICAL_PROCESSES_H_