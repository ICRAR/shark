/*
 * star_formation.h
 *
 *  Created on: 17May,2017
 *      Author: clagos
 */

//
// star formation classes
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


#ifndef INCLUDE_STAR_FORMATION_H_
#define INCLUDE_STAR_FORMATION_H_

#include "options.h"

namespace shark {

class StarFormationParameters : public Options {

public:
	StarFormationParameters(const std::string &filename);

	int Molecular_BR_Law;
	double nu_sf;
	double Po;
	double beta_press;
	double Accuracy_SFeqs;
};


class StarFormation {

public:
	StarFormation(StarFormationParameters parameters);

	double star_formation_rate(double mcold, double mstars, double rgas, double rstars);

private:
	StarFormationParameters parameters;
};

}  // namespace shark

#endif /* INCLUDE_STAR_FORMATION_H_ */

