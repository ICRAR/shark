//
// Parameters used as input for shark
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

#ifndef SHARK_PARAMETERS_H_
#define SHARK_PARAMETERS_H_

#include <vector>

namespace shark {

class StarFormationParameters {
public:
	int Molecular_BR_Law;
	double nu_sf;
	double Po;
	double beta_press;
	double Accuracy_SFeqs;
};


class reionisation_parameters {
public:
	double vcut;
	double zcut;
};

class RecyclingParameters {
public:
	double yield;
	double recycle;
	double zsun;
};

class Parameters {
public:
	std::vector<int> writing_outputs;
	//flags of physical processes that are on.
	int ReionizationOn;
	int SupernovaFeedbackOn;
	int DiskInstabilityOn;
	int AGNFeedbackOn;
	int SFprescription;
};

}  // namespace shark

#endif // SHARK_PARAMETERS_H_
