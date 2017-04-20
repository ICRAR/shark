//
// Cosmological parameters used as inputs for SHArk
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
/*
 * Cosmology.h
 *
 *  Created on: 10Apr.,2017
 *      Author: clagos
 */

#ifndef SHARK_COSMOLOGY_H_
#define SHARK_COSMOLOGY_H_

#include <vector>

namespace shark {

/**
 * An element of the power spectrum
 */
struct PowerSpectrumElement {
	float k;
	float p;
};

/**
 * A power spectrum is simply a vector of power spectrum elements
 */
typedef std::vector<PowerSpectrumElement> PowerSpectrum;

/**
 * A set of cosmological parameters
 */
class CosmologicalParameters {

public:
	float OmegaM;
	float OmegaB;
	float OmegaL;
	float n_s;
	float sigma8;
	float Hubble_h;
	PowerSpectrum power_spectrum;
};

}  // namespace shark

#endif // SHARK_COSMOLOGY_H_