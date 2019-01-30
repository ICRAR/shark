//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2018
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

#ifndef SHARK_REIONISATION_H_
#define SHARK_REIONISATION_H_

#include <memory>
#include <sstream>
#include <utility>

#include "exceptions.h"
#include "options.h"

namespace shark {

class ReionisationParameters {

public:
	explicit ReionisationParameters(const Options &options);

	enum ReionisationModel {
		LACEY16 = 0,
		SOBACCHI13
	};

	double zcut = 0;
	double vcut = 0;
	double alpha_v = 0;
	ReionisationModel model = LACEY16;
};


/**
 * A class that checks whether halos are affected by reionisation or not
 */
class Reionisation {

public:
	explicit Reionisation(const ReionisationParameters &parameters);
	virtual ~Reionisation();

	/// Checks whether a halo of viral velocity @p v and redshift @p z is affected by reionisation
	virtual bool reionised_halo (double v, double z) const = 0;

protected:
	const ReionisationParameters &get_reionisation_params() const { return parameters; }

private:
	ReionisationParameters parameters;
};

/// The Lacey16 model of reionisation
class Lacey16Reionisation : public Reionisation {
public:
	using Reionisation::Reionisation;
	bool reionised_halo (double v, double z) const override;
};

/// The Sobacchi13 model of reionisation
class Sobacchi13Reionisation : public Reionisation {
public:
	using Reionisation::Reionisation;
	bool reionised_halo (double v, double z) const override;
};

using ReionisationPtr = std::shared_ptr<Reionisation>;

template <typename ...Ts>
ReionisationPtr make_reionisation(const ReionisationParameters &parameters, Ts&&...ts)
{
	if (parameters.model == ReionisationParameters::LACEY16) {
		return std::make_shared<Lacey16Reionisation>(parameters, std::forward<Ts>(ts)...);
	}
	else if (parameters.model == ReionisationParameters::SOBACCHI13) {
		return std::make_shared<Sobacchi13Reionisation>(parameters, std::forward<Ts>(ts)...);
	}

	std::ostringstream os;
	os << "Reionisation model " << parameters.model << " not currently supported";
	throw invalid_argument(os.str());
}

}//end namespace shark

#endif //SHARK_REIONISATION_H_