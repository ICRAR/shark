/*
 * reionisation.h
 *
 *  Created on: 14Jun.,2017
 *      Author: clagos
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
	ReionisationParameters(const Options &options);

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
	Reionisation(const ReionisationParameters &parameters);
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
	virtual bool reionised_halo (double v, double z) const;
};

/// The Sobacchi13 model of reionisation
class Sobacchi13Reionisation : public Reionisation {
public:
	using Reionisation::Reionisation;
	virtual bool reionised_halo (double v, double z) const;
};

typedef std::shared_ptr<Reionisation> ReionisationPtr;

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