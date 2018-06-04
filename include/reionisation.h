/*
 * reionisation.h
 *
 *  Created on: 14Jun.,2017
 *      Author: clagos
 */

#ifndef SHARK_REIONISATION_H_
#define SHARK_REIONISATION_H_

#include <memory>
#include <utility>

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


class Reionisation {

public:
	Reionisation(const ReionisationParameters &parameters);

	bool reionised_halo (double v, double z);

private:

	ReionisationParameters parameters;

};

typedef std::shared_ptr<Reionisation> ReionisationPtr;

template <typename ...Ts>
ReionisationPtr make_reionisation(Ts&&...ts)
{
	return std::make_shared<Reionisation>(std::forward<Ts>(ts)...);
}

}//end namespace shark


#endif //SHARK_REIONISATION_H_