//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2018
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

/**
 * @file
 */

#include <regex>
#include <sstream>

#include "exceptions.h"
#include "naming_convention.h"

namespace shark {

static std::regex snake_case_regex("[a-z0-9]+(_[a-z0-9]+)*", std::regex::ECMAScript);
static std::regex camel_case_regex("([0-9]+)?[A-Z][a-z0-9]*([A-Z0-9][a-z0-9]*)*", std::regex::ECMAScript);
static std::regex lower_camel_case_regex("[a-z0-9]+([A-Z][a-z0-9]*)*", std::regex::ECMAScript);

bool follows_convention(const std::string &word, const naming_convention convention)
{
	if (convention == naming_convention::NONE) {
		return true;
	}
	else if (convention == naming_convention::SNAKE_CASE) {
		return std::regex_match(word, snake_case_regex);
	}
	else if (convention == naming_convention::CAMEL_CASE) {
		return std::regex_match(word, camel_case_regex);
	}
	else if (convention == naming_convention::LOWER_CAMEL_CASE) {
		return std::regex_match(word, lower_camel_case_regex);
	}

	std::ostringstream os;
	os << "Unsupported naming convention: " << convention;
	throw invalid_argument(os.str());
}

}  // namespace shark