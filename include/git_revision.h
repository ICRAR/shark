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
 *
 * Stores the git revision of shark
 */

#ifndef INCLUDE_GIT_REVISION_H_
#define INCLUDE_GIT_REVISION_H_

#include <string>

namespace shark
{

/**
 * Returns the git SHA1 of the current source commit.
 *
 * @return The git SHA1 value for the current source code commit.
 */
std::string get_git_sha1();

} // namespace shark

#endif /* INCLUDE_GIT_REVISION_H_ */