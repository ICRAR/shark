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
 *
 * Functions to retrieve static data files used by shark
 */

#ifndef SHARK_DATA_H_
#define SHARK_DATA_H_

#include <string>

namespace shark
{

/**
 * Given a relative filename @a fname, return the actual path in the filesystem
 * where this file can be found. @a fname must be relative to the location of
 * the @p data directory of the shark source code repository, even though the
 * actual filename returned by this routine does not belong to the repository
 * (e.g., if shark is running from a globally installed location).
 *
 * @param fname The file to find
 * @return The actual path to the file on disk
 */
std::string get_static_data_filepath(const std::string &fname);

}  // namespace shark



#endif /* SHARK_DATA_H_ */
