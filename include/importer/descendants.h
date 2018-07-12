//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2017
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
 * Descendant-related class definitions
 */

#ifndef SHARK_IMPORTER_DESCENDANTS
#define SHARK_IMPORTER_DESCENDANTS

#include <string>
#include <vector>

namespace shark {

namespace importer {

/**
 * The type of data that one finds in the descendants file
 */
struct descendants_data_t {
	long halo_id;
	int  halo_snapshot;
	long descendant_id;
	int  descendant_snapshot;
};

class DescendantReader {

public:
	DescendantReader(const std::string &filename);
	virtual ~DescendantReader();
	virtual std::vector<descendants_data_t> read_whole() = 0;

protected:
	std::string filename;
};

class AsciiDescendantReader : public DescendantReader {

public:
	AsciiDescendantReader(const std::string &filename);
	std::vector<descendants_data_t> read_whole() override;

};

class HDF5DescendantReader : public DescendantReader {

public:
	HDF5DescendantReader(const std::string &filename);
	std::vector<descendants_data_t> read_whole() override;

};

}  // namespace importer

}  // namespace shark

#endif // SHARK_IMPORTER_DESCENDANTS