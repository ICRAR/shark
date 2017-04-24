//
// Implementation of the IOBase class for HDF5 handling
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

#include <sstream>
#include <string>
#include <stdexcept>

#include "hdf5/iobase.h"

using namespace std;

namespace shark {

namespace hdf5 {

IOBase::IOBase(const string &filename, unsigned int flags) :
	hdf5_file(filename, flags)
{
	// no-op
}

IOBase::~IOBase()
{
	hdf5_file.close();
}

const string IOBase::get_filename() const
{
	return hdf5_file.getFileName();
}

H5::DataSpace IOBase::get_1d_dataspace(const H5::DataSet &dataset) const {
	H5::DataSpace space = dataset.getSpace();
	int ndims = space.getSimpleExtentNdims();
	if ( ndims != 1 ) {
		ostringstream os;
		os << "More than one dimension found in dataset " << dataset.getObjName();
		throw runtime_error(os.str());
	}
	return space;
}

hsize_t IOBase::get_1d_dimsize(const H5::DataSpace &space) const {
	hsize_t dim_size;
	space.getSimpleExtentDims(&dim_size, NULL);
	return dim_size;
}

}  // namespace hdf5

}  // namespace shark