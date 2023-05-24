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
 * HDF5-related traits
 */

#include <stdexcept>
#include <utility>
#include <sstream>
#include "hdf5/api.h"
#include "logging.h"

namespace shark {
namespace hdf5 {

Resource::Resource(H5I_type_t expectedType, hid_t handle) : handle(handle) {
	if (!isValid()) {
		throw std::runtime_error("Invalid handle");
	}

	auto type = H5Iget_type(handle);
	if (type == H5I_BADID) {
		throw std::runtime_error("Unable to determine resource type");
	} else if (type != expectedType) {
		std::ostringstream os;
		os << "Expected resource of type " << expectedType << ", got " << type;
		throw std::runtime_error(os.str());
	}
}

Resource::Resource(const Resource& other) : handle(other.handle) {
	H5Iinc_ref(handle);
}

Resource& Resource::operator=(const Resource& rhs) {
	H5Idec_ref(handle);
	handle = rhs.handle;
	H5Iinc_ref(handle);
	return *this;
}

// Needed to satisfy the linker, but we want to make overriding required for subclasses
// As a result we don't put this in the class definition
// See here for more details: https://stackoverflow.com/a/11437551
Resource::~Resource() = default;

bool Resource::isValid() const {
	auto isValid = H5Iis_valid(handle);
	if (isValid < 0) {
		throw std::runtime_error("Error in H5Iis_valid");
	}

	return isValid > 0;
}

void Resource::setComment(const std::string& comment) {
	H5Oset_comment(getHandle(), comment.c_str());
}

hid_t Resource::getHandle() const {
	return handle;
}

DataType::DataType(hid_t handle) : Resource(H5I_DATATYPE, handle) {
	// TODO check actually is datatype
}

DataType::~DataType() {
	H5Tclose(getHandle());
}

// Always copy the predefined type so that the deconstructor will succeed!
// (H5Tclose() will fail for the predefined types as they are considered immutable)
PredefinedDataType::PredefinedDataType(hid_t handle) : DataType(H5Tcopy(handle)) {}

const PredefinedDataType& PredefinedDataType::C_S1() {
	static const PredefinedDataType t(H5T_C_S1);
	return t;
}

const PredefinedDataType& PredefinedDataType::NATIVE_INT8() {
	static const PredefinedDataType t(H5T_NATIVE_INT8);
	return t;
}

const PredefinedDataType& PredefinedDataType::NATIVE_FLOAT() {
	static const PredefinedDataType t(H5T_NATIVE_FLOAT);
	return t;
}

const PredefinedDataType& PredefinedDataType::NATIVE_DOUBLE() {
	static const PredefinedDataType t(H5T_NATIVE_DOUBLE);
	return t;
}

const PredefinedDataType& PredefinedDataType::NATIVE_INT() {
	static const PredefinedDataType t(H5T_NATIVE_INT);
	return t;
}

const PredefinedDataType& PredefinedDataType::NATIVE_INT32() {
	static const PredefinedDataType t(H5T_NATIVE_INT32);
	return t;
}

const PredefinedDataType& PredefinedDataType::NATIVE_UINT() {
	static const PredefinedDataType t(H5T_NATIVE_UINT);
	return t;
}

const PredefinedDataType& PredefinedDataType::NATIVE_UINT32() {
	static const PredefinedDataType t(H5T_NATIVE_UINT32);
	return t;
}

const PredefinedDataType& PredefinedDataType::NATIVE_INT64() {
	static const PredefinedDataType t(H5T_NATIVE_INT64);
	return t;
}

StringDataType::StringDataType(size_t size) : DataType(H5Tcopy(PredefinedDataType::C_S1().getHandle())) {
	H5Tset_size(getHandle(), size);
}

Location::Location(H5I_type_t expectedType, hid_t handle) : Resource(expectedType, handle) {}

bool Location::attributeExists(const std::string& name) const {
	auto exists = H5Aexists(getHandle(), name.c_str());
	if (exists < 0) {
		throw std::runtime_error("Error on H5Aexists");
	}
	return exists > 0;
}

Attribute Location::openAttribute(const std::string& name) const {
	return Attribute(*this, name);
}

Attribute Location::createAttribute(const std::string& name, const DataType& dataType, const DataSpace& dataSpace) {
	return Attribute::create(*this, name, dataType, dataSpace);
}

AbstractGroup::AbstractGroup(H5I_type_t expectedType, hid_t handle) : Location(expectedType, handle) {}

hsize_t AbstractGroup::getNumObjs() const {
	hsize_t num;
	H5Gget_num_objs(getHandle(), &num);
	return num;
}

std::string AbstractGroup::getObjnameByIdx(hsize_t idx) const {
	auto size = H5Gget_objname_by_idx(getHandle(), idx, nullptr, 0);

	// Ensure there's enough space for the null terminator that the C API will write
	std::string name;
	name.resize(size + 1);

	H5Gget_objname_by_idx(getHandle(), idx, &name[0], name.size());

	// Now resize back down to the actual size
	name.resize(size);

	return name;
}

H5G_obj_t AbstractGroup::getObjTypeByIdx(hsize_t idx) const {
	return H5Gget_objtype_by_idx(getHandle(), idx);
}

Group AbstractGroup::openGroup(const std::string& name) const {
	return Group(*this, name);
}

DataSet AbstractGroup::openDataSet(const std::string& name) const {
	return DataSet(*this, name);
}

Group AbstractGroup::createGroup(const std::string& name) {
	return Group::create(*this, name);
}

DataSet AbstractGroup::createDataSet(const std::string& name, const DataType& dataType, const DataSpace& dataSpace) {
	return DataSet::create(*this, name, dataType, dataSpace);
}

File::File(const std::string& filename, const FileOpenMethod& openMethod) :
		AbstractGroup(H5I_FILE, openOrCreate(filename, openMethod)),
		filename(filename) {
}

File::~File() {
	H5Fclose(getHandle());
}

hid_t File::openOrCreate(const std::string& filename, const FileOpenMethod& openMethod) {
	switch (openMethod) {
		case FileOpenMethod::Read:
			return H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		case FileOpenMethod::ReadWrite:
			return H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		case FileOpenMethod::Truncate:
			return H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		case FileOpenMethod::ExclusiveWrite:
			return H5Fcreate(filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
		default:
			throw std::runtime_error("Unknown open type");
	}
}

const std::string& File::getFileName() const {
	return filename;
}

Group::Group(hid_t handle) : AbstractGroup(H5I_GROUP, handle) {}

Group::Group(const AbstractGroup& parent, const std::string& name) :
		Group(H5Gopen(parent.getHandle(), name.c_str(), H5P_DEFAULT)) {
}

Group::~Group() {
	H5Gclose(getHandle());
}

Group Group::create(AbstractGroup& parent, const std::string& name) {
	return Group(H5Gcreate(parent.getHandle(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
}

DataSet::DataSet(std::string name, hid_t handle) : Location(H5I_DATASET, handle), objName(std::move(name)) {}

DataSet::DataSet(const AbstractGroup& file, const std::string& name) :
		Location(H5I_DATASET, H5Dopen(file.getHandle(), name.c_str(), H5P_DEFAULT)),
		objName(name) {
}

DataSet::~DataSet() {
	H5Dclose(getHandle());
}

DataSet DataSet::create(shark::hdf5::AbstractGroup& parent, const std::string& name, const DataType& dataType,
                        const shark::hdf5::DataSpace& dataSpace) {
	return DataSet(name,
	               H5Dcreate(parent.getHandle(), name.c_str(), dataType.getHandle(), dataSpace.getHandle(), H5P_DEFAULT,
	                         H5P_DEFAULT, H5P_DEFAULT));
}

DataSpace DataSet::getSpace() const {
	return DataSpace(*this);
}

const std::string& DataSet::getObjName() const {
	return objName;
}

hid_t DataSet::getDataType() const {
	return H5Dget_type(getHandle());
}

herr_t DataSet::read(void *buf, hid_t dataType, const DataSpace& memSpace, const DataSpace& fileSpace) const {
	return H5Dread(getHandle(), dataType, memSpace.getHandle(), fileSpace.getHandle(), H5P_DEFAULT, buf);
}

herr_t
DataSet::write(const void *buf, const DataType& memDataType, const DataSpace& memSpace, const DataSpace& fileSpace) {
	return H5Dwrite(getHandle(), memDataType.getHandle(), memSpace.getHandle(), fileSpace.getHandle(), H5P_DEFAULT,
	                buf);
}

DataSpace::DataSpace(hid_t handle) : Resource(H5I_DATASPACE, handle) {}

DataSpace::DataSpace(const DataSet& dataSet) : Resource(H5I_DATASPACE, H5Dget_space(dataSet.getHandle())) {
}

DataSpace::~DataSpace() {
	H5Sclose(getHandle());
}

DataSpace DataSpace::create(const DataSpaceType& type) {
	return DataSpace(H5Screate(static_cast<H5S_class_t>(type)));
}

DataSpace DataSpace::create(std::vector<hsize_t> dimensions) {
	auto rank = static_cast<int>(dimensions.size());
	return DataSpace(H5Screate_simple(rank, dimensions.data(), nullptr));
}

int DataSpace::getSimpleExtentNdims() const {
	return H5Sget_simple_extent_ndims(getHandle());
}

std::vector<hsize_t> DataSpace::getSimpleExtentDims() const {
	auto ndims = getSimpleExtentNdims();
	std::vector<hsize_t> dim_sizes(ndims);
	H5Sget_simple_extent_dims(getHandle(), dim_sizes.data(), nullptr);
	return dim_sizes;
}

std::vector<hsize_t> DataSpace::getSimpleExtentMaxDims() const {
	auto ndims = getSimpleExtentNdims();
	std::vector<hsize_t> dim_sizes(ndims);
	std::vector<hsize_t> max_dim_sizes(ndims);
	H5Sget_simple_extent_dims(getHandle(), dim_sizes.data(), max_dim_sizes.data());
	return max_dim_sizes;
}

void DataSpace::selectHyperslab(const HyperslabSelection& op, const std::vector<hsize_t>& start,
                                const std::vector<hsize_t>& stride, const std::vector<hsize_t>& count,
                                const std::vector<hsize_t>& block) {
	// Assert all same length & same as rank
	H5Sselect_hyperslab(getHandle(), static_cast<H5S_seloper_t>(op), start.data(), stride.data(), count.data(),
	                    block.data());
}

Attribute::Attribute(hid_t handle) : Resource(H5I_ATTR, handle) {}

Attribute::Attribute(const Location& group, const std::string& name) :
		Attribute(H5Aopen(group.getHandle(), name.c_str(), H5P_DEFAULT)) {
}

Attribute::~Attribute() {
	H5Aclose(getHandle());
}

Attribute Attribute::create(Location& location, const std::string& name, const DataType& dataType,
                            const DataSpace& dataSpace) {
	return Attribute(
			H5Acreate(location.getHandle(), name.c_str(), dataType.getHandle(), dataSpace.getHandle(), H5P_DEFAULT,
			          H5P_DEFAULT));
}

template<>
std::string Attribute::read<std::string>() const {
	H5A_info_t info;
	H5Aget_info(getHandle(), &info);
	std::string val;
	val.resize(info.data_size); // Size includes null-terminator
	auto type = H5Aget_type(getHandle());
	H5Aread(getHandle(), type, &val[0]);
	val.resize(info.data_size - 1); // Trim included null-terminator
	return val;
}

template<>
void Attribute::write<std::string>(const DataType& dataType, const std::string& val) {
	H5Awrite(getHandle(), dataType.getHandle(), val.data());
}

} // namespace hdf5
} // namespace shark