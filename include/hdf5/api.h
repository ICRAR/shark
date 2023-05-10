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
 * Type traits for hdf5 data handling
 */

#ifndef SHARK_HDF5_API_H
#define SHARK_HDF5_API_H

#include <string>
#include <vector>
#include <hdf5.h>

namespace shark {
namespace hdf5 {

class Group;
class DataSet;
class DataSpace;
class Attribute;

class Resource {
public:
    Resource(H5I_type_t expectedType, hid_t handle);
    Resource (const Resource& other);
    Resource& operator= (const Resource& rhs);
    Resource(Resource&&) = default;
    Resource& operator=(Resource&&) = default;
    virtual ~Resource();

    void setComment(const std::string& comment);

    hid_t getHandle() const;

private:
    hid_t handle;

	bool isValid() const;
};

class DataType : public Resource {
public:
    explicit DataType(hid_t handle);
    ~DataType() override;
};

class PredefinedDataType : public DataType {
public:
    static const PredefinedDataType& C_S1();
    static const PredefinedDataType& NATIVE_INT8();
    static const PredefinedDataType& NATIVE_FLOAT();
    static const PredefinedDataType& NATIVE_DOUBLE();
    static const PredefinedDataType& NATIVE_INT();
    static const PredefinedDataType& NATIVE_INT32();
    static const PredefinedDataType& NATIVE_UINT();
    static const PredefinedDataType& NATIVE_UINT32();
    static const PredefinedDataType& NATIVE_INT64();

private:
    explicit PredefinedDataType(hid_t handle);
};

class StringDataType : public DataType {
public:
    explicit StringDataType(size_t size = H5T_VARIABLE);
};

class Location : public Resource {
public:
    // A file, group, dataset
    Location(H5I_type_t expectedType, hid_t handle);

    bool attributeExists(const std::string& name) const;
    Attribute openAttribute(const std::string& name) const;
    Attribute createAttribute(const std::string& name, const DataType& dataType, const DataSpace& dataSpace);
};

class AbstractGroup : public Location {
public:
    AbstractGroup(H5I_type_t expectedType, hid_t handle);

    hsize_t getNumObjs() const;
    std::string getObjnameByIdx(hsize_t idx) const;
    H5G_obj_t getObjTypeByIdx(hsize_t idx) const;
    Group openGroup(const std::string& name) const;
    DataSet openDataSet(const std::string& name) const;
    Group createGroup(const std::string& name) ;
    DataSet createDataSet(const std::string& name, const DataType& dataType, const DataSpace& dataSpace);
};

enum class FileOpenMethod {
	// Open file for reading only
	Read,
	// Open file for reading and writing
	ReadWrite,
	// Create empty file, deleting it if it already exists
	Truncate,
	// Create empty file, failing if it already exists
	ExclusiveWrite,
};

class File : public AbstractGroup {
public:
    File(const std::string& filename, const FileOpenMethod& openMethod);
    ~File() override;

    const std::string& getFileName() const;

private:
    std::string filename;

	static hid_t openOrCreate(const std::string& filename, const FileOpenMethod& openMethod);
};

class Group : public AbstractGroup {
public:
    Group(const AbstractGroup& parent, const std::string& name);
    ~Group() override;

    static Group create(AbstractGroup& parent, const std::string& name);

private:
    explicit Group(hid_t handle);
};

class DataSet : public Location {
public:
    DataSet(const AbstractGroup& file, const std::string& name);
    ~DataSet() override;

    static DataSet create(AbstractGroup& parent, const std::string& name, const DataType& dataType, const DataSpace& dataSpace);

    const std::string& getObjName() const;
    hid_t getDataType() const;
    DataSpace getSpace() const;

    herr_t read(void* buf, hid_t dataType, const DataSpace& memSpace, const DataSpace& fileSpace) const;
    herr_t write(const void* buf, const DataType& memDataType, const DataSpace& memSpace, const DataSpace& fileSpace);

private:
    explicit DataSet(std::string name, hid_t handle);
    std::string objName;
};

enum class DataSpaceType {
    Null = H5S_NULL,
    Scalar = H5S_SCALAR,
    Simple = H5S_SIMPLE,
};

enum class HyperslabSelection {
    // More available
    // See https://docs.hdfgroup.org/hdf5/v1_14/group___h5_s.html#ga6adfdf1b95dc108a65bf66e97d38536d
    Set = H5S_SELECT_SET,
};

class DataSpace : public Resource {
public:
    explicit DataSpace(const DataSet& dataSet);
    ~DataSpace() override;

    static DataSpace create(const DataSpaceType& type);
    static DataSpace create(std::vector<hsize_t> dimensions);

    int getSimpleExtentNdims() const;
    std::vector<hsize_t> getSimpleExtentDims() const;
    std::vector<hsize_t> getSimpleExtentMaxDims() const;
    void selectHyperslab(const HyperslabSelection& op, const std::vector<hsize_t>& start, const std::vector<hsize_t>& stride, const std::vector<hsize_t>& count, const std::vector<hsize_t>& block);

private:
    explicit DataSpace(hid_t handle);
};

class Attribute : public Resource {
public:
    Attribute(const Location& location, const std::string& name);
    ~Attribute() override;

    static Attribute create(Location& location, const std::string& name, const DataType& dataType, const DataSpace& dataSpace);

    template<typename T>
    T read() const {
        hid_t type = H5Aget_type(getHandle());
        T val;
        H5Aread(getHandle(), type, &val);
        return val;
    }

    template<typename T>
    void write(const DataType& dataType, const T& val) {
        H5Awrite(getHandle(), dataType.getHandle(), &val);
    }

private:
    explicit Attribute(hid_t handle);
};

} // namespace hdf5
} // namespace shark

#endif //SHARK_HDF5_API_H
