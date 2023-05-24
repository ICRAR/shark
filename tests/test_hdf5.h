//
// HDF5 simple unit tests
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

#include <utility>

#include <cxxtest/TestSuite.h>

#include <boost/filesystem.hpp>
#include "hdf5/io/reader.h"
#include "hdf5/io/writer.h"

using namespace shark;
namespace fs = boost::filesystem;

class TestHDF5 : public CxxTest::TestSuite {

private:
	template<typename ... Ts>
	hdf5::Writer get_writer(Ts&& ...args) {
		return hdf5::Writer("test.hdf5", std::forward<Ts>(args)...);
	}

	hdf5::Reader get_reader() {
		return hdf5::Reader("test.hdf5");
	}

	void remove_file() {
		fs::path path("test.hdf5");
		if (fs::exists(path)) {
			fs::remove(path);
		}
	}

	template<typename T>
	void _assert_invalid_names(T assert_func) {
		assert_func("MyGroup", naming_convention::SNAKE_CASE);
		assert_func("myGroup", naming_convention::SNAKE_CASE);
		assert_func("3dPrinter", naming_convention::SNAKE_CASE);
		assert_func("3DPrinter", naming_convention::SNAKE_CASE);
		assert_func("snake_case_name", naming_convention::CAMEL_CASE);
		assert_func("comment", naming_convention::CAMEL_CASE);
		assert_func("3dprinter", naming_convention::CAMEL_CASE);
	}

public:

	virtual void tearDown() {
		remove_file();
	}

	void test_write_dataset_scalars() {
		// Reference data
		int an_integer = 1;
		float a_float = 1;
		double a_double = 1;

		// Write first
		{
			auto writer = get_writer();
			writer.write_dataset("integer", an_integer);
			writer.write_dataset("float", a_float);
			writer.write_dataset("double", a_double);
		}

		// Read and verify
		auto reader = get_reader();
		auto hdf5_integer = reader.read_dataset<int>("integer");
		auto hdf5_float = reader.read_dataset<float>("float");
		auto hdf5_double = reader.read_dataset<double>("double");

		TS_ASSERT_EQUALS(an_integer, hdf5_integer);
		TS_ASSERT_EQUALS(a_float, hdf5_float);
		TS_ASSERT_EQUALS(a_double, hdf5_double);
	}

	void test_write_dataset_vectors() {
		// Reference data
		std::vector<int> integers{1, 2, 3, 4};
		std::vector<float> floats{1, 2, 3, 4};
		std::vector<double> doubles{1, 2, 3, 4};

		// Write first
		{
			auto writer = get_writer();
			writer.write_dataset("integers", integers);
			writer.write_dataset("floats", floats);
			writer.write_dataset("doubles", doubles);
		}

		// Read and verify
		auto reader = get_reader();
		auto hdf5_integers = reader.read_dataset_v<int>("integers");
		auto hdf5_floats = reader.read_dataset_v<float>("floats");
		auto hdf5_doubles = reader.read_dataset_v<double>("doubles");

		TS_ASSERT_EQUALS(integers, hdf5_integers);
		TS_ASSERT_EQUALS(floats, hdf5_floats);
		TS_ASSERT_EQUALS(doubles, hdf5_doubles);
	}

	void test_wrong_attribute_writes() {
		// Single-named attributes are not supported
		auto writer = get_writer();
		TS_ASSERT_THROWS(writer.write_attribute("/attr_name", 1), invalid_argument &);
		TS_ASSERT_THROWS(writer.write_attribute("/attr_name", std::string("1")), invalid_argument &);
	}

	template<typename T>
	void _test_attribute_writes(const T& val) {
		// Write/read an attribute both in a dataset and in a group
		{
			auto writer = get_writer();
			writer.write_dataset("/group/integers", std::vector<int>{1, 2, 3, 4});
			writer.write_attribute("/group/attribute", val);
			writer.write_attribute("/group/integers/attribute", val);
		}

		auto reader = get_reader();
		TS_ASSERT_EQUALS(reader.read_attribute<T>("/group/attribute"), val);
		TS_ASSERT_EQUALS(reader.read_attribute<T>("/group/integers/attribute"), val);
	}

	void test_attribute_writes() {
		_test_attribute_writes(1);
		_test_attribute_writes(1.f);
		_test_attribute_writes(1.);
		_test_attribute_writes(std::string("abc"));
	}

	void test_invalid_group_names() {
		_assert_invalid_names([this](const std::string& name, naming_convention convention) {
			auto writer = get_writer(false, convention, naming_convention::SNAKE_CASE, naming_convention::SNAKE_CASE);
			auto goodname = std::string(convention == naming_convention::CAMEL_CASE ? "Group" : "group");
			auto basename = std::string("/") + name;
			auto repeated = basename + std::string("/") + name;
			auto tailname = goodname + "/" + goodname + "/" + name;
			TS_ASSERT_THROWS(writer.write_attribute(basename + "/my_attribute", 1), invalid_argument &);
			TS_ASSERT_THROWS(writer.write_attribute(repeated + "/my_attribute", 1), invalid_argument &);
			TS_ASSERT_THROWS(writer.write_attribute(tailname + "/my_attribute", 1), invalid_argument &);
			writer.close();
			remove_file();
		});
	}

	void test_invalid_dataset_names() {
		_assert_invalid_names([this](const std::string& name, naming_convention convention) {
			auto writer = get_writer(false, naming_convention::SNAKE_CASE, convention, naming_convention::SNAKE_CASE);
			TS_ASSERT_THROWS(writer.write_dataset("/" + name, std::vector<int>{1, 2, 3, 4}), invalid_argument &);
			TS_ASSERT_THROWS(writer.write_dataset("/group/" + name, std::vector<int>{1, 2, 3, 4}), invalid_argument &);
			TS_ASSERT_THROWS(writer.write_dataset("/group1/group2/" + name, std::vector<int>{1, 2, 3, 4}),
			                 invalid_argument &);
			writer.close();
			remove_file();
		});
	}

	void test_invalid_attribute_names() {
		_assert_invalid_names([this](const std::string& name, naming_convention convention) {
			auto writer = get_writer(false, naming_convention::SNAKE_CASE, naming_convention::SNAKE_CASE, convention);
			writer.write_dataset("/group/integers", std::vector<int>{1, 2, 3, 4});
			TS_ASSERT_THROWS(writer.write_attribute(std::string("/group/") + name, 1), invalid_argument &);
			TS_ASSERT_THROWS(writer.write_attribute(std::string("/group/integers") + name, 1), invalid_argument &);
			writer.close();
			remove_file();
		});
	}


};