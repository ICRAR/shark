//
// Options class unit tests
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

#include <cxxtest/TestSuite.h>

#include "options.h"

using namespace shark;

class TestOptions : public CxxTest::TestSuite {

public:

	void _test_parse_option(const std::string optspec, const std::string &exp_name, const std::string &exp_value)
	{
		std::string name, value;
		Options::parse_option(optspec, name, value);
		TS_ASSERT_EQUALS(name, exp_name);
		TS_ASSERT_EQUALS(value, exp_value);
	}

	void test_parse_options() {
		_test_parse_option("key = value", "key", "value");
		_test_parse_option("key= value", "key", "value");
		_test_parse_option("key =value", "key", "value");
		_test_parse_option(" key = value", "key", "value");
		_test_parse_option("key = value ", "key", "value");
		_test_parse_option("key   =    value", "key", "value");
		_test_parse_option("key=value", "key", "value");
		_test_parse_option("a = 1 2 3", "a", "1 2 3");
	}

	void _test_invalid_parse_option(const std::string &optspec) {
		std::string name, value;
		TS_ASSERT_THROWS(Options::parse_option(optspec, name, value), invalid_option &);
	}

	void test_parse_options_invalid() {
		_test_invalid_parse_option("key=");
		_test_invalid_parse_option("key=   ");
		_test_invalid_parse_option("=value");
	}

	template <typename T>
	void _test_load(const std::string &optspec, const std::string &name, T exp_vals)
	{
		T values;
		Options opts;
		opts.add(optspec);
		opts.load(name, values);
		TS_ASSERT_EQUALS(exp_vals, values);
	}

	void test_int_vectors() {
		_test_load<std::vector<int>>("a = 1", "a", {1});
		_test_load<std::vector<int>>("a = 1 2 3", "a", {1, 2, 3});
		_test_load<std::vector<int>>("a = 1 2", "a", {1, 2});
		_test_load<std::vector<int>>("a = 3 2 1", "a", {3, 2, 1});
		_test_load<std::vector<int>>("a = 1 2 3 4", "a", {1, 2, 3, 4});
		_test_load<std::vector<int>>("a = 1-4", "a", {1, 2, 3, 4});
		_test_load<std::vector<int>>("a = 4-1", "a", {4, 3, 2, 1});
		_test_load<std::vector<int>>("a = 4-1 1 1", "a", {4, 3, 2, 1, 1, 1});
	}

	void test_int_sets() {
		_test_load<std::set<int>>("a = 1", "a", {1});
		_test_load<std::set<int>>("a = 1 2 3", "a", {1, 2, 3});
		_test_load<std::set<int>>("a = 1 2", "a", {1, 2});
		_test_load<std::set<int>>("a = 3 2 1", "a", {3, 2, 1});
		_test_load<std::set<int>>("a = 1 2 3 4", "a", {1, 2, 3, 4});
		_test_load<std::set<int>>("a = 1-4", "a", {1, 2, 3, 4});
		_test_load<std::set<int>>("a = 4-1 1 1", "a", {4, 3, 2, 1});
		_test_load<std::set<int>>("a = 4-1 1 1", "a", {1, 2, 3, 4});
	}

	void test_valid_option_names()
	{
		// Just the fact that these run means we are find
		auto _test_valid_option_name = [](const std::string &optspec) {
			Options opts;
			opts.add(optspec + " = value");
		};
		_test_valid_option_name("snake_case");
		_test_valid_option_name("group1.snake_case");
		_test_valid_option_name("group1_name.snake_case");
		_test_valid_option_name("group1_name.3d_properties");
	}

	void test_invalid_option_names()
	{
		auto _test_invalid_option_name = [](const std::string &optspec) {
			Options opt;
			TS_ASSERT_THROWS(opt.add(optspec + " = value"), invalid_option &);
		};
		_test_invalid_option_name("CamelCase");
		_test_invalid_option_name("group.CamelCase");
		_test_invalid_option_name("group.lowerCamelCase");
		_test_invalid_option_name("group.lowerCamelCase");
		_test_invalid_option_name("Group.snake_case");
	}

};