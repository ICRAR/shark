//
// naming convention unit tests
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

#include <cxxtest/TestSuite.h>

#include <string>

#include "naming_convention.h"

using namespace shark;

class TestOptions : public CxxTest::TestSuite {

private:
	void _assert_convention(const std::string &word, naming_convention convention, bool expected)
	{
		TS_ASSERT_EQUALS(follows_convention(word, convention), expected);
	}

	void _assert_convention(const std::vector<std::string> &ok, const std::vector<std::string> &not_ok, naming_convention convention)
	{
		for (auto &w: ok) {
			_assert_convention(w, convention, true);
		}
		for (auto &w: not_ok) {
			_assert_convention(w, convention, false);
		}
	}

public:

	void test_snake_case()
	{
		std::vector<std::string> ok =
		{
			"snake_case", "singleword", "a_b", "a_b_c_d",
			"a0_b1_c1_are_ok_actually"
		};
		std::vector<std::string> not_ok =
		{
			"separated-with-dash", "separated with spaces", "Uppercaseword",
			"ALL_UPPERCASE", "Some_Uppercase", "_starts_with_underscore",
			"ends_with_underscore_", "has_number_0_but_also_has_uppercasE",
			"CamelCaseOfCourseIsWrong"
		};
		_assert_convention(ok, not_ok, naming_convention::SNAKE_CASE);
	}

	void test_camel_case()
	{
		std::vector<std::string> ok =
		{
			"ThisIsGreat", "TBD", "AGNGalaxies", "3DPrinter", "F1Winner"
		};
		std::vector<std::string> not_ok =
		{
			"snake_case", "singleword", "a_b", "a_b_c_d", "ThisIs_Not_SoGreat",
			"a0_b1_c1_are_not_ok_actually", "lowerCamel", "a", "3dprinter"
		};
		_assert_convention(ok, not_ok, naming_convention::CAMEL_CASE);
	}

	void test_lower_camel_case(const std::string &word, bool follows_convention)
	{
		std::vector<std::string> ok =
		{
			"thisIsGreat", "runInfo", "agnGalaxies", "3dPrinter", "f1Winner",
			"singleword", "a"
		};
		std::vector<std::string> not_ok =
		{
			"snake_case", "CamelCase", "singleword", "a_b", "a_b_c_d",
			"a0_b1_c1_are_not_ok_actually", "lowerCamel", "a",
			"thisIs_Not_SoGreat"
		};
		_assert_convention(ok, not_ok, naming_convention::LOWER_CAMEL_CASE);
	}

	void test_unsupported()
	{
		auto _test_unsupported = [](const std::string &word) {
			TS_ASSERT(!follows_convention(word, naming_convention::SNAKE_CASE));
			TS_ASSERT(!follows_convention(word, naming_convention::CAMEL_CASE));
			TS_ASSERT(!follows_convention(word, naming_convention::LOWER_CAMEL_CASE));
			TS_ASSERT(follows_convention(word, naming_convention::NONE));
		};

		_test_unsupported("");
		_test_unsupported("AB_CD");
		_test_unsupported("Hey_You");
		_test_unsupported("hi-there");
		_test_unsupported("with spaces");
	}
};