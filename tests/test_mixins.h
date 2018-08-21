//
// Mixins unit tests
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

#include <cxxtest/TestSuite.h>

#include <functional>
#include <random>

#include "mixins.h"
#include "exceptions.h"

using namespace shark;

class TestXYZ : public CxxTest::TestSuite {

public:

	void test_addition()
	{
		xyz<float> x1 {1, 2, 3}, x2 {4, 5, 6};

		auto assert_add_xyz = [](const xyz<float> &x) {
			TS_ASSERT_DELTA(x.x, 5., 1e-8);
			TS_ASSERT_DELTA(x.y, 7., 1e-8);
			TS_ASSERT_DELTA(x.z, 9., 1e-8);
		};
		auto assert_add_scalar = [](const xyz<float> &x) {
			TS_ASSERT_DELTA(x.x, 7., 1e-8);
			TS_ASSERT_DELTA(x.y, 8., 1e-8);
			TS_ASSERT_DELTA(x.z, 9., 1e-8);
		};
		assert_add_xyz(x1 + x2);
		x2 += x1;
		assert_add_xyz(x2);

		assert_add_scalar(x1 + 6);
		x1 += 6;
		assert_add_scalar(x1);
	}

	void test_multiplication()
	{
		xyz<float> x {1, 2, 3};

		auto assert_mul = [](const xyz<float> &x) {
			TS_ASSERT_DELTA(x.x, 3., 1e-8);
			TS_ASSERT_DELTA(x.y, 6., 1e-8);
			TS_ASSERT_DELTA(x.z, 9., 1e-8);
		};
		assert_mul(x * 3);
		x *= 3;
		assert_mul(x);
	}

	void test_division()
	{
		xyz<float> x {1, 2, 3};

		auto assert_div = [](const xyz<float> &x) {
			TS_ASSERT_DELTA(x.x, 1/3., 1e-5);
			TS_ASSERT_DELTA(x.y, 2/3., 1e-5);
			TS_ASSERT_DELTA(x.z, 1., 1e-5);
		};
		assert_div(x / 3);
		x /= 3;
		assert_div(x);
	}

	void test_unit()
	{
		std::default_random_engine engine{std::random_device{}()};
		std::uniform_real_distribution<float> unidist(0, 1000);
		for(int i = 0; i != 100; i++) {
			xyz<float> x {unidist(engine), unidist(engine), unidist(engine)};
			TS_ASSERT_DELTA(x.unit().norm(), 1, 1e-5);
		}
	}

};