//
// Components unit tests
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

#include "exceptions.h"
#include "galaxy.h"
#include "halo.h"
#include "merger_tree.h"
#include "subhalo.h"

using namespace shark;

class TestBaryons : public CxxTest::TestSuite {

public:

	void _init_baryons(Baryon &b1, Baryon &b2) {
		b1.mass = 1.;
		b1.mass_metals = 2.;
		b2.mass = 3.;
		b2.mass_metals = 4.;
	}

	void _assert_addition(const Baryon &b) {
		TS_ASSERT_DELTA(b.mass, 4., 1e-8);
		TS_ASSERT_DELTA(b.mass_metals, 6., 1e-8);
	}

	void test_baryons_compound_addition() {
		Baryon b1, b2;
		_init_baryons(b1, b2);
		b2 += b1;
		_assert_addition(b2);
	}

	void test_baryons_addition() {
		Baryon b1, b2;
		_init_baryons(b1, b2);
		Baryon b3 = b2 + b1;
		_assert_addition(b3);
	}

};

class TestSubhalos : public CxxTest::TestSuite
{
private:

	SubhaloPtr make_subhalo(const std::string &types, Subhalo::subhalo_type_t subhalo_type, Galaxy::id_t id=0)
	{
		SubhaloPtr subhalo = std::make_shared<Subhalo>(0, 0);
		subhalo->subhalo_type = subhalo_type;
		for(auto t: types) {
			auto &g = subhalo->emplace_galaxy(id++);
			if (t == 'C') {
				g.galaxy_type = Galaxy::CENTRAL;
			}
			else if (t == '1') {
				g.galaxy_type = Galaxy::TYPE1;
			}
			else if (t == '2') {
				g.galaxy_type = Galaxy::TYPE2;
			}
		}
		TS_ASSERT_EQUALS(types.size(), subhalo->galaxy_count());
		return subhalo;
	}

	template <typename SpecificCheck>
	void _test_valid_galaxy_composition(const std::string &galaxy_types, Subhalo::subhalo_type_t subhalo_type, SpecificCheck specific_check, bool valid)
	{
		auto subhalo = make_subhalo(galaxy_types, subhalo_type);
		if (!valid) {
			TSM_ASSERT_THROWS(galaxy_types, subhalo->check_subhalo_galaxy_composition(), invalid_data);
			TS_ASSERT_THROWS(specific_check(subhalo), invalid_data);
		}
		else {
			subhalo->check_subhalo_galaxy_composition();
			specific_check(subhalo);
		}
	}

	void _test_valid_central_galaxy_composition(const std::string &types, bool valid)
	{
		_test_valid_galaxy_composition(types, Subhalo::CENTRAL, std::mem_fn(&Subhalo::check_central_subhalo_galaxy_composition), valid);
	}

	void _test_valid_satellite_galaxy_composition(const std::string &types, bool valid)
	{
		_test_valid_galaxy_composition(types, Subhalo::SATELLITE, std::mem_fn(&Subhalo::check_satellite_subhalo_galaxy_composition), valid);
	}

public:

	void test_valid_central_galaxy_composition()
	{
		_test_valid_central_galaxy_composition("C", true);
		_test_valid_central_galaxy_composition("1", false);
		_test_valid_central_galaxy_composition("2", false);

		_test_valid_central_galaxy_composition("CC", false);
		_test_valid_central_galaxy_composition("C1", false);
		_test_valid_central_galaxy_composition("C2", true);
		_test_valid_central_galaxy_composition("11", false);
		_test_valid_central_galaxy_composition("12", false);
		_test_valid_central_galaxy_composition("22", false);

		_test_valid_central_galaxy_composition("C22", true);
		_test_valid_central_galaxy_composition("C22222", true);
		_test_valid_central_galaxy_composition("C222221", false);

		_test_valid_central_galaxy_composition("122", false);
		_test_valid_central_galaxy_composition("122222", false);
		_test_valid_central_galaxy_composition("122222C", false);
	}

	void test_valid_satellite_galaxy_composition()
	{
		_test_valid_satellite_galaxy_composition("C", false);
		_test_valid_satellite_galaxy_composition("1", true);
		_test_valid_satellite_galaxy_composition("2", true);

		_test_valid_satellite_galaxy_composition("CC", false);
		_test_valid_satellite_galaxy_composition("C1", false);
		_test_valid_satellite_galaxy_composition("C2", false);
		_test_valid_satellite_galaxy_composition("11", false);
		_test_valid_satellite_galaxy_composition("12", true);
		_test_valid_satellite_galaxy_composition("22", true);

		_test_valid_satellite_galaxy_composition("C22", false);
		_test_valid_satellite_galaxy_composition("C22222", false);
		_test_valid_satellite_galaxy_composition("C222221", false);

		_test_valid_satellite_galaxy_composition("122", true);
		_test_valid_satellite_galaxy_composition("122222", true);
		_test_valid_satellite_galaxy_composition("122222C", false);
	}

	void test_galaxy_finding()
	{
		auto subhalo = make_subhalo("222C222122", Subhalo::SATELLITE, 0);
		TS_ASSERT_EQUALS(3, subhalo->central_galaxy()->id);
		TS_ASSERT_EQUALS(7, subhalo->type1_galaxy()->id);
		TS_ASSERT_EQUALS(8, subhalo->type2_galaxies_count());
	}

	void test_galaxy_movement()
	{
		auto subhalo1 = make_subhalo("C12", Subhalo::CENTRAL, 0);
		auto subhalo2 = make_subhalo("222", Subhalo::SATELLITE, 3);
		subhalo2->transfer_type2galaxies_to(subhalo1);
		TS_ASSERT_EQUALS(6, subhalo1->galaxy_count());
		TS_ASSERT_EQUALS(4, subhalo1->type2_galaxies_count());
		TS_ASSERT_EQUALS(0, subhalo2->type2_galaxies_count());
	}

};

class TestHalos : public CxxTest::TestSuite
{

private:
	template <typename ... Ts>
	HaloPtr make_halo(Ts && ... ts)
	{
		return std::make_shared<Halo>(std::forward<Ts>(ts)...);
	}

	template <typename ... Ts>
	MergerTreePtr make_merger_tree(Ts && ... ts)
	{
		return std::make_shared<MergerTree>(std::forward<Ts>(ts)...);
	}

public:

	void test_roots()
	{

		auto tree = make_merger_tree(1);

		// Merger Tree is:
		//
		// 1 --> 3 --> 5
		//       ^     ^
		// 2 ----|     |
		//             |
		//       4 ----|
		//
		// Roots should be 1, 2 and 4
		auto halo1 = make_halo(1, 1);
		auto halo2 = make_halo(2, 1);
		auto halo3 = make_halo(3, 2);
		auto halo4 = make_halo(4, 2);
		auto halo5 = make_halo(5, 3);
		tree->add_halo(halo5);
		halo5->merger_tree = tree;
		add_parent(halo5, halo3);
		add_parent(halo5, halo4);
		add_parent(halo3, halo2);
		add_parent(halo3, halo1);
		tree->consolidate();

		auto roots = tree->roots();
		TS_ASSERT_EQUALS(roots.size(), 3);
		std::set<Halo::id_t> ids;
		for (auto &root: roots) {
			ids.insert(root->id);
		}
		TS_ASSERT_EQUALS(ids, std::set<Halo::id_t>({1, 2, 4}));

		// 2, 2, and 1 halos for snapshots 1, 2 and 3 respectivelly
		TS_ASSERT_EQUALS(2, tree->halos_at(1).size());
		TS_ASSERT_EQUALS(2, tree->halos_at(2).size());
		TS_ASSERT_EQUALS(1, tree->halos_at(3).size());
	}

};
