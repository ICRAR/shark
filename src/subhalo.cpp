//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2019
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

#include <algorithm>
#include <iterator>
#include <sstream>

#include "exceptions.h"
#include "logging.h"
#include "galaxy.h"
#include "subhalo.h"

namespace shark {

SubhaloPtr Subhalo::main() const
{
	for (auto &sub: ascendants) {
		if (sub->main_progenitor) {
			return sub;
		}
	}
	return SubhaloPtr();
}

static
ConstGalaxyPtr pointer_to(const std::vector<Galaxy> &galaxies, Galaxy::galaxy_type_t type)
{
	auto galaxy = std::find_if(galaxies.begin(), galaxies.end(), [type](const Galaxy &galaxy) {
		return galaxy.galaxy_type == type;
	});
	if (galaxy == galaxies.end()) {
		return nullptr;
	}
	return &(*galaxy);
}

ConstGalaxyPtr Subhalo::central_galaxy() const
{
	return pointer_to(galaxies, Galaxy::CENTRAL);
}

GalaxyPtr Subhalo::central_galaxy()
{
	return const_cast<GalaxyPtr>(
		static_cast<const Subhalo &>(*this).central_galaxy());
}

ConstGalaxyPtr Subhalo::type1_galaxy() const
{
	return pointer_to(galaxies, Galaxy::TYPE1);
}

GalaxyPtr Subhalo::type1_galaxy()
{
	return const_cast<GalaxyPtr>(
		static_cast<const Subhalo &>(*this).type1_galaxy());
}

void Subhalo::check_subhalo_galaxy_composition() const
{
	if (subhalo_type == Subhalo::SATELLITE) {
		do_check_satellite_subhalo_galaxy_composition();
	}
	else {
		do_check_central_subhalo_galaxy_composition();
	}
}

void Subhalo::check_central_subhalo_galaxy_composition() const
{
	if (subhalo_type != Subhalo::CENTRAL) {
		std::ostringstream os;
		os << *this << " was expected to be of type central, but has type " << subhalo_type;
		throw invalid_data(os.str());
	}
	do_check_central_subhalo_galaxy_composition();
}

void Subhalo::check_satellite_subhalo_galaxy_composition() const
{
	if (subhalo_type != Subhalo::SATELLITE) {
		std::ostringstream os;
		os << *this << " was expected to be of type satellite, but has type " << subhalo_type;
		throw invalid_data(os.str());
	}
	do_check_satellite_subhalo_galaxy_composition();
}

void Subhalo::do_check_satellite_subhalo_galaxy_composition() const
{
	// If satellite subhalos have one or more galaxies, exactly one must be of
	// type TYPE1, and no TYPE1 galaxies at all

	auto i = 0;
	for (const auto &g: galaxies){
		if (g.galaxy_type == Galaxy::CENTRAL) {
			std::ostringstream os;
			os << "Satellite subhalo " << *this << " has at least one central galaxy";
			throw invalid_data(os.str());
		}
		if (g.galaxy_type == Galaxy::TYPE1) {
			i++;
		}
	}
	if (i > 1) {
		std::ostringstream os;
		os << "Satellite Subhalo " << *this << " has " << i <<" type 1 galaxies (should be <= 1)";
		throw invalid_data(os.str());
	}
}

void Subhalo::do_check_central_subhalo_galaxy_composition() const
{
	// If central subhalos have one or more galaxies, exactly one must be of
	// type CENTRAL, and no TYPE1 galaxies at all

	auto n_central = 0;
	for(const auto &g: galaxies) {
		if (g.galaxy_type == Galaxy::TYPE1) {
			std::ostringstream os;
			os << "Central subhalo " << *this << " has at least one type 1 galaxy";
			throw invalid_data(os.str());
		}
		else if (g.galaxy_type == Galaxy::CENTRAL) {
			n_central++;
		}
	}

	if (n_central == 0 && galaxy_count() > 0) {
		std::ostringstream os;
		os << "Central Subhalo " << *this << " has no central galaxy";
		throw invalid_data(os.str());
	}
	else if (n_central > 1) {
		std::vector<float> mbaryon(galaxies.size());
		std::transform(galaxies.begin(), galaxies.end(), mbaryon.begin(), std::mem_fn(&Galaxy::baryon_mass));
		std::ostringstream os;
		os << "Central Subhalo " << *this << " has " << n_central <<" central galaxies";
		os << "Baryon masses of galaxies: ";
		std::copy(mbaryon.begin(), mbaryon.end(), std::ostream_iterator<float>(os, ", "));
		throw invalid_data(os.str());
	}
}

const Subhalo::type2_galaxies_view Subhalo::type2_galaxies() const
{
	return {galaxies};
}

Subhalo::type2_galaxies_view Subhalo::type2_galaxies()
{
	return {galaxies};
}

void Subhalo::transfer_galaxies_to(SubhaloPtr &target)
{
	auto gals_before = target->galaxy_count();
	auto our_gals = galaxies.size();
	if (LOG_ENABLED(trace)) {
		LOG(trace) << "Transferring " << our_gals << " galaxies from " << *this << " to " << target << " (currently " << gals_before << " galaxies)";
	}
	if (gals_before == 0) {
		std::swap(galaxies, target->galaxies);
	}
	else {
		std::move(galaxies.begin(), galaxies.end(), std::back_inserter(target->galaxies));
		galaxies.clear();
		galaxies.shrink_to_fit();
	}
	assert(galaxies.empty());
	assert(gals_before + our_gals == target->galaxy_count());
}

void Subhalo::transfer_type2galaxies_to(SubhaloPtr &target)
{
	auto our_type2_gals = type2_galaxies_count();
	if (our_type2_gals == 0) {
		return;
	}
	auto gals_before = target->galaxy_count();
	if (LOG_ENABLED(trace)) {
		LOG(trace) << "Transferring " << our_type2_gals << " galaxies from " << *this << " to " << target << " (currently " << gals_before << " galaxies)";
	}
	std::vector<Galaxy::id_t> ids;
	auto type2_gals = type2_galaxies();
	std::transform(type2_gals.begin(), type2_gals.end(), std::back_inserter(ids), [](const Galaxy &g) { return g.id; });
	std::move(type2_gals.begin(), type2_gals.end(), std::back_inserter(target->galaxies));
	remove_galaxies(ids);
	galaxies.shrink_to_fit();
	assert(gals_before + our_type2_gals == target->galaxy_count());
}

void Subhalo::remove_galaxies(const std::vector<Galaxy::id_t> &to_remove)
{
	// TODO: Maybe not most efficiently, but it will do for now
	for(const auto &galaxy_id: to_remove) {
		auto it = std::find_if(galaxies.begin(), galaxies.end(), [galaxy_id](const Galaxy &galaxy) {
			return galaxy.id == galaxy_id;
		});
		if (it == galaxies.end()) {
			LOG(warning) << "Trying to remove galaxy with id=" << galaxy_id << " which is not in subhalo " << *this << ", ignoring";
			continue;
		}
		if (LOG_ENABLED(debug)) {
			LOG(debug) << "Removing galaxy " << *it << " from subhalo " << *this;
		}
		galaxies.erase(it);
	}
}


void Subhalo::transfer_halo_gas_to(SubhaloPtr &target)
{

	target->hot_halo_gas += hot_halo_gas;
	target->cold_halo_gas += cold_halo_gas;

	// After transferring, we make these baryon components zero.
	hot_halo_gas.restore_baryon();
	cold_halo_gas.restore_baryon();

}

double Subhalo::total_baryon_mass() const
{
	double mass= 0.0;

	// add subhalo components.
	mass += hot_halo_gas.mass + cold_halo_gas.mass + ejected_galaxy_gas.mass + lost_galaxy_gas.mass;

	for (const auto &galaxy: galaxies){
		mass += galaxy.baryon_mass() + galaxy.smbh.mass;
	}

	return mass;
}

double Subhalo::inside_subhalo_baryon_mass() const
{
	double mass= 0.0;

	// add subhalo components.
	mass += hot_halo_gas.mass + cold_halo_gas.mass;

	for (const auto &galaxy: galaxies){
		mass += galaxy.baryon_mass() + galaxy.smbh.mass;
	}

	return mass;
}


}  // namespace shark
