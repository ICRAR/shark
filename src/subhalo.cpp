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

GalaxyPtr Subhalo::central_galaxy() const
{
	for (auto &galaxy: galaxies){
		if(galaxy->galaxy_type == Galaxy::CENTRAL){
			return galaxy;
		}
	}
	return GalaxyPtr();
}

GalaxyPtr Subhalo::type1_galaxy() const
{
	for (auto &galaxy: galaxies){
		if(galaxy->galaxy_type == Galaxy::TYPE1){
			return galaxy;
		}
	}
	return GalaxyPtr();
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
	for (auto &g: galaxies){
		if (g->galaxy_type == Galaxy::CENTRAL) {
			std::ostringstream os;
			os << "Satellite subhalo " << *this << " has at least one central galaxy";
			throw invalid_data(os.str());
		}
		if (g->galaxy_type == Galaxy::TYPE1) {
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
	for(auto &g: galaxies) {
		if (g->galaxy_type == Galaxy::TYPE1) {
			std::ostringstream os;
			os << "Central subhalo " << *this << " has at least one type 1 galaxy";
			throw invalid_data(os.str());
		}
		else if (g->galaxy_type == Galaxy::CENTRAL) {
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

std::vector<GalaxyPtr> Subhalo::all_type2_galaxies() const
{

	std::vector<GalaxyPtr> all;

	for (auto &galaxy: galaxies){
		if(galaxy->galaxy_type == Galaxy::TYPE2){
			all.push_back(galaxy);
		}
	}

	return all;
}

void Subhalo::copy_galaxies_to(SubhaloPtr &target, const std::vector<GalaxyPtr> &gals) const
{
	target->galaxies.insert(target->galaxies.end(), gals.begin(), gals.end());
}

void Subhalo::transfer_galaxies_to(SubhaloPtr &target)
{
	auto gals_before = target->galaxy_count();
	auto our_gals = galaxies.size();
	if (LOG_ENABLED(trace)) {
		LOG(trace) << "Transferring " << our_gals << " galaxies from " << *this << " to " << target << " (currently " << gals_before << " galaxies)";
	}

	copy_galaxies_to(target, galaxies);
	galaxies.clear();

	assert(gals_before + our_gals == target->galaxy_count());
}

void Subhalo::transfer_type2galaxies_to(SubhaloPtr &target)
{
	auto type2_gals = all_type2_galaxies();
	auto gals_before = target->galaxy_count();
	auto our_gals = type2_gals.size();
	if (LOG_ENABLED(trace)) {
		LOG(trace) << "Transferring " << our_gals << " galaxies from " << *this << " to " << target << " (currently " << gals_before << " galaxies)";
	}

	copy_galaxies_to(target, type2_gals);
	remove_galaxies(type2_gals);

	assert(gals_before + our_gals == target->galaxy_count());
}


void Subhalo::remove_galaxies(const std::vector<GalaxyPtr> &to_remove)
{
	// TODO: Maybe not most efficiently, but it will do for now
	for(auto &galaxy: to_remove) {
		auto it = std::find(galaxies.begin(), galaxies.end(), galaxy);
		if (it == galaxies.end()) {
			LOG(warning) << "Trying to remove galaxy " << galaxy << " which is not in subhalo " << *this << ", ignoring";
			continue;
		}
		if (LOG_ENABLED(debug)) {
			LOG(debug) << "Removing galaxy " << galaxy << " from subhalo " << *this;
		}
		galaxies.erase(it);
	}
}

double Subhalo::total_baryon_mass() const
{
	double mass= 0.0;

	// add subhalo components.
	mass += hot_halo_gas.mass + cold_halo_gas.mass + ejected_galaxy_gas.mass + lost_galaxy_gas.mass;

	for (auto &galaxy: galaxies){
		mass += galaxy->baryon_mass() + galaxy->smbh.mass;
	}

	return mass;
}
}  // namespace shark