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

/**
 * @file
 */

#include <algorithm>
#include <cassert>
#include <functional>
#include <iterator>
#include <numeric>
#include <sstream>

#include "components.h"
#include "exceptions.h"
#include "logging.h"

namespace shark {

std::vector<HaloPtr> MergerTree::NONE;

SubhaloPtr Subhalo::main() const
{
	for (auto &sub: ascendants) {
		if (sub->main_progenitor) {
			return sub;
		}
	}
	return SubhaloPtr();
}

HaloPtr Halo::main_progenitor() const
{
	auto prog_cen_subh = central_subhalo->main();
	if (prog_cen_subh) {
		return prog_cen_subh->host_halo;
	}

	return HaloPtr();
}

HaloPtr Halo::final_halo() const
{
	auto final_halos = merger_tree->halos_at_last_snapshot();
	assert(final_halos.size() == 1);
	return final_halos[0];
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

std::vector<SubhaloPtr> Halo::all_subhalos() const
{

	std::vector<SubhaloPtr> all;

	if (central_subhalo) {
		all.push_back(central_subhalo);
	}
	all.insert(all.end(), satellite_subhalos.begin(), satellite_subhalos.end());

	// If there are more than one subhalo, then return them ordered by mass in decreasing order.
	if(all.size() > 1){
		std::sort(all.begin(), all.end(), [](const SubhaloPtr &lhs, const SubhaloPtr &rhs) {
			return lhs->Mvir > rhs->Mvir;
		});
	}

	assert(all.size() == satellite_subhalos.size() + (central_subhalo ? 1 : 0));
	return all;
}

void add_parent(const HaloPtr &halo, const HaloPtr &parent)
{
	auto result = halo->ascendants.insert(parent);
	auto halos_linked = std::get<1>(result);

	// Fail if a halo has more than one descendant
	if (parent->descendant && parent->descendant->id != halo->id) {
		std::ostringstream os;
		os << parent << " already has a descendant " << parent->descendant;
		os << " but " << halo << " is claiming to be its descendant as well";
		throw invalid_data(os.str());
	}
	parent->descendant = halo;

	// Link this halo to merger tree and back
	if (!halo->merger_tree) {
		std::ostringstream os;
		os << "Descendant " << halo << " has no MergerTree associated to it";
		throw invalid_data(os.str());
	}
	parent->merger_tree = halo->merger_tree;
	if (halos_linked) {
		parent->merger_tree->add_halo(parent);
	}
}

void Halo::add_subhalo(SubhaloPtr &&subhalo)
{
	// Add subhalo mass to halo
	Mvir += subhalo->Mvir;
	Mgas += subhalo->Mgas;

	// Assign subhalo to proper member
	if (subhalo->subhalo_type == Subhalo::CENTRAL) {
		central_subhalo = std::move(subhalo);
	}
	else {
		satellite_subhalos.emplace_back(std::move(subhalo));
	}
}

void Halo::remove_subhalo(const SubhaloPtr &subhalo)
{
	if (subhalo == central_subhalo) {
		central_subhalo.reset();
		return;
	}

	auto it = std::find(satellite_subhalos.begin(), satellite_subhalos.end(), subhalo);
	if (it == satellite_subhalos.end()) {
		throw subhalo_not_found("subhalo not in satellites", subhalo->id);
	}
	satellite_subhalos.erase(it);
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

double Halo::total_baryon_mass() const
{
	double mass= 0.0;

	for (auto &subhalo: all_subhalos()){
		mass += subhalo->total_baryon_mass();
	}

	return mass;
}

galaxies_size_type Halo::galaxy_count() const
{
	galaxies_size_type count = 0;
	if (central_subhalo) {
		count = central_subhalo->galaxy_count();
	}

	const auto &sats = satellite_subhalos;
	return std::accumulate(sats.begin(), sats.end(), count, [](galaxies_size_type galaxy_count, const SubhaloPtr &subhalo) {
		return galaxy_count + subhalo->galaxy_count();
	});
}

std::vector<double> TotalBaryon::get_masses (const std::vector<BaryonBase> &B) const
{
	std::vector<double> masses(B.size());
	std::transform(B.begin(), B.end(), masses.begin(), [](const BaryonBase &b) {
		return b.mass;
	});
	return masses;
}

std::vector<double> TotalBaryon::get_metals (const std::vector<BaryonBase> &B) const
{
	std::vector<double> masses(B.size());
	std::transform(B.begin(), B.end(), masses.begin(), [](const BaryonBase &b) {
		return b.mass_metals;
	});
	return masses;
}

}  // namespace shark
