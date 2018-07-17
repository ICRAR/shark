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
#include <numeric>

#include "components.h"
#include "exceptions.h"
#include "logging.h"

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
	for (auto galaxy: galaxies){
		if(galaxy->galaxy_type == Galaxy::CENTRAL){
			return galaxy;
		}
	}
	return GalaxyPtr();
}

void Subhalo::copy_galaxies_to(SubhaloPtr &target) const
{
	target->galaxies.insert(target->galaxies.end(), galaxies.begin(), galaxies.end());
}

void Subhalo::transfer_galaxies_to(SubhaloPtr &target)
{
	auto gals_before = target->galaxy_count();
	auto our_gals = galaxies.size();
	LOG(trace) << "Transferring " << our_gals << " galaxies from " << *this << " to " << target << " (currently " << gals_before << " galaxies)";

	copy_galaxies_to(target);
	galaxies.clear();

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
		LOG(debug) << "Removing galaxy " << galaxy << " from subhalo " << *this;
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

	return all;
}

void Halo::add_subhalo(const SubhaloPtr &&subhalo)
{
	// Add subhalo mass to halo
	Mvir += subhalo->Mvir;

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
	mass += hot_halo_gas.mass + cold_halo_gas.mass + ejected_galaxy_gas.mass;

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

unsigned long Halo::galaxy_count() const
{
	unsigned long count = 0;
	if (central_subhalo) {
		count = central_subhalo->galaxy_count();
	}

	const auto &sats = satellite_subhalos;
	return std::accumulate(sats.begin(), sats.end(), count, [](unsigned long galaxy_count, const SubhaloPtr &subhalo) {
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
