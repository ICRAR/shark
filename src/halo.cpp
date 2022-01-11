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
#include <cassert>
#include <numeric>
#include <sstream>
#include <vector>

#include "exceptions.h"
#include "halo.h"
#include "merger_tree.h"
#include "subhalo.h"

namespace shark {

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
	return *final_halos.begin();
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

double Halo::total_baryon_mass() const
{
	double mass= 0.0;

	for (auto &subhalo: all_subhalos()){
		mass += subhalo->total_baryon_mass();
	}

	return mass;
}

double Halo::inside_halo_baryon_mass() const
{
	double mass= 0.0;

	for (auto &subhalo: all_subhalos()){
		mass += subhalo->inside_subhalo_baryon_mass();
	}

	return mass;
}

double Halo::total_mass_ascendants() const
{
	double mass= 0.0;

	for (auto &halo: ascendants){
		mass += halo->Mvir;
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

void add_parent(const HaloPtr &halo, const HaloPtr &parent)
{
	auto parent_in_ascendants = std::find(halo->ascendants.begin(), halo->ascendants.end(), parent);
	auto parent_previously_linked = parent_in_ascendants != halo->ascendants.end();
	if (!parent_previously_linked) {
		halo->ascendants.push_back(parent);
	}

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
	if (!parent_previously_linked) {
		parent->merger_tree->add_halo(parent);
	}
}

}  // namespace shark
