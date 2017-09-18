//
// Header file for utility routines
//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2017
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

#ifndef SHARK_UTILS
#define SHARK_UTILS

#include <algorithm>
#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace shark {

/**
 * Breaks down a string into substrings delimited by delims
 *
 * Taken from:
 *  http://oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html
 */
std::vector<std::string> tokenize(const std::string &s, const std::string &delims);

/**
 * Changes string `s` to be all lower-case
 * @param s A string
 */
void lower(std::string &s);

/**
 * Changes string `s` to be all upper-case
 * @param s A string
 */
void upper(std::string &s);

/**
 * Trims string `s` on both sides
 * @param s A string
 */
void trim(std::string &s);

/**
 * Performs a binary search and returns the iterator pointing to the element
 * found by the search.
 *
 * Based on the code found on:
 * https://stackoverflow.com/questions/446296/where-can-i-get-a-useful-c-binary-search-algorithm
 *
 * @param begin An iterator
 * @param end Another iterator
 * @param val The search term
 * @param comp A functor to use for equality comparisons.
 *
 * @return An iterator pointing to the element that matches `val`
 */
template<class Iter, class T, typename Comp>
Iter binary_find(Iter begin, Iter end, T val, Comp comp)
{
	// Finds the lower bound in at most log(last - first) + 1 comparisons
	Iter i = std::lower_bound(begin, end, val, comp);

	if (i != end && comp(*i, val)) {
		return i; // found
	}
	return end; // not found
}

template<typename K, typename V>
std::vector<K> get_keys(const std::map<K, V> m)
{
	std::vector<K> keys;
	std::transform(m.begin(), m.end(), std::back_inserter(keys), [](const std::pair<K, V> kv) {
		return kv.first;
	});
	return keys;
}

/**
 * Opens the given file and returns an ifstream for it. If there is a problem
 * while opening the file an exception is throw.
 *
 * @param name The filename
 * @return An ifstream for the opened file
 */
std::ifstream open_file(const std::string &name);

/**
 * Returns \code{true} if the string is empty or corresponds to a line of
 * comment (i.e., starts with the hash (#) character) from a script or
 * configuration file.
 *
 * @param s The string
 * @return Whether the string is empty or is a comment
 */
bool empty_or_comment(const std::string &s);

}  // namespace shark

#endif // SHARK_UTILS