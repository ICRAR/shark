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

/**
 * @file
 *
 * Header file for utility routines
 */

#ifndef SHARK_UTILS
#define SHARK_UTILS

#include <algorithm>
#include <chrono>
#include <fstream>
#include <functional>
#include <iomanip>
#include <map>
#include <numeric>
#include <ostream>
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
 * Like lower(std::string &), but returns the lowered string instead
 *
 * @param s A string
 * @return the lower-cased string
 */
std::string lower(const std::string &s);

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


namespace detail {

	template <int N, typename T>
	struct _fixed {
		T _val;
	};

	template <typename T, int N, typename VT>
	inline
	std::basic_ostream<T> &operator<<(std::basic_ostream<T> &os, detail::_fixed<N, VT> v)
	{
		os << std::setprecision(N) << std::fixed << v._val;
		return os;
	}

} // namespace detail

///
/// Sent to a stream object, this manipulator will print the given value with a
/// precision of N decimal places.
///
/// @param v The value to send to the stream
///
template <int N, typename T>
inline
detail::_fixed<N, T> fixed(T v) {
	return {v};
}

namespace detail {

	struct _memory_amount {
		std::size_t _val;
	};

	struct _nanoseconds_amount {
		std::chrono::nanoseconds::rep _val;
	};

	template <typename T>
	inline
	std::basic_ostream<T> &operator<<(std::basic_ostream<T> &os, const detail::_memory_amount &m)
	{

		if (m._val < 1024) {
			os << m._val << " [B]";
			return os;
		}

		float v = m._val / 1024.;
		const char *suffix = " [KB]";

		if (v > 1024) {
			v /= 1024;
			suffix = " [MB]";
		}
		if (v > 1024) {
			v /= 1024;
			suffix = " [GB]";
		}
		if (v > 1024) {
			v /= 1024;
			suffix = " [TB]";
		}
		if (v > 1024) {
			v /= 1024;
			suffix = " [PB]";
		}
		if (v > 1024) {
			v /= 1024;
			suffix = " [EB]";
		}
		// that should be enough...

		os << fixed<3>(v) << suffix;
		return os;
	}

	template <typename T>
	inline
	std::basic_ostream<T> &operator<<(std::basic_ostream<T> &os, const detail::_nanoseconds_amount &t)
	{
		auto time = t._val;
		if (time < 1000) {
			os << time << " [ns]";
			return os;
		}

		time /= 1000;
		if (time < 1000) {
			os << time << " [us]";
			return os;
		}

		time /= 1000;
		if (time < 1000) {
			os << time << " [ms]";
			return os;
		}

		float ftime = time / 1000.f;
		const char *prefix = " [s]";
		if (ftime > 60) {
			ftime /= 60;
			prefix = " [min]";
			if (ftime > 60) {
				ftime /= 60;
				prefix = " [h]";
				if (ftime > 24) {
					ftime /= 24;
					prefix = " [d]";
				}
			}
		}
		// that should be enough...

		os << fixed<3>(ftime) << prefix;
		return os;
	}

} // namespace detail

///
/// Sent to a stream object, this manipulator will print the given amount of
/// memory using the correct suffix and 3 decimal places.
///
/// @param v The value to send to the stream
///
inline
detail::_memory_amount memory_amount(std::size_t amount) {
	return {amount};
}

///
/// Sent to a stream object, this manipulator will print the given amount of
/// nanoseconds using the correct suffix and 3 decimal places.
///
/// @param v The value to send to the stream
///
inline
detail::_nanoseconds_amount ns_time(std::chrono::nanoseconds::rep amount) {
	return {amount};
}

/// Returns the name of the computer executing this program
std::string gethostname();

/// A class template for deleters that use a function to delete objects
template<typename T, void (*F)(T *)>
class deleter {
public:
	void operator()(T *x)
	{ F(x); }
};

/// Sum all elements in a collection of type T and return the result
template<typename Container>
typename Container::value_type sum(Container &c)
{
	using value_type = typename Container::value_type;
	return std::accumulate(c.begin(), c.end(), value_type{}, std::plus<value_type>{});
}

/// Returns the maximum amount of memory used by this process
std::size_t peak_rss();

}  // namespace shark

#endif // SHARK_UTILS