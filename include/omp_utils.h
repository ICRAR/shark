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
 *
 * OpenMP-related utilities. These allow users to easily write "parallel for"
 * constructs without having to deal with the OpenMP pragmas themselves, and
 * ensuring they don't break compatibility with the minimum supported version
 * of OpenMP we ensure with shark.
 *
 * At the moment, we offer compatibility with OpenMP 2.0. This is a pretty low
 * version (5.0 is about to come out as of this writing), but on the one hand
 * we don't really need more at the moment, and it also means we can hit more
 * compilers.
 */

#ifndef SHARK_OMP_UTILS_H_
#define SHARK_OMP_UTILS_H_

#include <type_traits>

#include "config.h"

#ifdef SHARK_OPENMP
#include <omp.h>
#endif // SHARK_OPENMP

namespace shark {

namespace detail {

inline unsigned int thread_index()
{
#ifdef SHARK_OPENMP
	return static_cast<unsigned int>(omp_get_thread_num());
#else
	return 0;
#endif // SHARK_OPENMP
}

template<typename Integer1, typename Integer2>
using omp_int_t = typename std::make_signed<typename std::common_type<Integer1, Integer2>::type>::type;

}  // namespace detail

/**
 * Utility function to call a function over individual members of a container
 * using an OpenMP parallel for with static scheduling. Using this function is
 * simpler than having to repeat this code, and less error prone. If OpenMP
 * support is not present then no parallelization takes place.
 *
 * @param first The first number of the range
 * @param last The last (exclusive) number of the range
 * @param num_threads The number of threads to use for parallelization
 * @param f A callable that takes a single number from the range @p [first,last)
 * and the thread index under which it is being executed
 */
template <typename Integer1, typename Integer2, typename Callable>
void omp_static_for(Integer1 first, Integer2 last, unsigned int num_threads, Callable &&f)
{
	using omp_int = detail::omp_int_t<Integer1, Integer2>;
	using common_type = typename std::common_type<Integer1, Integer2>::type;
	omp_int _first = omp_int(first);
	omp_int _last = omp_int(last);
#ifdef SHARK_OPENMP
	#pragma omp parallel for num_threads(num_threads) schedule(static)
#endif // SHARK_OPENMP
	for (omp_int i = _first; i < _last; i++) {
		f(common_type(i), detail::thread_index());
	}
}

/**
 * Like omp_static_for<Integer, Callable>(Integer, Integer, int, Callable),
 * but takes a container reference instead of two integers, and iterates over
 * the elements of the container. The container's iterator needs to be a random
 * iterator.
 *
 * @param container A container of elements (set, vector, etc.)
 * @param num_threads The number of threads to use for parallelization
 * @param f A callable that takes a single item from the container and the
 * thread index under which it is being executed
 */
template <typename Container, typename Callable>
void omp_static_for(Container &&container, unsigned int num_threads, Callable &&f)
{
	using size_type = typename std::decay<Container>::type::size_type;
	omp_static_for(0, container.size(), num_threads, [&](size_type i, unsigned int thread_num) {
		f(container[i], thread_num);
	});
}

/**
 * Like omp_static_for<Integer, Callable>(Integer, Integer, int, Callable),
 * but using dynamic scheduling with a given chunk size.
 *
 * @param first The first number of the range
 * @param last The last (exclusive) number of the range
 * @param num_threads The number of threads to use for parallelization
 * @param chunk The chunk size for the dynamic scheduling
 * @param f A callable that takes a single number from the range @p [first,last)
 * and the thread index under which it is being executed
 */
template <typename Integer1, typename Integer2, typename Callable>
void omp_dynamic_for(Integer1 first, Integer2 last, unsigned int num_threads, int chunk, Callable &&f)
{
	using omp_int = detail::omp_int_t<Integer1, Integer2>;
	using common_type = typename std::common_type<Integer1, Integer2>::type;
	omp_int _first = omp_int(first);
	omp_int _last = omp_int(last);
#ifdef SHARK_OPENMP
	#pragma omp parallel for num_threads(num_threads) schedule(dynamic, chunk)
#endif // SHARK_OPENMP
	for (omp_int i = _first; i < _last; i++) {
		f(common_type(i), detail::thread_index());
	}
}

/**
 * Like omp_dynamic_for<Integer, Callable>(Integer, Integer, int, int, Callable),
 * but takes a container reference instead of two integers, and iterates over
 * the elements of the container. The container's iterator needs to be a random
 * iterator.
 *
 * @param container A container of elements (set, vector, etc.)
 * @param num_threads The number of threads to use for parallelization
 * @param chunk The chunk size for the dynamic scheduling
 * @param f A callable that takes a single item from the container and the
 * thread index under which it is being executed
 */
template <typename Container, typename Callable>
void omp_dynamic_for(Container &&container, unsigned int num_threads, int chunk, Callable &&f)
{
	using size_type = typename std::decay<Container>::type::size_type;
	omp_dynamic_for(0, container.size(), num_threads, chunk, [&](size_type i, unsigned int thread_num) {
		f(container[i], thread_num);
	});
}

}  // namespace shark

#endif /* SHARK_OMP_UTILS_H_ */