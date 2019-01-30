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
 * Simple timer class
 */

#ifndef SHARK_TIMER_H_
#define SHARK_TIMER_H_

#include <ostream>
#include <chrono>

#include "utils.h"

namespace shark {

/**
 * A simple timer class that starts measuring time when created and returns
 * the elapsed time when requested.
 */
class Timer {

public:

	using duration = typename std::chrono::milliseconds::rep;

	/**
	 * Returns the number of milliseconds elapsed since the creation
	 * of the timer
	 *
	 * @return The time elapsed since the creation of the timer, in [ms]
	 */
	inline
	duration get() const {
		return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t0).count();
	}

private:
	std::chrono::steady_clock::time_point t0 {std::chrono::steady_clock::now()};

};

template <typename T>
inline
std::basic_ostream<T> &operator<<(std::basic_ostream<T> &os, const Timer &t) {

	auto time = t.get();
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

}  // namespace shark

#endif // SHARK_TIMER_H_