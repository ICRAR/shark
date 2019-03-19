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
 * The NFW randon number distribution
 */

#ifndef INCLUDE_NFW_DISTRIBUTION
#define INCLUDE_NFW_DISTRIBUTION

#include <ostream>

namespace shark {

/**
 * A class implementing the Lambert W0 function.
 *
 * Template specializations of this class must implement the () operator, which
 * calculates the Lambert W0 function value for the given input. Both the parameter
 * and the return value must be of type @ref FT.
 *
 * @tparam FT The floating point type supported by this functor
 */
template <typename FT=void>
struct lambert_w0
{
};

/**
 * A random number distribution that follows the NFW profile
 * to distribute its randomly generated numbers.
 *
 * Only floating point types are supported. On top of that, users must provide
 * a function object that can calculate the Lambert W0 function for the given
 * floating point type. Because there is no implementation of this function
 * in the @p std namespace users need to provide their own implementation by
 * specializing the lambert_w0 class.
 *
 * This class satisfies the RandomNumberDistribution C++11 concept, so it can
 * be safely used like any other distribution class from the @p std namespace.
 *
 * This class uses the analytic form of the quantile function described
 * in the work by Robotham & Howlett (2018), arXiv:1805.09550
 */
template <typename FT=double>
class nfw_distribution {

	static_assert(std::is_floating_point<FT>::value, "nfw_distribution must use a floating point type");

public:

	using result_type = FT;

	class param_type {
		const result_type _c;
		const result_type _a;
		const result_type _norm;

	public:

		using distribution_type = nfw_distribution;

		param_type(result_type c=1) : _c(c), _a(1 / c),
		    _norm(std::log((_a + 1) / _a) - 1 / (_a + 1))
		{ }

		result_type c() const { return _c; }
		result_type a() const { return _a; }
		result_type norm() const { return _norm; }

		friend bool operator==(const param_type& x, const param_type& y)
		{ return x._c == y._c; }

		friend
		bool operator!=(const param_type& __x, const param_type& __y)
		{return !(__x == __y);}
	};

private:
	std::uniform_real_distribution<FT> uniform;
	param_type _p;

public:

	// Constructors
	nfw_distribution() : uniform(0, 1), _p() { }
	explicit nfw_distribution(const param_type &p) : uniform(0, 1), _p(p) { }

	// Reset
	void reset() { }

	// Param setter and getter
	param_type param() const { return _p; }
	void param(const param_type &p) { _p = p; }

	// Properties
	result_type c() const { return _p.c(); }
	result_type a() const { return _p.a(); }
	result_type norm() const { return _p.norm(); }

	// Functor operators
	template <typename G>
	result_type operator()(G &g)
	{
		return (*this)(g, _p);
	}

	template <typename G>
	result_type operator()(G &g, const param_type &param)
	{
		auto p = uniform(g);
		auto m = p * param.norm();
		auto wm = - 1.0 / (std::exp(m + 1));
		auto w = lambert_w0<FT>()(wm);
		return param.a() * (std::exp(w + m + 1) - 1);
	}

	// Min/max
	result_type min() const { return 0; }
	result_type max() const { return 1; }

	// Comparison operators
	friend
	bool operator==(const nfw_distribution &lhs, const nfw_distribution &rhs)
	{ return lhs._p == rhs._p; }

	friend
	bool operator!=(const nfw_distribution &lhs, const nfw_distribution &rhs)
	{ return lhs._p != rhs._p; }

};

// Output/Input support
namespace detail {
	template <typename CharT, typename Traits>
	class flag_saver {
	public:
		explicit flag_saver(std::basic_ios<CharT, Traits> &stream) : stream(stream), flags(stream.flags()) {}
		~flag_saver() { stream.flags(flags); }
	private:
		std::basic_ios<CharT, Traits> &stream;
		typename std::basic_ios<CharT, Traits>::fmtflags flags;
	};
} // namespace detail

template <typename CharT, typename Traits, typename RealType>
std::basic_ostream<CharT, Traits> &
operator<<(std::basic_ostream<CharT, Traits> &os, const nfw_distribution<RealType> &x)
{
	detail::flag_saver<CharT, Traits> saver(os);
	os.flags(std::ios_base::dec);
	CharT sp = os.widen(' ');
	return os << x.c();
}

template <typename CharT, typename Traits, typename RealType>
std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits>& is, nfw_distribution<RealType>& x)
{
	using result_type = typename nfw_distribution<RealType>::result_type;
	using param_type = typename nfw_distribution<RealType>::param_type;

	detail::flag_saver<CharT, Traits> saver(is);
	is.flags(std::ios_base::dec | std::ios_base::skipws);
	result_type c;
	is >> c;
	if (!is.fail())
		x.param(param_type(c));
	return is;
}

} // namespace shark

#endif // INCLUDE_NFW_DISTRIBUTION