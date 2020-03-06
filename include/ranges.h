//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2019
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

/**
 * @file
 *
 * Utility classes to easily iterate over subranges of a container.
 */

#ifndef SHARK_RANGES_H_
#define SHARK_RANGES_H_

#include <cassert>
#include <iterator>
#include <type_traits>
#include <tuple>

namespace shark {

/**
 * A range-like filter over a range-like object.
 *
 * This class takes a range object (something with a begin() and an end()
 * method) and a unary predicate function object and allows easy iteration over
 * the filtered elements of the range for which the predicate is true. Iteration
 * is provided via begin() and end() iterators, and therefore this class is also
 * a range.
 */
template <typename Range, typename UnaryPredicate>
class range_filter {

public:

	/**
	 * Iterator class providing filtering vased on the UnaryPredicate.
	 * It is templated on its constness so it can offer a const version of its
	 * dereference operator when required.
	 */
	template <bool is_const>
	class _iterator {

	public:
		using iterator_category = std::input_iterator_tag;
		using value_type = typename Range::value_type;
		using difference_type = std::ptrdiff_t;
		using reference = typename std::conditional<is_const, const value_type &, value_type &>::type;
		using pointer = typename std::conditional<is_const, const value_type *, value_type *>::type;

	private:
		typename Range::const_iterator m_pos;
		typename Range::const_iterator m_end;
		const range_filter &m_range_filter;

	public:
		_iterator(const range_filter &t_range_filter, bool is_end)
		 : m_pos(is_end ? t_range_filter.m_range.end() : t_range_filter.m_range.begin()),
		   m_end(t_range_filter.m_range.end()), m_range_filter(t_range_filter)
		{
			if (!is_end) {
				advance();
			}
		}

		bool operator==(const _iterator &other)
		{
			return m_pos == other.m_pos;
		}

		bool operator!=(const _iterator &other)
		{
			return !(*this == other);
		}

		_iterator &operator++()
		{
			m_pos++;
			advance();
			return *this;
		}

		template <bool _is_const=is_const>
		typename std::enable_if<!_is_const, reference>::type
		operator*()
		{
			return const_cast<reference>(*m_pos);
		}

		template <bool _is_const=is_const>
		typename std::enable_if<_is_const, reference>::type
		operator*() const
		{
			return *m_pos;
		}

	private:
		void advance() {
			while (m_pos != m_end && !m_range_filter.m_filter(*m_pos)) {
				m_pos++;
			}
		}
	};

	using const_iterator = _iterator<true>;
	using iterator = _iterator<false>;

	/**
	 * Creates a subrange of @p range filtering items with the @p filter
	 * predicate
	 * @param range The original range to iterate over
	 * @param filter The unary predicate used to filter elements out of the range
	 */
	range_filter(const Range &range, const UnaryPredicate &filter=UnaryPredicate())
	 : m_range(range), m_filter(filter)
	{ }

	const_iterator begin() const
	{
		return {*this, false};
	}

	const_iterator end() const
	{
		return {*this, true};
	}

	iterator begin()
	{
		return {*this, false};
	}

	iterator end()
	{
		return {*this, true};
	}

	std::size_t size() const
	{
		auto distance = std::distance(begin(), end());
		assert(distance >= 0);
		return std::size_t(distance);
	}

private:
	const Range &m_range;
	const UnaryPredicate m_filter;
};

/**
 * Easy construction of range_filter with type deduction
 *
 * @param range The range to filter
 * @param filter The predicate used for filtering
 * @return The range_filter object providing iteration over filtered values
 */
template <typename Range, typename UnaryPredicate>
range_filter<Range, UnaryPredicate> make_range_filter(const Range &range, const UnaryPredicate &filter)
{
	return range_filter<Range, UnaryPredicate>(range, filter);
}


/**
 * A range-like object
 *
 * This class takes two iterators, which are not necessarily the begin() and
 * end() of a particular range, and makes them look as such. This makes this new
 * range object easy to iterate over range-for loops.
 */
template <typename Iterator>
class range : public std::tuple<Iterator, Iterator>
{
public:
	using iterator = Iterator;
	range(iterator begin, iterator end)
	 : std::tuple<iterator, iterator>(begin, end)
	{ }

	iterator begin() const
	{
		return std::get<0>(*this);
	}

	iterator end() const
	{
		return std::get<1>(*this);
	}

	std::size_t size() const
	{
		auto distance = std::distance(begin(), end());
		assert(distance >= 0);
		return std::size_t(distance);
	}
};

}  // namespace shark



#endif /* SHARK_RANGES_H_ */