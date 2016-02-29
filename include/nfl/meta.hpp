#ifndef NFL_META_HPP
#define NFL_META_HPP

#include <cstddef>

namespace nfl {

// This already exists in boost, but we avoid the dependency doing so!

namespace impl {

template <size_t N>
struct _log2
{
	static constexpr size_t value = 1 + _log2<N/2>::value;
};

template <>
struct _log2<1>
{
	static constexpr size_t value = 0;
};

} // internals

template <size_t N>
struct static_log2
{
	static constexpr size_t value = impl::_log2<N>::value;
};

template <>
struct static_log2<0>
{ };

template <size_t...> struct seq {};
template <size_t N, size_t... S> struct gens : gens<N - 1, N - 1, S...> {};

template <size_t... S> struct gens<0, S...> {
  typedef seq<S...> type;
};


}

#endif
