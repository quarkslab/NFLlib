#ifndef NFL_PERMUT_HPP
#define NFL_PERMUT_HPP

#define PERMUT_LIMIT_UNROLL 1024

#include <limits>

namespace nfl {

namespace details {

// meta-compilation of the bit-reversing permutation
template<size_t H, size_t degree, size_t R, size_t I>
struct r_loop {
  static constexpr size_t value = r_loop<(H << 1), degree, (R << 1) | (I & 1), (I >> 1)>::value;
};
template<size_t degree, size_t R, size_t I>
struct r_loop<degree, degree, R, I> {
  static constexpr size_t value = R;
};
template <size_t I, size_t J, size_t degree>
struct r_set {
  template<class value_type>
  void operator()(value_type *y, value_type const *x) {
    r_set<2*I,2*J, degree>{}(y, x);
    r_set<2*I + 1,2*J, degree>{}(y, x);
  }
};
template <size_t I, size_t degree>
struct r_set<I, degree, degree> {
  template<class value_type>
  void operator()(value_type *y, value_type const *x) {
    constexpr size_t r = r_loop<1, degree, 0u, I>::value;
    y[r] = x[I];
  }
};

template <uint64_t N>
struct uint_value_t
{
  using type = typename std::conditional<
    N <= std::numeric_limits<uint8_t>::max(),
      uint8_t,
      typename std::conditional<
        N <= std::numeric_limits<uint16_t>::max(),
        uint16_t,
        typename std::conditional<
          N <= std::numeric_limits<uint32_t>::max(),
          uint32_t,
          uint64_t>::type
      >::type
   >::type;
};

template <size_t degree>
struct permut_compute
{
  using idx_type = typename uint_value_t<degree>::type;
  idx_type data_[degree];

  permut_compute()
  {
    for (idx_type i = 0; i < degree; ++i) {
      idx_type ii = i;
      idx_type r = 0; 

      for (idx_type h = 1; h < degree; h=h<<1)
      {    
        r = (r << 1) | (ii & 1);
        ii >>= 1;
      }    

      data_[i] = r; 
    }
  }
  
  inline idx_type operator()(size_t i) const {
    assert(i < degree);
    return data_[i];
  }
};

template <size_t degree, bool unroll>
struct permut;

template <size_t degree>
struct permut<degree, true>
{
  template <class V>
  static inline void compute(V* y, V const* x)
  {
    r_set<0, 1, degree>{}(y, x);
  }
};

template <size_t degree>
struct permut<degree, false>
{
  static permut_compute<degree> P;

  template <class V>
  static inline void compute(V* y, V const* x)
  {
    for (size_t i = 0; i < degree; ++i) {
      y[i] = x[P(i)];
    }
  }
};

template <size_t degree>
permut_compute<degree> permut<degree,false>::P;

} // details

template <size_t degree>
struct permut: public details::permut<degree, degree <= PERMUT_LIMIT_UNROLL>
{ };


} // nfl

#endif
