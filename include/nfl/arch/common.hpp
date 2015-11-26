#ifndef NFL_ARCH_COMMON_HPP
#define NFL_ARCH_COMMON_HPP

#include <immintrin.h>
#include <utility>

namespace nfl {

namespace simd {

struct serial
{
  template <class T>
  static inline T load(T const* p) { return *p; }

  template <class T>
  static inline void store(T* p, T const v) { *p = v; }

  template <class T>
  struct elt_count
  {
    static constexpr size_t value = 1;
  };

  static constexpr int mode = 0;
};


} // simd

template<class... M> struct common_mode;
template<class M> struct common_mode<M> { using type = M; };
template<class M0, class... M> struct common_mode<M0, M...> { using type = typename common_mode<M0, typename common_mode<M...>::type>::type; };
template<class M> struct common_mode<M, M> { using type = M;};




template <class Body, class T, class... Args>
static inline void loop_body_runner(Body& body, T* ret, Args&& ... args)
{
  auto ret_ = body(Body::simd_type::load(std::forward<Args>(args))...);
  Body::simd_type::store(ret, ret_);
}

#ifdef CHECK_STRICTMOD
template <class SIMD, class Integer, class Vec, class F>
void assert_vec_values(Vec const v, F const& f)
{
  alignas(32) Integer v_[SIMD::template elt_count<Integer>::value];
  SIMD::store(&v_[0], v);
  for (Integer const i: v_) {
    assert(f(i));
  }
}

template <class SIMD, class Integer, class Vec>
void assert_addmod_input(Vec const x, Vec const y, Integer const p)
{
  static constexpr size_t n = SIMD::template elt_count<Integer>::value;
  alignas(32) Integer x_[n];
  alignas(32) Integer y_[n];
  SIMD::store(&x_[0], x);
  SIMD::store(&y_[0], y);
  for (size_t i = 0; i < n; i++) {
    assert(((size_t)x_[i] + y_[i]) < (2*p));
  }
}
#else
template <class SIMD, class Integer, class Vec, class F>
void assert_vec_values(Vec const, F const&)
{ }

template <class SIMD, class Integer, class Vec>
void assert_addmod_input(Vec const x, Vec const y, Integer const p)
{ }
#endif

} // nfl

#endif
