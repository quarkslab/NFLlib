#ifndef NFL_ARCH_SSE_HPP
#define NFL_ARCH_SSE_HPP

#include "nfl/arch/common.hpp"
#include "nfl/algos.hpp"
#include <pmmintrin.h>
#include "nfl/arch/common.hpp"
#include <iostream>

namespace nfl {

namespace simd {

struct sse
{
  template <class T>
  static inline __m128i load(T const* p) { return _mm_load_si128((__m128i const*) p); }

  template <class T>
  static inline void store(T* p, __m128i const v) { _mm_store_si128((__m128i*) p, v); }

  template <class T>
  struct elt_count
  {
    static constexpr size_t value = 16/sizeof(T);
  };

  static constexpr int mode = 1;
};

} // simd

template<> struct common_mode<simd::sse, simd::serial> { using type = simd::serial;};
template<> struct common_mode<simd::serial, simd::sse> { using type = simd::serial;};


namespace ops {

//
// TOOLS
//

static inline __m128i mulhi_epu32(__m128i sse_a, __m128i sse_b)
{
  const __m128i mullow = _mm_srli_epi64(_mm_mul_epu32(sse_a, sse_b), 32);
  sse_a = _mm_shuffle_epi32(sse_a, 1 | (0 << 2) | (3 << 4) | (2 << 6));
  sse_b = _mm_shuffle_epi32(sse_b, 1 | (0 << 2) | (3 << 4) | (2 << 6));
  const __m128i mulhigh = _mm_mul_epu32(sse_a, sse_b);
  return reinterpret_cast<__m128i>(_mm_blend_ps(reinterpret_cast<__m128>(mullow), reinterpret_cast<__m128>(mulhigh), 0b1010));
}

static inline __m128i mulhi_epu16(__m128i sse_a, __m128i sse_b)
{
  return _mm_mulhi_epu16(sse_a, sse_b);
}

template <class Integer>
static inline void assert_strict_mod_sse(__m128i const sse_v, Integer const p)
{
  assert_vec_values<simd::sse, Integer>(sse_v, [p](Integer const v) { return v < p; });
}

template <class Integer>
static inline void assert_addmod_input_sse(__m128i const x, __m128i const y, Integer const p)
{
  assert_addmod_input<simd::sse>(x, y, p);
}

//
// BASIC OPS
//
template<class T>
struct addmod<T, simd::sse> : addmod<T, simd::serial> {};

template<>
struct addmod<uint32_t, simd::sse>
{
  using simd_mode = simd::sse;
  __m128i operator()(__m128i const x, __m128i const y, size_t const cm) const {
    auto const p = params<uint32_t>::P[cm];
    assert_addmod_input_sse<uint32_t>(x, y, p);

    __m128i sse_p = _mm_set1_epi32(p);
    __m128i sse_pc = _mm_set1_epi32(p - 0x80000000 - 1);
    __m128i sse_80 = _mm_set1_epi32(0x80000000);
    const __m128i z = _mm_add_epi32(x, y);
    const __m128i sse_cmp = _mm_cmpgt_epi32(_mm_sub_epi32(z, sse_80), sse_pc);
    const __m128i sse_res = _mm_sub_epi32(z, _mm_and_si128(sse_cmp, sse_p));

    assert_strict_mod_sse<uint32_t>(sse_res, p);
    return sse_res;
  }
};

template<>
struct addmod<uint16_t, simd::sse>
{
  using simd_mode = simd::sse;
  __m128i operator()(__m128i x, __m128i y, size_t const cm) const {
    auto const p = params<uint16_t>::P[cm];
    assert_addmod_input_sse<uint16_t>(x, y, p);

    __m128i sse_p = _mm_set1_epi16(p);
    __m128i sse_pc = _mm_set1_epi16(p - 0x8000 - 1);
    __m128i sse_80 = _mm_set1_epi16(0x8000);
    const __m128i z = _mm_add_epi16(x, y);
    const __m128i sse_cmp = _mm_cmpgt_epi16(_mm_sub_epi16(z, sse_80), sse_pc);
    const __m128i sse_res = _mm_sub_epi16(z, _mm_and_si128(sse_cmp, sse_p));

    assert_strict_mod_sse<uint16_t>(sse_res, p);
    return sse_res;
  }
};

template<class T>
struct mulmod<T, simd::sse> : mulmod<T, simd::serial> {};

template<class T>
struct submod<T, simd::sse> : submod<T, simd::serial> {};

template<>
struct submod<uint16_t, simd::sse>
{
  using simd_mode = simd::sse;
  __m128i operator()(__m128i const x, __m128i const y, size_t const cm) const {
    auto const p = params<uint16_t>::P[cm];
    assert_strict_mod_sse<uint16_t>(x, p);
    assert_strict_mod_sse<uint16_t>(y, p);

    auto const sse_p = _mm_set1_epi16(p);
    auto const sse_res = addmod<uint16_t, simd::sse>{}(x, _mm_sub_epi16(sse_p, y), cm);

    assert_strict_mod_sse<uint16_t>(sse_res, p);
    return sse_res;
  }
};

template<>
struct submod<uint32_t, simd::sse>
{
  using simd_mode = simd::sse;
  __m128i operator()(__m128i const x, __m128i const y, size_t const cm) const {
    auto const p = params<uint32_t>::P[cm];
    assert_strict_mod_sse<uint32_t>(x, p);
    assert_strict_mod_sse<uint32_t>(y, p);

    auto const sse_p = _mm_set1_epi32(p);
    auto const sse_res = addmod<uint32_t, simd::sse>{}(x, _mm_sub_epi32(sse_p, y), cm);

    assert_strict_mod_sse<uint32_t>(sse_res, p);
    return sse_res;
  }
};

template<class T>
struct muladd<T, simd::sse> : muladd<T, simd::serial> {};

//
// NTT
//

template<class poly>
struct ntt_loop_body<simd::sse, poly, uint32_t>
{
  using value_type = typename poly::value_type;
  using greater_value_type = typename poly::greater_value_type;
  using simd_mode = simd::sse;

  ntt_loop_body(value_type const p)
  {
    _sse_2p = _mm_set1_epi32(p << 1);
    _sse_p = _mm_set1_epi32(p);
    _sse_80 = _mm_set1_epi32(0x80000000);
    _sse_2pc = _mm_set1_epi32((2*p) - 0x80000000 - 1);
  }

  inline void operator()(value_type* x0, value_type* x1, value_type const* winvtab, value_type const* wtab) const
  {
    __m128i sse_u0, sse_u1, sse_winvtab, sse_wtab, sse_t0, sse_t1, sse_q, sse_t2, sse_cmp;
    sse_u0 = _mm_load_si128((__m128i const*) x0);
    sse_u1 = _mm_load_si128((__m128i const*) x1);
    sse_winvtab = _mm_load_si128((__m128i const*) winvtab);
    sse_wtab = _mm_load_si128((__m128i const*) wtab);

    sse_t1 = _mm_add_epi32(_sse_2p, _mm_sub_epi32(sse_u0, sse_u1));

    sse_q = mulhi_epu32(sse_t1, sse_winvtab);
    sse_t2 = _mm_sub_epi32(_mm_mullo_epi32(sse_t1, sse_wtab),
        _mm_mullo_epi32(sse_q, _sse_p));

    sse_t0 = _mm_add_epi32(sse_u0, sse_u1);
    sse_cmp = _mm_cmpgt_epi32(_mm_sub_epi32(sse_t0, _sse_80), _sse_2pc);
    sse_t0 = _mm_sub_epi32(sse_t0, _mm_and_si128(sse_cmp, _sse_2p));

    _mm_store_si128((__m128i*) x0, sse_t0);
    _mm_store_si128((__m128i*) x1, sse_t2);
  }

  __m128i _sse_2p;
  __m128i _sse_p;
  __m128i _sse_80;
  __m128i _sse_2pc;
};

template<class poly>
struct ntt_loop_body<simd::sse, poly, uint16_t>
{
  using value_type = typename poly::value_type;
  using greater_value_type = typename poly::greater_value_type;
  using simd_mode = simd::sse;

  ntt_loop_body(value_type const p)
  {
    _sse_2p = _mm_set1_epi16(2*p);
    _sse_p = _mm_set1_epi16(p);
    _sse_80 = _mm_set1_epi16(0x8000);
    _sse_2pc = _mm_set1_epi16((2*p) - 0x8000 - 1);
  }

  inline void operator()(value_type* x0, value_type* x1, value_type const* winvtab, value_type const* wtab) const
  {
    __m128i sse_u0, sse_u1, sse_winvtab, sse_wtab, sse_t0, sse_t1, sse_q, sse_t2, sse_cmp;
    sse_u0 = _mm_load_si128((__m128i const*) x0);
    sse_u1 = _mm_load_si128((__m128i const*) x1);
    sse_winvtab = _mm_load_si128((__m128i const*) winvtab);
    sse_wtab = _mm_load_si128((__m128i const*) wtab);

    sse_t1 = _mm_add_epi16(_sse_2p, _mm_sub_epi16(sse_u0, sse_u1));
    sse_q = mulhi_epu16(sse_t1, sse_winvtab);

    sse_t2 = _mm_sub_epi16(
        _mm_mullo_epi16(sse_t1, sse_wtab),
        _mm_mullo_epi16(sse_q, _sse_p));

    sse_t0 = _mm_add_epi16(sse_u0, sse_u1);
    sse_cmp = _mm_cmpgt_epi16(_mm_sub_epi16(sse_t0, _sse_80), _sse_2pc);
    sse_t0 = _mm_sub_epi16(sse_t0, _mm_and_si128(sse_cmp, _sse_2p));

    _mm_store_si128((__m128i*) x0, sse_t0);
    _mm_store_si128((__m128i*) x1, sse_t2);
  }

  __m128i _sse_2p;
  __m128i _sse_p;
  __m128i _sse_80;
  __m128i _sse_2pc;
};

// Default implementation
template <class poly, class T>
struct ntt_loop<simd::sse, poly, T>: ntt_loop<simd::serial, poly, T> {};

template<class poly>
struct ntt_loop_sse_unrolled {
  using value_type = typename poly::value_type;
  using greater_value_type = typename poly::greater_value_type;
  using simd_mode = simd::sse;

  static constexpr size_t J = static_log2<poly::degree>::value-2;

#ifdef __GNUC__
  __attribute__((optimize("no-tree-vectorize")))
#endif
  static size_t run(value_type* x, const value_type* &wtab, const value_type* &winvtab, const value_type p) {
    ntt_loop_body<simd::sse, poly, value_type> body_sse(p);
    ntt_loop_body<simd::serial, poly, value_type> body(p);

    for (size_t w = 0; w < J-1; w++) {
      const size_t M = 1 << w;
      const size_t N = poly::degree >> w;
      for (size_t r = 0; r < M; r++) {
        for (size_t i = 0; i < N/2; i += simd_mode::elt_count<value_type>::value) {
          body_sse(&x[N * r + i], &x[N * r + i + N/2], &winvtab[i], &wtab[i]);
        }
      }
      wtab += N / 2;
      winvtab += N / 2;
    }

    // Special case for w=J-1: not able to fill a register of data!
    constexpr size_t w = J-1;
    constexpr size_t M = 1 << w;
    constexpr size_t N = poly::degree >> w;
    for (size_t r = 0; r < M; r++) {
      for (size_t i = 0; i < N/2; i++) {
        body(&x[N * r + i], &x[N * r + i + N/2], &winvtab[i], &wtab[i]);
      }
    }
    wtab += N / 2;
    winvtab += N / 2;

    return 1<<J;
  }
};
template<class poly>
struct ntt_loop<simd::sse, poly, uint32_t>: public ntt_loop_sse_unrolled<poly>
{ };

template<class poly>
struct ntt_loop<simd::sse, poly, uint16_t>: public ntt_loop_sse_unrolled<poly>
{ };

//
// MULMOD_SHOUP
//
template<class T>
struct mulmod_shoup<T, simd::sse> : mulmod_shoup<T, simd::serial> {};

template<>
struct mulmod_shoup<uint32_t, simd::sse>
{
  using value_type = uint32_t;
  using simd_mode = simd::sse;

  constexpr static inline size_t elt_count() { return simd_mode::elt_count<value_type>::value; }

  static inline __m128i finish(__m128i sse_x, __m128i sse_y, __m128i sse_q, __m128i sse_p_32, __m128i sse_p_64, __m128i sse_pc_64, __m128i sse_80)
  {
      //greater_value_type res;
      //res = x * y - q * p;
      //res -= ((res>=p) ? p : 0);
    const __m128i sse_res = _mm_sub_epi64(
        _mm_mul_epu32(sse_x, sse_y),
        _mm_mul_epu32(sse_q, sse_p_32)
      );

    const __m128i sse_cmp = _mm_cmpgt_epi64(_mm_sub_epi64(sse_res, sse_80), sse_pc_64);
    return _mm_sub_epi64(sse_res, _mm_and_si128(sse_cmp, sse_p_64));
  }


  static inline __m128i shuffle_lh(__m128i const v)
  {
    return _mm_shuffle_epi32(v, 1 | (0 << 2) | (3 << 4) | (2 << 6));
  }

  inline __m128i operator()(__m128i const sse_x, __m128i const sse_y, __m128i const sse_yprime, size_t const cm) const
  {
    auto const p = params<value_type>::P[cm];
    assert_strict_mod_sse<uint32_t>(sse_x, p);
    assert_strict_mod_sse<uint32_t>(sse_y, p);

    __m128i sse_p_32;
    __m128i sse_p_64;
    __m128i sse_pc_64;
    __m128i sse_80;
    __m128i sse_q, res1, res2;
    sse_p_32 = _mm_set1_epi32(p);
    sse_p_64 = _mm_set1_epi64x(p);
    sse_pc_64 = _mm_set1_epi64x((uint64_t)(p)-0x8000000000000000ULL-1);
    sse_80 = _mm_set1_epi64x(0x8000000000000000ULL);
    sse_q = mulhi_epu32(sse_x, sse_yprime);

    res1 = finish(sse_x, sse_y, sse_q, sse_p_32, sse_p_64, sse_pc_64, sse_80);
    res2 = finish(shuffle_lh(sse_x), shuffle_lh(sse_y), shuffle_lh(sse_q), sse_p_32, sse_p_64, sse_pc_64, sse_80);
	res2 = _mm_slli_epi64(res2, 32);
    const __m128i sse_res = reinterpret_cast<__m128i>(_mm_blend_ps(reinterpret_cast<__m128>(res1), reinterpret_cast<__m128>(res2), 0b1010));

    assert_strict_mod_sse<uint32_t>(sse_res, p);
    return sse_res;
  }

};

template<>
struct mulmod_shoup<uint16_t, simd::sse>
{
  using value_type = uint16_t;
  using simd_mode = simd::sse;

  constexpr static inline size_t elt_count() { return simd_mode::elt_count<value_type>::value; }

  static inline __m128i finish(__m128i sse_x, __m128i sse_y, __m128i sse_q, __m128i sse_p_32, __m128i sse_pc_32, __m128i sse_80)
  {
      //greater_value_type res;
      //res = x * y - q * p;
      //res -= ((res>=p) ? p : 0);

      const __m128i sse_x_32l = _mm_cvtepu16_epi32(sse_x);
      const __m128i sse_y_32l = _mm_cvtepu16_epi32(sse_y);
      const __m128i sse_q_32l = _mm_cvtepu16_epi32(sse_q);

      const __m128i sse_res = _mm_sub_epi32(
        _mm_mullo_epi32(sse_x_32l, sse_y_32l),
        _mm_mullo_epi32(sse_q_32l, sse_p_32)
      );

      const __m128i sse_cmp = _mm_cmpgt_epi32(_mm_sub_epi32(sse_res, sse_80), sse_pc_32);
      return _mm_sub_epi32(sse_res, _mm_and_si128(sse_cmp, sse_p_32));
  }

  static inline __m128i shift8(__m128i const v)
  {
    return _mm_srli_si128(v, 8);
  }

  inline __m128i operator()(__m128i const sse_x, __m128i const sse_y, __m128i const sse_yprime, size_t const cm) const
  {
    auto const p = params<value_type>::P[cm];
    assert_strict_mod_sse<uint16_t>(sse_x, p);
    assert_strict_mod_sse<uint16_t>(sse_y, p);

    __m128i sse_p_32;
    __m128i sse_pc_32;
    __m128i sse_80;
    __m128i sse_q, res1, res2;
    sse_p_32 = _mm_set1_epi32(p);
    sse_pc_32 = _mm_set1_epi32((uint32_t)(p)-0x80000000-1);
    sse_80 = _mm_set1_epi32(0x80000000);
    sse_q = mulhi_epu16(sse_x, sse_yprime);

    // Cast to 32-bits and use mullo_epi32. Do it twice and merge!
    res1 = finish(sse_x, sse_y, sse_q, sse_p_32, sse_pc_32, sse_80);
    res2 = finish(shift8(sse_x), shift8(sse_y), shift8(sse_q), sse_p_32, sse_pc_32, sse_80);

    const __m128i sse_res = _mm_packus_epi32(res1, res2);
    assert_strict_mod_sse<uint16_t>(sse_res, p);
    return sse_res;
  }

};

//
// MULADD_SHOUP
//

template<class T>
struct muladd_shoup<T, simd::sse> : muladd_shoup<T, simd::serial> {};

template<>
struct muladd_shoup<uint16_t, simd::sse>
{
  using value_type = uint16_t;
  using simd_mode = simd::sse;

  constexpr static inline size_t elt_count() { return simd_mode::elt_count<value_type>::value; }

  static inline __m128i finish(__m128i sse_rop, __m128i sse_x, __m128i sse_y, __m128i sse_q, __m128i sse_p_32, __m128i sse_pc_32, __m128i sse_80)
  {
      //greater_value_type res;
      //res = x * y - q * p;
      //res -= ((res>=p) ? p : 0);

      const __m128i sse_x_32l = _mm_cvtepu16_epi32(sse_x);
      const __m128i sse_y_32l = _mm_cvtepu16_epi32(sse_y);
      const __m128i sse_q_32l = _mm_cvtepu16_epi32(sse_q);
      const __m128i sse_rop_32l = _mm_cvtepu16_epi32(sse_rop);

      const __m128i sse_res = _mm_add_epi32(sse_rop_32l,
        _mm_sub_epi32(
          _mm_mullo_epi32(sse_x_32l, sse_y_32l),
          _mm_mullo_epi32(sse_q_32l, sse_p_32)
        )
      );

      const __m128i sse_cmp = _mm_cmpgt_epi32(_mm_add_epi32(sse_res, sse_80), sse_pc_32);
      return _mm_sub_epi32(sse_res, _mm_and_si128(sse_cmp, sse_p_32));
  }

  static inline __m128i shift8(__m128i const v)
  {
    return _mm_srli_si128(v, 8);
  }

  inline __m128i operator()(__m128i const sse_rop, __m128i const sse_x, __m128i const sse_y, __m128i const sse_yprime, size_t const cm) const
  {
    auto const p = params<uint16_t>::P[cm];
    assert_strict_mod_sse<uint16_t>(sse_x, p);
    assert_strict_mod_sse<uint16_t>(sse_y, p);
    assert_strict_mod_sse<uint16_t>(sse_rop, p);

    __m128i _sse_p_32;
    __m128i _sse_pc_32;
    __m128i _sse_80;
    _sse_p_32 = _mm_set1_epi32(p);
    _sse_pc_32 = _mm_set1_epi32((uint32_t)(p)-0x80000000-1);
    _sse_80 = _mm_set1_epi32(0x80000000);
    __m128i sse_q, res1, res2;
    sse_q = mulhi_epu16(sse_x, sse_yprime);

    // Cast to 32-bits and use mullo_epi32. Do it twice and merge!
    res1 = finish(sse_rop, sse_x, sse_y, sse_q, _sse_p_32, _sse_pc_32, _sse_80);
    res2 = finish(shift8(sse_rop), shift8(sse_x), shift8(sse_y), shift8(sse_q), _sse_p_32, _sse_pc_32, _sse_80);

    const __m128i sse_res = _mm_packus_epi32(res1, res2);
    assert_strict_mod_sse<uint16_t>(sse_res, p);
    return sse_res;
  }

};

template<class T>
struct shoup<T, simd::sse> {
  using simd_mode = simd::sse;
};

template<class T>
struct compute_shoup<T, simd::sse> : compute_shoup<T, simd::serial> {};

} // ops

} // nfl

#endif
