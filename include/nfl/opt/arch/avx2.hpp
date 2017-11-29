#ifndef DEF_NFLImplems_avx2
#define DEF_NFLImplems_avx2

#include <immintrin.h>
#include "nfl/opt/arch/sse.hpp"

namespace nfl {

namespace simd {

struct avx2 {
  template <class T>
  static inline __m256i load(T const* p) { return _mm256_load_si256((__m256i const*) p); }

  template <class T>
  static inline void store(T* p, __m256i const v) { _mm256_store_si256((__m256i*) p, v); }

  template <class T>
  struct elt_count
  {
    static constexpr size_t value = 32/sizeof(T);
  };

  static constexpr int mode = 2;
};

} // simd

template<> struct common_mode<simd::avx2, simd::serial> { using type = simd::serial;};
template<> struct common_mode<simd::avx2, simd::sse> { using type = simd::sse;};
template<> struct common_mode<simd::sse, simd::avx2> { using type = simd::sse;};
template<> struct common_mode<simd::serial, simd::avx2> { using type = simd::serial;};

namespace ops {

//
// TOOLS
//

static inline __m256i avx2_mulhi_epu32(__m256i avx_a, __m256i avx_b)
{
  const __m256i mullow = _mm256_srli_epi64(_mm256_mul_epu32(avx_a, avx_b), 32);
  avx_a = _mm256_shuffle_epi32(avx_a, 1 | (0 << 2) | (3 << 4) | (2 << 6));
  avx_b = _mm256_shuffle_epi32(avx_b, 1 | (0 << 2) | (3 << 4) | (2 << 6));
  const __m256i mulhigh = _mm256_mul_epu32(avx_a, avx_b);
  return reinterpret_cast<__m256i>(_mm256_blend_ps(reinterpret_cast<__m256>(mullow), reinterpret_cast<__m256>(mulhigh), 0b10101010));
}

template <class Integer>
static inline void assert_strict_mod_avx2(__m256i const avx_v, Integer const p)
{
  assert_vec_values<simd::avx2, Integer>(avx_v, [p](Integer const v) { return v < p; });
}

template <class Integer>
static inline void assert_addmod_input_avx2(__m256i const x, __m256i const y, Integer const p)
{
  assert_addmod_input<simd::avx2>(x, y, p);
}


//
// BASIC OPS
//

template<class T>
struct addmod<T, simd::avx2> : addmod<T, simd::sse> {};

template<>
struct addmod<uint32_t, simd::avx2>
{
  using simd_mode = simd::avx2;
  __m256i operator()(__m256i x, __m256i y, size_t const cm) const {
    auto const p = params<uint32_t>::P[cm];
    assert_addmod_input_avx2<uint32_t>(x, y, p);

    __m256i avx_p = _mm256_set1_epi32(p);
    __m256i avx_pc = _mm256_set1_epi32(p - 0x80000000 - 1);
    __m256i avx_80 = _mm256_set1_epi32(0x80000000);
    const __m256i z = _mm256_add_epi32(x, y);
    const __m256i avx_cmp = _mm256_cmpgt_epi32(_mm256_sub_epi32(z, avx_80), avx_pc);
    const __m256i avx_res = _mm256_sub_epi32(z, _mm256_and_si256(avx_cmp, avx_p));

    assert_strict_mod_avx2<uint32_t>(avx_res, p);
    return avx_res;
  }
};

template<>
struct addmod<uint16_t, simd::avx2>
{
  using simd_mode = simd::avx2;
  __m256i operator()(__m256i x, __m256i y, size_t const cm) const {
    auto const p = params<uint16_t>::P[cm];
    assert_addmod_input_avx2<uint16_t>(x, y, p);

    __m256i avx_p = _mm256_set1_epi16(p);
    __m256i avx_pc = _mm256_set1_epi16(p - 0x8000 - 1);
    __m256i avx_80 = _mm256_set1_epi16(0x8000);
    const __m256i z = _mm256_add_epi16(x, y);
    const __m256i avx_cmp = _mm256_cmpgt_epi16(_mm256_sub_epi16(z, avx_80), avx_pc);
    const __m256i avx_res = _mm256_sub_epi16(z, _mm256_and_si256(avx_cmp, avx_p));

    assert_strict_mod_avx2<uint16_t>(avx_res, p);
    return avx_res;
  }
};

template<class T>
struct mulmod<T, simd::avx2> : mulmod<T, simd::sse> {};

template<class T>
struct submod<T, simd::avx2> : submod<T, simd::sse> {};

template<>
struct submod<uint16_t, simd::avx2>
{
  using simd_mode = simd::avx2;
  __m256i operator()(__m256i const x, __m256i const y, size_t const cm) const {
    auto const p = params<uint16_t>::P[cm];
    assert_strict_mod_avx2<uint16_t>(x, p);
    assert_strict_mod_avx2<uint16_t>(y, p);

    auto const avx_p = _mm256_set1_epi16(p);
    auto const avx_res = addmod<uint16_t, simd::avx2>{}(x, _mm256_sub_epi16(avx_p, y), cm);

    assert_strict_mod_avx2<uint16_t>(avx_res, p);
    return avx_res;
  }
};

template<>
struct submod<uint32_t, simd::avx2>
{
  using simd_mode = simd::avx2;
  __m256i operator()(__m256i const x, __m256i const y, size_t const cm) const {
    auto const p = params<uint32_t>::P[cm];
    assert_strict_mod_avx2<uint32_t>(x, p);
    assert_strict_mod_avx2<uint32_t>(y, p);

    auto const avx_p = _mm256_set1_epi32(p);
    auto const avx_res = addmod<uint32_t, simd::avx2>{}(x, _mm256_sub_epi32(avx_p, y), cm);

    assert_strict_mod_avx2<uint32_t>(avx_res, p);
    return avx_res;
  }
};

template<class T>
struct muladd<T, simd::avx2> : muladd<T, simd::sse> {};

//
// NTT
//

template <class poly, class T>
struct ntt_loop<simd::avx2, poly, T>: public ntt_loop<simd::sse, poly, T>
{ };

template<class poly>
struct ntt_loop_body<simd::avx2, poly, uint32_t> {
  using value_type = typename poly::value_type;
  using greater_value_type = typename poly::greater_value_type;
  using simd_type = simd::avx2;

  ntt_loop_body(value_type const p)
  {
    _avx_2p = _mm256_set1_epi32(p << 1);
    _avx_p = _mm256_set1_epi32(p);
    _avx_80 = _mm256_set1_epi32(0x80000000);
    _avx_2pc = _mm256_set1_epi32((2*p) - 0x80000000 - 1);
  }

  inline void operator()(value_type* x0, value_type* x1, value_type const* winvtab, value_type const* wtab) const
  {
    __m256i avx_u0, avx_u1, avx_winvtab, avx_wtab, avx_t0, avx_t1, avx_q, avx_t2, avx_cmp;
    avx_u0 = _mm256_load_si256((__m256i const*) x0);
    avx_u1 = _mm256_load_si256((__m256i const*) x1);
    avx_winvtab = _mm256_load_si256((__m256i const*) winvtab);
    avx_wtab = _mm256_load_si256((__m256i const*) wtab);

    avx_t1 = _mm256_add_epi32(_avx_2p, _mm256_sub_epi32(avx_u0, avx_u1));

    avx_q = avx2_mulhi_epu32(avx_t1, avx_winvtab);
    avx_t2 = _mm256_sub_epi32(_mm256_mullo_epi32(avx_t1, avx_wtab),
        _mm256_mullo_epi32(avx_q, _avx_p));

    avx_t0 = _mm256_add_epi32(avx_u0, avx_u1);
    avx_cmp = _mm256_cmpgt_epi32(_mm256_sub_epi32(avx_t0, _avx_80), _avx_2pc);
    avx_t0 = _mm256_sub_epi32(avx_t0, _mm256_and_si256(avx_cmp, _avx_2p));

    _mm256_store_si256((__m256i*) x0, avx_t0);
    _mm256_store_si256((__m256i*) x1, avx_t2);
  }

  __m256i _avx_2p;
  __m256i _avx_p;
  __m256i _avx_80;
  __m256i _avx_2pc;
};

template<class poly>
struct ntt_loop_body<simd::avx2, poly, uint16_t> {
  using value_type = typename poly::value_type;
  using greater_value_type = typename poly::greater_value_type;
  using simd_type = simd::avx2;

  ntt_loop_body(value_type const p)
  {
    _avx_2p = _mm256_set1_epi16(p << 1);
    _avx_p = _mm256_set1_epi16(p);
    _avx_80 = _mm256_set1_epi16(0x8000);
    _avx_2pc = _mm256_set1_epi16((2*p) - 0x8000 - 1);
  }

  inline void operator()(value_type* x0, value_type* x1, value_type const* winvtab, value_type const* wtab) const
  {
    __m256i avx_u0, avx_u1, avx_winvtab, avx_wtab, avx_t0, avx_t1, avx_q, avx_t2, avx_cmp;
    avx_u0 = _mm256_load_si256((__m256i const*) x0);
    avx_u1 = _mm256_load_si256((__m256i const*) x1);
    avx_winvtab = _mm256_load_si256((__m256i const*) winvtab);
    avx_wtab = _mm256_load_si256((__m256i const*) wtab);

    avx_t1 = _mm256_add_epi16(_avx_2p, _mm256_sub_epi16(avx_u0, avx_u1));

    avx_q = _mm256_mulhi_epu16(avx_t1, avx_winvtab);
    avx_t2 = _mm256_sub_epi16(_mm256_mullo_epi16(avx_t1, avx_wtab),
        _mm256_mullo_epi16(avx_q, _avx_p));

    avx_t0 = _mm256_add_epi16(avx_u0, avx_u1);
    avx_cmp = _mm256_cmpgt_epi16(_mm256_sub_epi16(avx_t0, _avx_80), _avx_2pc);
    avx_t0 = _mm256_sub_epi16(avx_t0, _mm256_and_si256(avx_cmp, _avx_2p));

    _mm256_store_si256((__m256i*) x0, avx_t0);
    _mm256_store_si256((__m256i*) x1, avx_t2);
  }

  __m256i _avx_2p;
  __m256i _avx_p;
  __m256i _avx_80;
  __m256i _avx_2pc;
};

template<class poly>
struct ntt_loop_avx2_unrolled {
  using value_type = typename poly::value_type;
  using greater_value_type = typename poly::greater_value_type;
  using simd_type = simd::avx2;

  static constexpr size_t J = static_log2<poly::degree>::value-2;
  static constexpr size_t elt_count = simd_type::elt_count<value_type>::value;
  static constexpr size_t elt_count_sse = simd::sse::elt_count<value_type>::value;

#ifdef __GNUC__
  __attribute__((optimize("no-tree-vectorize")))
#endif
  static size_t run(value_type* x, const value_type* &wtab,
      const value_type* &winvtab, const value_type p) {
    ntt_loop_body<simd::avx2, poly, value_type> body_avx2(p);
    ntt_loop_body<simd::sse, poly, value_type> body_sse(p);
    ntt_loop_body<simd::serial, poly, value_type> body(p);

    for (size_t w = 0; w < J-1; w++) {
      const size_t M = 1 << w;
      const size_t N = poly::degree >> w;
      for (size_t r = 0; r < M; r++) {
        const size_t Navx2 = ((N/2)/elt_count)*elt_count;
        for (size_t i = 0; i < Navx2; i += elt_count) {
          body_avx2(&x[N * r + i], &x[N * r + i + N/2], &winvtab[i], &wtab[i]);
        }
        if (Navx2 != (N/2)) {
          const size_t i = Navx2;
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
struct ntt_loop<simd::avx2, poly, uint32_t>: public ntt_loop_avx2_unrolled<poly>
{ };

template<class poly>
struct ntt_loop<simd::avx2, poly, uint16_t>: public ntt_loop_avx2_unrolled<poly>
{ };

//
// MULMOD_SHOUP
//

template<class T>
struct mulmod_shoup<T, simd::avx2> : mulmod_shoup<T, simd::sse> {};

template<>
struct mulmod_shoup<uint16_t, simd::avx2>
{
  using value_type = uint16_t;
  // AG: this is wanted, as loads/stores are still done using SSE registers
  using simd_mode = simd::sse;

  constexpr static inline size_t elt_count() { return simd_mode::elt_count<value_type>::value; }

  static inline __m128i finish(__m128i sse_x, __m128i sse_y, __m128i sse_q, __m256i avx_p_32, __m256i avx_pc_32, __m256i avx_80)
  {
    const __m256i avx_x_32 = _mm256_cvtepu16_epi32(sse_x);
    const __m256i avx_y_32 = _mm256_cvtepu16_epi32(sse_y);
    const __m256i avx_q_32 = _mm256_cvtepu16_epi32(sse_q);

    const __m256i avx_res = _mm256_sub_epi32(
        _mm256_mullo_epi32(avx_x_32, avx_y_32),
        _mm256_mullo_epi32(avx_q_32, avx_p_32)
        );

    const __m256i avx_cmp = _mm256_cmpgt_epi32(_mm256_sub_epi32(avx_res, avx_80), avx_pc_32);
    const __m256i tmp1 = _mm256_sub_epi32(avx_res, _mm256_and_si256(avx_cmp, avx_p_32));
    const __m256i tmp2 = _mm256_permute2x128_si256(tmp1, tmp1, 1);
    const __m128i ret = _mm_packus_epi32(_mm256_castsi256_si128(tmp1), _mm256_castsi256_si128(tmp2));
    return ret;
  }

#if __GNUC__ == 4 && __GNUC_MINOR__ == 8
  // AG: there is a weird bug with GCC 4.8 that seems to generate incorrect code here.
  // adding a noinline attribute fixes the bug, but we need to get to the bottom of this...
  // This works fine with GCC 4.9 and clang 3.5! 
  __attribute__((noinline))
#endif
  __m128i operator()(__m128i const sse_x, __m128i const sse_y, __m128i const sse_yprime, size_t const cm) const
  {
    auto const p = params<value_type>::P[cm];
    assert_strict_mod_sse<value_type>(sse_x, p);
    assert_strict_mod_sse<value_type>(sse_y, p);

    __m256i avx_p_32;
    __m256i avx_pc_32;
    __m256i avx_80;
    avx_p_32 = _mm256_set1_epi32(p);
    avx_pc_32 = _mm256_set1_epi32((uint32_t)(p)-0x80000000-1);
    avx_80 = _mm256_set1_epi32(0x80000000);

    __m128i sse_q = mulhi_epu16(sse_x, sse_yprime);
    __m128i sse_res = finish(sse_x, sse_y, sse_q, avx_p_32, avx_pc_32, avx_80);

    assert_strict_mod_sse<value_type>(sse_res, p);
    return sse_res;
  }

};

//
// MULADD_SHOUP
//

template<class T>
struct muladd_shoup<T, simd::avx2> : muladd_shoup<T, simd::sse> {};

template<>
struct muladd_shoup<uint16_t, simd::avx2>
{
  using value_type = uint16_t;
  // AG: this is wanted, as loads/stores are still done using SSE registers
  using simd_mode = simd::sse;

  constexpr static inline size_t elt_count() { return simd_mode::elt_count<value_type>::value; }

  static inline __m128i finish(__m128i sse_rop, __m128i sse_x, __m128i sse_y, __m128i sse_q, __m256i avx_p_32, __m256i avx_pc_32, __m256i avx_80)
  {
    const __m256i avx_x_32 = _mm256_cvtepu16_epi32(sse_x);
    const __m256i avx_y_32 = _mm256_cvtepu16_epi32(sse_y);
    const __m256i avx_q_32 = _mm256_cvtepu16_epi32(sse_q);
    const __m256i avx_rop_32 = _mm256_cvtepu16_epi32(sse_rop);

    const __m256i avx_res = _mm256_add_epi32(avx_rop_32,
      _mm256_sub_epi32(
        _mm256_mullo_epi32(avx_x_32, avx_y_32),
        _mm256_mullo_epi32(avx_q_32, avx_p_32)
      )
    );

    const __m256i avx_cmp = _mm256_cmpgt_epi32(_mm256_sub_epi32(avx_res, avx_80), avx_pc_32);
    const __m256i tmp1 = _mm256_sub_epi32(avx_res, _mm256_and_si256(avx_cmp, avx_p_32));
    const __m256i tmp2 = _mm256_permute2x128_si256(tmp1, tmp1, 1);
    const __m128i ret = _mm_packus_epi32(_mm256_castsi256_si128(tmp1), _mm256_castsi256_si128(tmp2));
    return ret;
  }

  inline __m128i operator()(__m128i const rop, __m128i const sse_x, __m128i const sse_y, __m128i const sse_yprime, size_t const cm) const
  {
    auto const p = params<value_type>::P[cm];
    assert_strict_mod_sse<value_type>(sse_x, p);
    assert_strict_mod_sse<value_type>(sse_y, p);
    assert_strict_mod_sse<value_type>(rop, p);

    __m256i avx_p_32;
    __m256i avx_pc_32;
    __m256i avx_80;
    avx_p_32 = _mm256_set1_epi32(p);
    avx_pc_32 = _mm256_set1_epi32((uint32_t)(p)-0x80000000-1);
    avx_80 = _mm256_set1_epi32(0x80000000);
    __m128i sse_q = mulhi_epu16(sse_x, sse_yprime);
    __m128i sse_res = finish(rop, sse_x, sse_y, sse_q, avx_p_32, avx_pc_32, avx_80);

    assert_strict_mod_sse<value_type>(sse_res, p);
    return sse_res;
  }

};

template<class T>
struct shoup<T, simd::avx2> {
  using simd_mode = simd::avx2;
};

template<class T>
struct compute_shoup<T, simd::avx2> : compute_shoup<T, simd::sse> {};

} // ops

} // ntt


#endif

// END

