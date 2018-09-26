/* Copyright (C) 2015  Carlos Aguilar, Tancrède Lepoint, Adrien Guinet and Serge Guelton
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#ifndef NFL_POLY_HPP
#define NFL_POLY_HPP

/***
 * Polynomial class for NFL
 *
 * The poly class is used to manipulate blah blah
 */

#include "nfl/debug.hpp"
#include "nfl/meta.hpp"
#include "nfl/params.hpp"
#include "nfl/ops.hpp"
#include "nfl/arch.hpp"
#include "nfl/prng/fastrandombytes.h"
#include "nfl/prng/FastGaussianNoise.hpp"

#include <iostream>
#include <algorithm>
#include <gmpxx.h>
#include <array>


namespace nfl {

/***
 * Generators to initialize random polynomials
 */
struct uniform {};

struct non_uniform {
  uint64_t upper_bound;
  uint64_t amplifier;
  non_uniform(uint64_t ub) : upper_bound{ub}, amplifier{1} {}
  non_uniform(uint64_t ub, uint64_t amp) : upper_bound{ub}, amplifier{amp} {}
};

struct hwt_dist { // hamming weight distribution.
  uint32_t hwt;
  hwt_dist(uint32_t hwt_) : hwt(hwt_) {}
};

struct ZO_dist { // zero distribution.
  uint8_t rho; // P(1) = P(-1) = (rho/0xFF)/2, P(0) = 1 - P(1) - P(-1)
  ZO_dist(uint8_t rho_ = 0x7F) : rho(rho_) {}
};

template<class in_class, class out_class, unsigned _lu_depth>
struct gaussian {
  FastGaussianNoise<in_class, out_class, _lu_depth> *fg_prng;
  uint64_t amplifier;
  gaussian(FastGaussianNoise<in_class, out_class, _lu_depth> *prng) : fg_prng{prng}, amplifier{1} {}
  gaussian(FastGaussianNoise<in_class, out_class, _lu_depth> *prng, uint64_t amp) : fg_prng{prng}, amplifier{amp} {}
};

// Forward declaration for proxy class used in tests to access poly
// protected/private function members
namespace tests {

template <class P>
class poly_tests_proxy;

} // tests

/* Core polynomial class an array of value_types
 * No indication of whether these are coefficients or values (in NTT form)
 * The developer must keep track of which representation each object is under
 */
template<class T, size_t Degree, size_t NbModuli>
class poly {

  template <class P> friend class tests::poly_tests_proxy;

  static constexpr size_t N = Degree * NbModuli;
  T _data[N] __attribute__((aligned(32)));

public:
  using value_type = typename params<T>::value_type;
  using greater_value_type = typename params<T>::greater_value_type;
  using signed_value_type = typename params<T>::signed_value_type;
  using pointer_type = T*;
  using const_pointer_type = T const*;

  using iterator = pointer_type;
  using const_iterator = const_pointer_type;

  using simd_mode = CC_SIMD;
  static constexpr size_t degree = Degree;
  static constexpr size_t nmoduli = NbModuli;
  static constexpr size_t nbits = params<T>::kModulusBitsize;
  static constexpr size_t aggregated_modulus_bit_size = NbModuli * nbits;

public:
  /* constructors
   */
  poly();
  poly(uniform const& mode);
  poly(non_uniform const& mode);
  poly(hwt_dist const& mode);
  poly(ZO_dist const& mode);
  poly(value_type v, bool reduce_coeffs = true);
  poly(std::initializer_list<value_type> values, bool reduce_coeffs = true);
  template <class It> poly(It first, It last, bool reduce_coeffs = true);
  template <class Op, class... Args> poly(ops::expr<Op, Args...> const& expr);
  template <class in_class, unsigned _lu_depth> poly(gaussian<in_class, T, _lu_depth> const& mode);

  void set(uniform const& mode);
  void set(non_uniform const& mode);
  void set(hwt_dist const& mode);
  void set(ZO_dist const& mode);
  void set(value_type v, bool reduce_coeffs = true);
  void set(std::initializer_list<value_type> values, bool reduce_coeffs = true);
  template <class It> void set(It first, It last, bool reduce_coeffs = true);
  template <class in_class, unsigned _lu_depth> void set(gaussian<in_class, T, _lu_depth> const& mode);
  
  /* assignment
   */
  poly& operator=(value_type v) { set(v); return *this; }
  poly& operator=(uniform const& mode) { set(mode); return *this; }
  poly& operator=(non_uniform const& mode) { set(mode); return *this; }
  poly& operator=(hwt_dist const& mode) { set(mode); return *this; }
  poly& operator=(ZO_dist const& mode) { set(mode); return *this; }
  poly& operator=(std::initializer_list<value_type> values) { set(values); return *this; }
  template <class in_class, unsigned _lu_depth> poly& operator=(gaussian<in_class, T, _lu_depth> const& mode) { set(mode); return *this; }
  template <class Op, class... Args> poly& operator=(ops::expr<Op, Args...> const& expr);

  /* conversion operators
   */
  explicit operator bool() const;

  /* iterators
   */
  iterator begin() { return std::begin(_data); }
  iterator end() { return std::end(_data); }
  const_iterator begin() const { return std::begin(_data); }
  const_iterator end() const { return std::end(_data); }
  const_iterator cbegin() const { return std::begin(_data); }
  const_iterator cend() const { return std::end(_data); }


  /* polynomial indexing
   */
  value_type const& operator()(size_t cm, size_t i) const { return _data[cm * degree + i]; }
  value_type& operator()(size_t cm, size_t i) { return _data[cm * degree + i]; }
  template<class M> auto load(size_t cm, size_t i) const -> decltype(M::load(&(this->operator()(cm, i)))) { return M::load(&(*this)(cm, i)); }

  /* misc
   */
  pointer_type data() { return _data; }
  static constexpr value_type get_modulus(size_t n) { return params<T>::P[n]; }

  /* ntt stuff - public API
   */
  void ntt_pow_phi() { base.ntt_pow_phi(*this);}
  void invntt_pow_invphi() { base.invntt_pow_invphi(*this); }

  // Serialization API
  //
  // Serialization of polynomials is not portable in terms of "cross
  // architecture portability".
  // Therefore, serialization should not be used in the rare cases a big-endian
  // machine will use NFLlib and communicate serialized data to a little-endian
  // machine.
  
  /* manual serializers
  */
  void serialize_manually(std::ostream& outputstream) {
    outputstream.write(reinterpret_cast<char*>(_data), N * sizeof(T));
  }
  void deserialize_manually(std::istream& inputstream) {
    inputstream.read(reinterpret_cast<char*>(_data), N * sizeof(T));
  }

  /* serializer (cereal)
  */
  template<class Archive> void serialize(Archive & archive) { 
    archive( _data ); // serialize coefficients by passing them to the archive
  }

  protected:
  // NTT-based Fast Lattice library main class
  // - degree is the degree of the quotient polynomial in the ring, it must
  // be lower or equal than params<T>::kMaxPolyDegree
  // - aggregatedModulusBitsize is the size of the composite modulus to be used,
  // it must be a multiple of params<T>::kModulusBitsize
  // - T is the limb for the submoduli, it must be uint16_t, uint32_t or uint64_t
  class core {

	  template <class P> friend class tests::poly_tests_proxy;

  public:

    /* Constructor
     */
    core();

    /* Number-Theoretic Functions
     */

    void ntt_pow_phi(poly &op);
    void invntt_pow_invphi(poly&);
    static bool ntt(value_type* x, const value_type* wtab, const value_type* winvtab, value_type  const p);
    static bool inv_ntt(value_type *x, const value_type* const inv_wtab, const value_type* const inv_winvtab, value_type invK, value_type const p);

  private:
    // NTT and inversse NTT related attributes
    // omega**i values (omega = polydegree-th primitive root of unity),
    // and their inverses
    // phi**i values (phi = square root of omega), and their inverses
    // multiplied by a constant
    // Shoup variables contain redundant data to speed up the process
    // NOTE : omega and derived values are ordered following David Harvey's
    // algorithm w**0 w**1 ... w**(degree/2) w**0 w**2 ...
    // w**(degree/2) w**0 w**4 ... w**(degree/2) etc. (degree values)
    value_type
      phis[nmoduli][degree] __attribute__((aligned(32))),
      shoupphis[nmoduli][degree]  __attribute__((aligned(32))),
      invpoly_times_invphis[nmoduli][degree]  __attribute__((aligned(32))),
      shoupinvpoly_times_invphis[nmoduli][degree]  __attribute__((aligned(32))),
      omegas[nmoduli][degree * 2]  __attribute__((aligned(32))),
      *shoupomegas[nmoduli]  __attribute__((aligned(32))),
      invomegas[nmoduli][2 * degree]  __attribute__((aligned(32))),
      *shoupinvomegas[nmoduli]  __attribute__((aligned(32))),
      invpolyDegree[nmoduli]  __attribute__((aligned(32)));

    /* Initializing function
     */

    void initialize();
    static void prep_wtab(value_type* wtab, value_type* winvtab, value_type w, size_t cm);

  } __attribute__((aligned(32)));

  static core base;

  protected:
  // NTT-based Fast Lattice library GMP class
  class GMP {

  public:

    mpz_t   moduli_product;
    mpz_t   modulus_shoup;
    size_t bits_in_moduli_product;
    size_t bits_in_modulus_shoup;
    size_t shift_modulus_shoup;
    std::array<mpz_t, nmoduli> lifting_integers;

    /* Constructor & Destructor
     */
    GMP();
    ~GMP();

    /* GMP Functions
     */

    std::array<mpz_t, Degree>  poly2mpz(poly const&);
    void    poly2mpz(std::array<mpz_t, Degree>&, poly const&);
    void    mpz2poly(poly&, std::array<mpz_t, Degree> const&);
  };

  static GMP gmp;

public:
  poly(mpz_t const& v);
  poly(mpz_class const& v);
  poly(std::array<mpz_t, Degree> const& values);
  poly(std::array<mpz_class, Degree> const& values);
  poly(std::initializer_list<mpz_t> const& values);
  poly(std::initializer_list<mpz_class> const& values);
  
  void set_mpz(mpz_t const& v);
  void set_mpz(mpz_class const& v);
  void set_mpz(std::array<mpz_t, Degree> const& values);
  void set_mpz(std::array<mpz_class, Degree> const& values);
  void set_mpz(std::initializer_list<mpz_t> const& values);
  void set_mpz(std::initializer_list<mpz_class> const& values);
  template<class It> void set_mpz(It first, It last);
  
  poly& operator=(mpz_t const& v) { set_mpz(v); return *this; }
  poly& operator=(mpz_class const& v) { set_mpz(v); return *this; }
  poly& operator=(std::array<mpz_t, Degree> const& values) { set_mpz(values); return *this; }
  poly& operator=(std::array<mpz_class, Degree> const& values) { set_mpz(values); return *this; }
  poly& operator=(std::initializer_list<mpz_t> const& values) { set_mpz(values); return *this; }
  poly& operator=(std::initializer_list<mpz_class> const& values) { set_mpz(values); return *this; }

  inline std::array<mpz_t, Degree> poly2mpz() { return gmp.poly2mpz(*this); };
  inline void poly2mpz(std::array<mpz_t, Degree> & array) { gmp.poly2mpz(array, *this); };
  inline void mpz2poly(std::array<mpz_t, Degree> const& array) { gmp.mpz2poly(*this, array); };
  
  inline static constexpr size_t bits_in_moduli_product() { return gmp.bits_in_moduli_product; };
  inline static constexpr mpz_t& moduli_product() { return gmp.moduli_product; };
  inline static constexpr mpz_t& modulus_shoup() { return gmp.modulus_shoup; };
  inline static constexpr std::array<mpz_t, nmoduli> lifting_integers() { return gmp.lifting_integers; };

}  __attribute__((aligned(32)));


/* High level operations on polynomials. They are just wrappers over the operator overloads
 */
template<class T, size_t Degree, size_t NbModuli>
void sub(poly<T, Degree, NbModuli> & out,
         poly<T, Degree, NbModuli> const &arg0, poly<T, Degree, NbModuli> const & arg1) {
  out = arg0 - arg1;
}

/**
 * Perform the modular addition of polynomial @p arg0 with polynomial @p arg1 and store the result in polynomial @p out.
 */
template<class T, size_t Degree, size_t NbModuli>
void add(poly<T, Degree, NbModuli> & out,
         poly<T, Degree, NbModuli> const &arg0, poly<T, Degree, NbModuli> const & arg1) {
  out = arg0 + arg1;
}
template<class T, size_t Degree, size_t NbModuli>
void mul(poly<T, Degree, NbModuli> & out,
         poly<T, Degree, NbModuli> const &arg0, poly<T, Degree, NbModuli> const & arg1) {
  out = arg0 * arg1;
}

/* misc type adaptor
 */
template<class T, size_t Degree, size_t AggregatedModulusBitSize>
using poly_from_modulus = poly<T, Degree, AggregatedModulusBitSize / params<T>::kModulusBitsize>;

/* stream operator
 */
template<class T, size_t Degree, size_t NbModuli>
std::ostream& operator<<(std::ostream& os, nfl::poly<T, Degree, NbModuli> const& p);

/* operator overloads - includes expression templates
 */
DECLARE_BINARY_OPERATOR(operator-, submod)
DECLARE_BINARY_OPERATOR(operator+, addmod)
DECLARE_BINARY_OPERATOR(operator==, eqmod)
DECLARE_BINARY_OPERATOR(operator!=, neqmod)
DECLARE_BINARY_OPERATOR(operator*, mulmod)
DECLARE_BINARY_OPERATOR(shoup, shoup)
DECLARE_UNARY_OPERATOR(compute_shoup, compute_shoup)

}

#include "nfl/core.hpp"
#include "nfl/gmp.hpp"

#endif
