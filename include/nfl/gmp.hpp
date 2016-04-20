#ifndef NFL_GMP_HPP
#define NFL_GMP_HPP

#include <vector>
#include "nfl/poly.hpp"

namespace nfl {

template <class T, size_t Degree, size_t NbModuli>
typename poly<T, Degree, NbModuli>::GMP poly<T, Degree, NbModuli>::gmp;

/*
 * Constructors
 */
template <class T, size_t Degree, size_t NbModuli>
poly<T, Degree, NbModuli>::poly(mpz_t const& v) {
  set_mpz(v);
}

template <class T, size_t Degree, size_t NbModuli>
poly<T, Degree, NbModuli>::poly(mpz_class const& v) {
  set_mpz({v}, (size_t)1);
}

template <class T, size_t Degree, size_t NbModuli>
poly<T, Degree, NbModuli>::poly(std::initializer_list<mpz_class> const& values, size_t length) {
  set_mpz(values, length);
}

template <class T, size_t Degree, size_t NbModuli>
poly<T, Degree, NbModuli>::poly(mpz_t* const& values, size_t length) {
  set_mpz(values, length);
}

template <class T, size_t Degree, size_t NbModuli>
poly<T, Degree, NbModuli>::poly(mpz_class* const& values, size_t length) {
  set_mpz(values, length);
}

/*
 * set_mpz's functions
 */

template <class T, size_t Degree, size_t NbModuli>
void poly<T, Degree, NbModuli>::set_mpz(mpz_t const& v) {
  set_mpz({mpz_class(v)}, (size_t)1);
}

template <class T, size_t Degree, size_t NbModuli>
void poly<T, Degree, NbModuli>::set_mpz(mpz_class const& v) {
  set_mpz({v}, (size_t)1);
}

template <class T, size_t Degree, size_t NbModuli>
void poly<T, Degree, NbModuli>::set_mpz(mpz_t *const &values, size_t length) {
  std::vector<mpz_class> v(length);
  for (size_t i = 0; i < length; i++) {
    v[i] = mpz_class(values[i]);
  }
  set_mpz(v.begin(), v.end(), length);
}

template <class T, size_t Degree, size_t NbModuli>
void poly<T, Degree, NbModuli>::set_mpz(mpz_class* const& values, size_t length) {
  std::vector<mpz_class> v(length);
  for (size_t i = 0; i < length; i++) {
    v[i] = values[i];
  }
  set_mpz(v.begin(), v.end(), length);
}

template <class T, size_t Degree, size_t NbModuli>
void poly<T, Degree, NbModuli>::set_mpz(
    std::initializer_list<mpz_class> const& values, size_t length) {
  set_mpz(values.begin(), values.end(), length);
}

template <class T, size_t Degree, size_t NbModuli>
template <class It>
void poly<T, Degree, NbModuli>::set_mpz(It first, It last, size_t length) {
  // CRITICAL: the object must be 32-bytes aligned to avoid vectorization issues
  assert((unsigned long)(this->_data) % 32 == 0);

  auto* iter = begin();
  auto viter = first;

  size_t size = (!length) ? std::distance(first, last) : length;
  // If the initializer has no more values than the polynomial degree use them
  // to initialize the associated coefficients for each sub-modulus
  if (size <= degree || size == degree * nmoduli) {
    for (size_t cm = 0; cm < nmoduli; ++cm) {
      viter = (size != degree * nmoduli) ? first : viter;

      for (size_t i = 0; i < degree; i++) {
        if (viter < last) {
          *iter++ = mpz_fdiv_ui(viter->get_mpz_t(), params<T>::P[cm]);
          viter++;
        } else
          *iter++ = 0;
      }
    }
  } else {
    throw std::runtime_error(
        "gmp: CRITICAL, initializer of size above degree but different "
        "from nmoduli*degree");
  }
}

/*
 * Nested GMP class
 */
template <class T, size_t Degree, size_t NbModuli>
poly<T, Degree, NbModuli>::GMP::GMP() {
  // Compute product of moduli
  mpz_init_set_ui(moduli_product, 1);
  for (size_t i = 0; i < nmoduli; i++) {
    mpz_mul_ui(moduli_product, moduli_product, params<T>::P[i]);
  }

  bits_in_moduli_product = mpz_sizeinbase(moduli_product, 2);

  // Compute Shoup value for optimized reduction modulo "moduli_product"
  auto const shift_modulus_shoup = bits_in_moduli_product +
                                   sizeof(T) * CHAR_BIT +
                                   static_log2<nmoduli>::value + 1;
  mpz_t one;
  mpz_init_set_ui(one, 1);
  mpz_init(modulus_shoup);
  mpz_mul_2exp(modulus_shoup, one, shift_modulus_shoup);
  mpz_tdiv_q(modulus_shoup, modulus_shoup, moduli_product);

  // Compute the lifting coefficients
  mpz_t quotient;
  mpz_init(quotient);
  lifting_integers = (mpz_t*)malloc(nmoduli * sizeof(mpz_t));

  for (size_t i = 0; i < nmoduli; i++) {
    // Current modulus
    mpz_t current_modulus;
    mpz_init_set_ui(current_modulus, params<T>::P[i]);

    // compute the product of primes except the current one
    mpz_divexact(quotient, moduli_product, current_modulus);

    // Compute the inverse of the product
    mpz_init(lifting_integers[i]);
    mpz_invert(lifting_integers[i], quotient, current_modulus);

    // Multiply by the quotient
    mpz_mul(lifting_integers[i], lifting_integers[i], quotient);
    mpz_clear(current_modulus);
  }

  // Clear
  mpz_clear(one);
  mpz_clear(quotient);
}

template <class T, size_t Degree, size_t NbModuli>
poly<T, Degree, NbModuli>::GMP::~GMP() {
  for (size_t i = 0; i < nmoduli; i++) {
    mpz_clear(lifting_integers[i]);
  }
  free(lifting_integers);
  mpz_clear(modulus_shoup);
  mpz_clear(moduli_product);
}

/**
 * Functions to convert a poly into an array of mpz_t, and back
 *
 * The mpz_t array is allocated with new[]
 */

template <class T, size_t Degree, size_t NbModuli>
mpz_t* poly<T, Degree, NbModuli>::GMP::poly2mpz(poly const& op) {
  // Assign and init
  mpz_t* poly_mpz = new mpz_t[degree];
  for (size_t i = 0; i < degree; i++) {
    mpz_init2(poly_mpz[i], bits_in_moduli_product);
  }

  // Fill
  poly2mpz(poly_mpz, op);
  return poly_mpz;
}

template <class T, size_t Degree, size_t NbModuli>
void poly<T, Degree, NbModuli>::GMP::poly2mpz(mpz_t* rop, poly const& op) {

  auto const shift_modulus_shoup = bits_in_moduli_product +
                                   sizeof(T) * CHAR_BIT +
                                   static_log2<nmoduli>::value + 1;

  mpz_t tmp;
  mpz_init2(tmp, shift_modulus_shoup + static_log2<nmoduli>::value);

  // Loop on each coefficient
  for (size_t i = 0; i < degree; i++) {
    mpz_set_ui(rop[i], 0);
    for (size_t cm = 0; cm < nmoduli; cm++) {
      if (op(cm, i) != 0) {
        mpz_mul_ui(tmp, lifting_integers[cm], op(cm, i));
      }
      mpz_add(rop[i], rop[i], tmp);
    }
    // Modular reduction using Shoup
    mpz_mul(tmp, rop[i], modulus_shoup);
    mpz_tdiv_q_2exp(tmp, tmp, shift_modulus_shoup);
    mpz_mul(tmp, tmp, moduli_product);
    mpz_sub(rop[i], rop[i], tmp);
    if (mpz_cmp(rop[i], moduli_product) >= 0) {
      mpz_sub(rop[i], rop[i], moduli_product);
    }
  }

  // Clean
  mpz_clear(tmp);
}

template <class T, size_t Degree, size_t NbModuli>
void poly<T, Degree, NbModuli>::GMP::mpz2poly(poly<T, Degree, NbModuli>& rop,
                                              mpz_t* const& poly_mpz) {
  for (size_t cm = 0; cm < nmoduli; cm++) {
    for (size_t i = 0; i < degree; i++) {
      rop(cm, i) = mpz_fdiv_ui(poly_mpz[i], params<T>::P[cm]);
    }
  }
}
}

#endif