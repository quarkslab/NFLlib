#include <iostream>
#include <nfl.hpp>
#include <string>
#include <vector>
#include "tools.h"

template <size_t degree, size_t modulus, class T>
bool run() {
  using poly_t = nfl::poly_from_modulus<T, degree, modulus>;

  bool ret_value = true;

  // Set a random number generator
  // (we don't need a cryptographically secure random number generator)
  gmp_randclass prng(gmp_randinit_default);
  prng.seed(0);

  // define a random polynomial
  poly_t& p0 = *alloc_aligned<poly_t, 32>(1, nfl::uniform());

  // get the corresponding array of mpz_t
  std::array<mpz_t, degree> coefficients = p0.poly2mpz();

  // construct a polynomial from mpz_t coefficients
  poly_t& p1 = *alloc_aligned<poly_t, 32>(1, 0);
  p1.mpz2poly(coefficients);

  // verify that the first and second polynomials are equal
  ret_value &= (p0 == p1);

  // construct a vector of mpz_class
  std::vector<mpz_class> coefficients_mpz_class(poly_t::degree);
  for (size_t i = 0; i < poly_t::degree; i++) {
    coefficients_mpz_class[i] = mpz_class(coefficients[i]);
  }

  // construct a polynomial from mpz_class coefficients
  poly_t& p2 = *alloc_aligned<poly_t, 32>(1, 0);
  p2.set_mpz(coefficients_mpz_class.begin(), coefficients_mpz_class.end());

  // verify that the first and third polynomials are equal
  ret_value &= (p0 == p2);

  // Generate a random vector of mpz_class
  std::vector<mpz_class> big_array;
  for (size_t i = 0; i < poly_t::nmoduli * poly_t::degree; i++) {
    // (say) 200 bits, it doesn't matter as long as it is > 64 bits for the
    // tests
    big_array.push_back(prng.get_z_bits(200));
  }

  // construct a polynomial from the big_array vector
  poly_t& p3 = *alloc_aligned<poly_t, 32>(1, 0);
  p3.set_mpz(big_array.begin(), big_array.end());

  // verify that the coefficients of the polynomial have been set correctly
  for (size_t cm = 0; cm < poly_t::nmoduli; cm++) {
    mpz_class modulus_mpz_class(std::to_string(poly_t::get_modulus(cm)));
    for (size_t i = 0; i < poly_t::degree; i++) {
      mpz_class coeff(std::to_string(p3(cm, i)));
      ret_value &=
          ((big_array[cm * poly_t::degree + i] % modulus_mpz_class) == coeff);
    }
  }

  // Cleaning
  for (size_t i = 0; i < poly_t::degree; i++) {
    mpz_clear(coefficients[i]);
  }
  free_aligned(1, &p0);
  free_aligned(1, &p1);
  free_aligned(1, &p2);
  free_aligned(1, &p3);

  return ret_value;
}

template <size_t degree, size_t modulus, class T>
bool run_p() {
  using poly_p = nfl::poly_p_from_modulus<T, degree, modulus>;

  bool ret_value = true;

  // Set a random number generator
  // (we don't need a cryptographically secure random number generator)
  gmp_randclass prng(gmp_randinit_default);
  prng.seed(0);

  // define a random polynomial
  poly_p p0{nfl::uniform()};

  // get the corresponding array of mpz_t
  std::array<mpz_t, degree> coefficients = p0.poly2mpz();

  // construct a polynomial from mpz_t coefficients
  poly_p p1;
  p1.mpz2poly(coefficients);

  // verify that the first and second polynomials are equal
  ret_value &= (p0 == p1);

  // construct a vector of mpz_class
  std::vector<mpz_class> coefficients_mpz_class(poly_p::degree);
  for (size_t i = 0; i < poly_p::degree; i++) {
    coefficients_mpz_class[i] = mpz_class(coefficients[i]);
  }

  // construct a polynomial from mpz_class coefficients
  poly_p p2;
  p2.set_mpz(coefficients_mpz_class.begin(), coefficients_mpz_class.end());

  // verify that the first and third polynomials are equal
  ret_value &= (p0 == p2);

  // Generate a random vector of mpz_class
  std::vector<mpz_class> big_array;
  for (size_t i = 0; i < poly_p::nmoduli * poly_p::degree; i++) {
    // (say) 200 bits, it doesn't matter as long as it is > 64 bits for the
    // tests
    big_array.push_back(prng.get_z_bits(200));
  }

  // construct a polynomial from the big_array vector
  poly_p p3;
  p3.set_mpz(big_array.begin(), big_array.end());

  // verify that the coefficients of the polynomial have been set correctly
  for (size_t cm = 0; cm < poly_p::nmoduli; cm++) {
    mpz_class modulus_mpz_class(std::to_string(poly_p::get_modulus(cm)));
    for (size_t i = 0; i < poly_p::degree; i++) {
      mpz_class coeff(std::to_string(p3(cm, i)));
      ret_value &=
          ((big_array[cm * poly_p::degree + i] % modulus_mpz_class) == coeff);
    }
  }

  // Cleaning
  for (size_t i = 0; i < poly_p::degree; i++) {
    mpz_clear(coefficients[i]);
  }

  return ret_value;
}

int main(int argc, char const* argv[]) {
  return not(run<CONFIG>() and run_p<CONFIG>());
}
