#include <iostream>
#include <nfl.hpp>
#include "tools.h"
#include <vector>

template <size_t degree, size_t modulus, class T>
bool run() {
  using poly_t = nfl::poly_from_modulus<T, degree, modulus>;

  bool ret_value = true;

  // Set a random number generator
  // (we don't need a cryptographically secure random number generator)
  std::srand(std::time(0));

  // define a random polynomial
  poly_t& p0 = *alloc_aligned<poly_t, 32>(1, nfl::uniform());

  // get all the coefficients in an array
  std::vector<T> big_array(poly_t::degree * poly_t::nmoduli);
  for (size_t cm = 0; cm < poly_t::nmoduli; cm++) {
    for (size_t i = 0; i < poly_t::degree; i++) {
      big_array[cm * poly_t::degree + i] = p0(cm, i);
    }
  }

  // define a polynomial from the array
  poly_t& p1 =
      *alloc_aligned<poly_t, 32>(1, big_array.begin(), big_array.end());

  // verify that the first and second polynomials are equal
  ret_value &= (p0 == p1);

  // Rerandomize the array
  std::generate(big_array.begin(), big_array.end(), std::rand);

  // define a polynomial from the array
  poly_t& p2 =
      *alloc_aligned<poly_t, 32>(1, big_array.begin(), big_array.end());

  // verify that the coefficients of the polynomial have been set correctly
  for (size_t cm = 0; cm < poly_t::nmoduli; cm++) {
    for (size_t i = 0; i < poly_t::degree; i++) {
      ret_value &=
          (p2(cm, i) ==
           big_array[cm * poly_t::degree + i] % poly_t::get_modulus(cm));
    }
  }

  // define a polynomial from the array without reducing the coefficients
  // reduce_coeffs = false
  poly_t& p3 =
      *alloc_aligned<poly_t, 32>(1, big_array.begin(), big_array.end(), false);

  // verify that the coefficients of the polynomial have been set correctly
  for (size_t cm = 0; cm < poly_t::nmoduli; cm++) {
    for (size_t i = 0; i < poly_t::degree; i++) {
      ret_value &= (p3(cm, i) == big_array[cm * poly_t::degree + i]);
    }
  }

  // construct an array of "degree" random coefficients
  std::vector<T> small_array(poly_t::degree);
  std::generate(small_array.begin(), small_array.end(), std::rand);

  // define a polynomial from the array
  poly_t& p4 =
      *alloc_aligned<poly_t, 32>(1, small_array.begin(), small_array.end());

  // verify that the coefficients of the polynomial have been set correctly
  for (size_t cm = 0; cm < poly_t::nmoduli; cm++) {
    for (size_t i = 0; i < poly_t::degree; i++) {
      ret_value &= (p4(cm, i) == small_array[i] % poly_t::get_modulus(cm));
    }
  }

  // define a polynomial from the array without reducing the coefficients
  // reduce_coeffs = false
  poly_t& p5 = *alloc_aligned<poly_t, 32>(1, small_array.begin(),
                                          small_array.end(), false);

  // verify that the coefficients of the polynomial have been set correctly
  for (size_t cm = 0; cm < poly_t::nmoduli; cm++) {
    for (size_t i = 0; i < poly_t::degree; i++) {
      ret_value &= (p5(cm, i) == small_array[i]);
    }
  }

  // define the zero polynomial
  poly_t& p6 = *alloc_aligned<poly_t, 32>(1, 0);

  // verify that the coefficients of the polynomial have been set correctly
  for (size_t cm = 0; cm < poly_t::nmoduli; cm++) {
    for (size_t i = 0; i < poly_t::degree; i++) {
      ret_value &= (p6(cm, i) == 0);
    }
  }

  // Cleaning
  free_aligned(1, &p0);
  free_aligned(1, &p1);
  free_aligned(1, &p2);
  free_aligned(1, &p3);
  free_aligned(1, &p4);
  free_aligned(1, &p5);
  free_aligned(1, &p6);

  return ret_value;
}

template <size_t degree, size_t modulus, class T>
bool run_p() {
  using poly_p = nfl::poly_p_from_modulus<T, degree, modulus>;

  bool ret_value = true;

  // Set a random number generator
  // (we don't need a cryptographically secure random number generator)
  std::srand(std::time(0));

  // define a random polynomial
  poly_p p0{nfl::uniform()};

  // get all the coefficients in an array
  std::vector<T> big_array(poly_p::degree * poly_p::nmoduli);
  for (size_t cm = 0; cm < poly_p::nmoduli; cm++) {
    for (size_t i = 0; i < poly_p::degree; i++) {
      big_array[cm * poly_p::degree + i] = p0(cm, i);
    }
  }

  // define a polynomial from the array
  poly_p p1{big_array.begin(), big_array.end()};

  // verify that the first and second polynomials are equal
  ret_value &= (p0 == p1);

  // Rerandomize the array
  std::generate(big_array.begin(), big_array.end(), std::rand);

  // define a polynomial from the array
  poly_p p2{big_array.begin(), big_array.end()};

  // verify that the coefficients of the polynomial have been set correctly
  for (size_t cm = 0; cm < poly_p::nmoduli; cm++) {
    for (size_t i = 0; i < poly_p::degree; i++) {
      ret_value &=
          (p2(cm, i) ==
           big_array[cm * poly_p::degree + i] % poly_p::get_modulus(cm));
    }
  }

  // define a polynomial from the array without reducing the coefficients
  // reduce_coeffs = false
  poly_p p3{big_array.begin(), big_array.end(), false};

  // verify that the coefficients of the polynomial have been set correctly
  for (size_t cm = 0; cm < poly_p::nmoduli; cm++) {
    for (size_t i = 0; i < poly_p::degree; i++) {
      ret_value &= (p3(cm, i) == big_array[cm * poly_p::degree + i]);
    }
  }

  // construct an array of "degree" random coefficients
  std::vector<T> small_array(poly_p::degree);
  std::generate(small_array.begin(), small_array.end(), std::rand);

  // define a polynomial from the array
  poly_p p4{small_array.begin(), small_array.end()};

  // verify that the coefficients of the polynomial have been set correctly
  for (size_t cm = 0; cm < poly_p::nmoduli; cm++) {
    for (size_t i = 0; i < poly_p::degree; i++) {
      ret_value &= (p4(cm, i) == small_array[i] % poly_p::get_modulus(cm));
    }
  }

  // define a polynomial from the array without reducing the coefficients
  // reduce_coeffs = false
  poly_p p5{small_array.begin(), small_array.end(), false};

  // verify that the coefficients of the polynomial have been set correctly
  for (size_t cm = 0; cm < poly_p::nmoduli; cm++) {
    for (size_t i = 0; i < poly_p::degree; i++) {
      ret_value &= (p5(cm, i) == small_array[i]);
    }
  }

  // define the zero polynomial
  poly_p p6{0};

  // verify that the coefficients of the polynomial have been set correctly
  for (size_t cm = 0; cm < poly_p::nmoduli; cm++) {
    for (size_t i = 0; i < poly_p::degree; i++) {
      ret_value &= (p6(cm, i) == 0);
    }
  }

  return ret_value;
}

int main(int argc, char const* argv[]) {
  return not(run<CONFIG>() and run_p<CONFIG>());
}
