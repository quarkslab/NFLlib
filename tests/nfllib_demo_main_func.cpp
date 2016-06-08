#include "nfllib_demo_main.hpp"
#include <chrono>
#include <iostream>
#include "tools.h"

#define REPETITIONS 10
#define NOISE_UB 4
#define NOISE_SK 2
#define LOOP 10

#define TEST_SETCOEFFS
#define TEST_UNIFORM_GENERATION
#define TEST_BOUNDED_GENERATION
#define TEST_GAUSSIAN_GENERATION
#define TEST_NTT
#define TEST_INVNTT
#define TEST_ADDITIONS
#define TEST_ADDITIONS_INPLACE
#define TEST_SUBTRACTIONS
#define TEST_MULTIPLICATIONS
#define TEST_MULTIPLICATIONS_SHOUP
#define TEST_FMAS
#define TEST_FMAS_SHOUP
#define TEST_LWE_SYMMETRIC

template<size_t degree, size_t modulus, class T>
int run()
{
  std::cout << std::endl << "Polynomials of degree " << degree << " with " << modulus << " bit coefficients and " << sizeof(T) * 8 << " bit limbs" << std::endl;
  std::cout << "======================================================================" << std::endl;

  using poly_t = nfl::poly_from_modulus<T, degree, modulus>;
  static_assert(sizeof(poly_t) % 32 == 0, "sizeof(poly_t) must be 32-bytes aligned");

  auto start = std::chrono::steady_clock::now();
  auto end = std::chrono::steady_clock::now();

  // Polynomial arrays to do the tests 
  start = std::chrono::steady_clock::now();
/* AG: on my system, this gives pointers non aligned on 32-bytes!
  poly_t *resa = new poly_t[REPETITIONS],
         *resb = new poly_t[REPETITIONS],
         *resc = new poly_t[REPETITIONS],
         *resd = new poly_t[REPETITIONS];
*/
  poly_t *resa = alloc_aligned<poly_t, 32>(REPETITIONS),
         *resb = alloc_aligned<poly_t, 32>(REPETITIONS),
         *resc = alloc_aligned<poly_t, 32>(REPETITIONS),
         *resd = alloc_aligned<poly_t, 32>(REPETITIONS);
  if ((((uintptr_t)resa % 32) != 0) ||
	  (((uintptr_t)resb % 32) != 0) ||
	  (((uintptr_t)resc % 32) != 0) ||
	  (((uintptr_t)resd % 32) != 0)) {
	  printf("fatal error: pointer unaligned!\n");
	  exit(1);
  }
  std::fill(resa, resa + REPETITIONS, 0);
  std::fill(resb, resb + REPETITIONS, 0);
  std::fill(resc, resc + REPETITIONS, 0);
  std::fill(resd, resd + REPETITIONS, 0);


#ifdef TEST_ADDITIONS_INPLACE
  std::fill(resa, resa + REPETITIONS, nfl::uniform());
  std::fill(resb, resb + REPETITIONS, nfl::uniform());
  start = std::chrono::steady_clock::now();
  for (unsigned i = 0; i < REPETITIONS ; i++)
  {
    nfl::add(resa[i], resa[i], resb[i]);
  }
  end = std::chrono::steady_clock::now();
  std::cout << "Time per polynomial in-place addition a+=b: " << get_time_us(start, end, REPETITIONS) << " us" << std::endl;
#endif

#ifdef TEST_ADDITIONS
  std::fill(resa, resa + REPETITIONS, nfl::uniform());
  std::fill(resb, resb + REPETITIONS, nfl::uniform());
  start = std::chrono::steady_clock::now();
  for (unsigned i = 0; i < REPETITIONS ; i++)
  {
    nfl::add(resc[i], resa[i], resb[i]);
  }
  end = std::chrono::steady_clock::now();
  std::cout << "Time per polynomial addition c=a+b: " << get_time_us(start, end, REPETITIONS) << " us" << std::endl;
#endif

#ifdef TEST_SUBTRACTIONS
  std::fill(resa, resa + REPETITIONS, nfl::uniform());
  std::fill(resb, resb + REPETITIONS, nfl::uniform());
  start = std::chrono::steady_clock::now();
  for (unsigned i = 0; i < REPETITIONS ; i++)
  {
    nfl::sub(resc[i], resa[i], resb[i]);
  }
  end = std::chrono::steady_clock::now();
  std::cout << "Time per polynomial subtraction c=a-b: " << get_time_us(start, end, REPETITIONS) << " us" << std::endl;
#endif

#ifdef TEST_MULTIPLICATIONS
  std::fill(resa, resa + REPETITIONS, nfl::uniform());
  std::fill(resb, resb + REPETITIONS, nfl::uniform());
  start = std::chrono::steady_clock::now();
  for (unsigned i = 0; i < REPETITIONS ; i++)
  {
    nfl::mul(resc[i], resa[i], resb[i]);
  }
  end = std::chrono::steady_clock::now();
  std::cout << "Time per polynomial multiplication (NTT form) c=a*b: " << get_time_us(start, end, REPETITIONS) << " us" << std::endl;
#endif

  // Cleaning
  free_aligned(REPETITIONS, resa);
  free_aligned(REPETITIONS, resb);
  free_aligned(REPETITIONS, resc);
  free_aligned(REPETITIONS, resd);

  return 0;
}

int main() {
  return run<CONFIG>();
}


