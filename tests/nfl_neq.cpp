#include <nfl.hpp>
#include "tools.h"

#define ITERATIONS 100


template<size_t degree, size_t modulus, class T>
bool run()
{
  using poly_t = nfl::poly_from_modulus<T, degree, modulus>;

  poly_t* resa = alloc_aligned<poly_t, 32>(ITERATIONS, nfl::uniform());
  poly_t* resb = alloc_aligned<poly_t, 32>(ITERATIONS, nfl::uniform());
  
  bool ret_value = true;
  // Randomized test
  // if resb[i] is the null polynomial this test fails but 
  // then the uniform generator is wrong or we are extremely unlucky
  for (size_t i = 0 ; i < ITERATIONS; i++)
    ret_value &= (resa[i] != (resa[i]+resb[i]));

  free_aligned(ITERATIONS, resa);
  free_aligned(ITERATIONS, resb);

  return ret_value;
}

template<size_t degree, size_t modulus, class T>
bool run_p()
{
  using poly_p = nfl::poly_p_from_modulus<T, degree, modulus>;

  poly_p* resa = new poly_p[ITERATIONS];
  poly_p* resb = new poly_p[ITERATIONS];
  for (size_t i = 0 ; i < ITERATIONS; i++) {
    resa[i] = nfl::uniform();
    resb[i] = nfl::uniform();
  }
  
  bool ret_value = true;
  // Randomized test
  // if resb[i] is the null polynomial this test fails but 
  // then the uniform generator is wrong or we are extremely unlucky
  for (size_t i = 0 ; i < ITERATIONS; i++)
    ret_value &= (resa[i] != (resa[i]+resb[i]));

  delete[] resa;
  delete[] resb;

  return ret_value;
}

int main(int argc, char const* argv[]) {
  return not (run<CONFIG>() and run_p<CONFIG>()) ;
}
