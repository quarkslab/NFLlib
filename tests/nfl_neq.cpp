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
  // Simple unitary test
  // AG: commented because it uses stack variables!
  //ret_value &= (poly_t{{4,3,2,1}} != poly_t{{5,3,2,1}});
  
  // Some more randomized tests
  // if resb[i] is the null polynomial this test fails but 
  // then the uniform generator is wrong or we are extremely unlucky
  for (size_t i = 0 ; i < ITERATIONS; i++)
    ret_value &= (resa[i] != (resa[i]+resb[i]));

  free_aligned(ITERATIONS, resa);
  free_aligned(ITERATIONS, resb);

  return ret_value;
}

int main(int argc, char const* argv[]) {
  return not run<CONFIG>() ;
}
