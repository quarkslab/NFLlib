#include <nfl.hpp>
#include "tools.h"

#define ITERATIONS 100

template<size_t degree, size_t modulus, class T>
bool run()
{
  using poly_t = nfl::poly_from_modulus<T, degree, modulus>;

  poly_t* resa = alloc_aligned<poly_t, 32>(ITERATIONS, nfl::uniform());
  
  bool ret_value = true;
  // Simple unitary test
  // AG: commented because it uses stack variables!
  //ret_value &= (poly_t{{1,2,3,4}} == poly_t{{1,2,3,4}}); 
  //ret_value &= (poly_t{{1,2,3,4}} == poly_t{{2,12,34,56}}); 
  // Some more randomized tests
  for (size_t i = 0 ; i < ITERATIONS; i++)
    ret_value &= (resa[i] == resa[i]);

  free_aligned(ITERATIONS, resa);

  return ret_value;
}

int main(int argc, char const* argv[]) {
  return not run<CONFIG>() ;
}
