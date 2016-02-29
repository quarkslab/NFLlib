#ifndef NFL_TEST_BINARY_OP_H
#define NFL_TEST_BINARY_OP_H

#include <nfl.hpp>
#include "tools.h"

#define ITERATIONS 2

template <class poly_t, class OpPoly, class OpValue>
bool test_binary_op(OpValue const& op_value, OpPoly const& op_poly)
{
  size_t nmoduli = poly_t::nmoduli;

  poly_t* resa = alloc_aligned<poly_t, 32>(ITERATIONS, nfl::uniform());
  poly_t* resb = alloc_aligned<poly_t, 32>(ITERATIONS, nfl::uniform());
 
  bool ret_value = true;
  poly_t& tmp = *alloc_aligned<poly_t, 32>(1);
  for (size_t i = 0 ; i < ITERATIONS; i++)
  {
    for (size_t cm = 0; cm < nmoduli; cm++)
      for (size_t j = 0 ; j < poly_t::degree ; j++)
        tmp(cm, j) = op_value(resa[i](cm, j), resb[i](cm, j), poly_t::get_modulus(cm));
    ret_value &= (tmp == op_poly(resa[i], resb[i]));
  }

  free_aligned(ITERATIONS, resa);
  free_aligned(ITERATIONS, resb);
  free_aligned(1, &tmp);

  return ret_value;
}

#endif
