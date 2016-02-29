#include <iostream>
#include <nfl/poly_p.hpp>
#include "tools.h"

template<size_t degree, size_t modulus, class T>
bool run()
{
  using poly_p = nfl::poly_p<T, degree, modulus>;
  using poly_t = nfl::poly<T, degree, modulus>;

  // Test simple addition
  poly_p a{nfl::uniform()};
  poly_p b{nfl::uniform()};

  poly_t &tmp = *alloc_aligned<poly_t, 32>(1, a.poly_obj()+b.poly_obj());
  poly_p add_p{a+b};

  bool ret = true;
  ret &= (add_p == tmp);

  // Test detach
  poly_p c{b};
  ret &= (c == b);

  c = {1};
  ret &= (c != b);

  // Shoup
  poly_p bshoup = nfl::compute_shoup(b);
  tmp = nfl::compute_shoup(b.poly_obj());
  ret &= (bshoup == tmp);

  poly_p res = nfl::shoup(a * b, bshoup);
  tmp = nfl::shoup(a.poly_obj() * b.poly_obj(), bshoup.poly_obj());
  ret &= (res == tmp);

  free_aligned(1, &tmp);

  return ret;
}

int main(int argc, char const* argv[]) {
  return not run<CONFIG>() ;
}
