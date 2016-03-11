#include <nfl.hpp>
#include "test_binary_op.h"

template <class P>
bool test_op()
{
  using T = typename P::value_type; 
  using greater_value_type = typename P::greater_value_type;

  return test_binary_op<P>(
      [](T a, T b, T m) { return ((greater_value_type)a*b) % m; },
      [](P const& a, P const& b) { return a*b; }
  );
}

template<size_t degree, size_t modulus, class T>
bool run()
{
  using poly_t = nfl::poly_from_modulus<T, degree, modulus>;
  using poly_p_t = nfl::poly_p_from_modulus<T, degree, modulus>;

  return test_op<poly_t>() and test_op<poly_p_t>();
}

int main(int argc, char const* argv[]) {
  return not run<CONFIG>() ;
}
