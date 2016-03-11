#include <nfl.hpp>
#include <sstream>

#include "tools.h"

template<size_t degree, size_t modulus, class T>
bool run()
{
  using poly_t = nfl::poly_from_modulus<T, degree, modulus>;
  poly_t& p = *alloc_aligned<poly_t, 32>(1);
  p = 1;
  std::ostringstream oss;
  oss << p << std::endl;
  free_aligned(1, &p);
  return oss.str().substr(0, 4) == "{ 1U";
}

template<size_t degree, size_t modulus, class T>
bool run_p()
{
  using poly_p = nfl::poly_p_from_modulus<T, degree, modulus>;

  poly_p p;
  p = 1;
  std::ostringstream oss;
  oss << p << std::endl;
  return oss.str().substr(0, 4) == "{ 1U";
}

int main(int argc, char const* argv[]) {
  return not (run<CONFIG>() and run_p<CONFIG>());
}
