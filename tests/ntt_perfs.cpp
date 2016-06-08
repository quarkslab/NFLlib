#include "nfllib_demo_main.hpp"
#include <chrono>
#include <iostream>
#include "tools.h"

#define REPETITIONS 50000

typedef unsigned int uint128_t __attribute__((mode(TI)));

__attribute__((noinline)) void ntt_new(uint64_t* x, uint64_t* wtab, uint64_t* winvtab,
    unsigned k, uint64_t p)
{
  if (k == 1)
    return;

  // special case
  if (k == 2)
  {
    uint64_t u0 = x[0];
    uint64_t u1 = x[1];
    uint64_t t0 = u0 + u1;
    uint64_t t1 = u0 - u1;
    t0 -= (t0 >= 2*p) ? (2*p) : 0;
    t1 += ((int64_t) t1 < 0) ? (2*p) : 0;
    x[0] = t0;
    x[1] = t1;
    return;
  }

  size_t N = k;           // size of block
  size_t M = 1;                         // number of blocks

  for (; N > 4; N /= 2, M *= 2)
  {
    uint64_t* x0 = x;
    uint64_t* x1 = x + N/2;
    for (size_t r = 0; r < M; r++, x0 += N, x1 += N)
    {
      ptrdiff_t i = N/2 - 2;
      do
      {
        {
          uint64_t u0 = x0[i+1];
          uint64_t u1 = x1[i+1];

          uint64_t t0 = u0 + u1;
          t0 -= ((t0 >= 2*p) ? (2*p) : 0);

          uint64_t t1 = u0 - u1 + 2*p;

          uint64_t q = ((uint128_t) t1 * winvtab[i+1]) >> 64;
          uint64_t t2 = t1 * wtab[i+1] - q * p;

          x0[i+1] = t0;
          x1[i+1] = t2;
        }
        {
          uint64_t u0 = x0[i];
          uint64_t u1 = x1[i];

          uint64_t t0 = u0 + u1;
          t0 -= ((t0 >= 2*p) ? (2*p) : 0);

          uint64_t t1 = u0 - u1 + 2*p;

          uint64_t q = ((uint128_t) t1 * winvtab[i]) >> 64;
          uint64_t t2 = t1 * wtab[i] - q * p;

          x0[i] = t0;
          x1[i] = t2;
        }

        i -= 2;
      }
      while (i >= 0);
    }
    wtab += N/2;
    winvtab += N/2;
  }

  // last two layers
  for (size_t r = 0; r < M; r++, x += 4)
  {
    uint64_t u0 = x[0];
    uint64_t u1 = x[1];
    uint64_t u2 = x[2];
    uint64_t u3 = x[3];

    uint64_t v0 = u0 + u2;
    v0 -= (v0 >= 2*p) ? (2*p) : 0;
    uint64_t v2 = u0 - u2;
    v2 += ((int64_t) v2 < 0) ? (2*p) : 0;

    uint64_t v1 = u1 + u3;
    v1 -= (v1 >= 2*p) ? (2*p) : 0;
    uint64_t t = u1 - u3 + 2*p;

    uint64_t q = ((uint128_t) t * winvtab[1]) >> 64;
    uint64_t v3 = t * wtab[1] - q * p;

    uint64_t z0 = v0 + v1;
    z0 -= (z0 >= 2*p) ? (2*p) : 0;
    uint64_t z1 = v0 - v1;
    z1 += ((int64_t) z1 < 0) ? (2*p) : 0;

    uint64_t z2 = v2 + v3;
    z2 -= (z2 >= 2*p) ? (2*p) : 0;
    uint64_t z3 = v2 - v3;
    z3 += ((int64_t) z3 < 0) ? (2*p) : 0;

    x[0] = z0;
    x[1] = z1;
    x[2] = z2;
    x[3] = z3;
  }

}

namespace nfl { namespace tests {

template <class P>
class poly_tests_proxy
{
  using value_type = typename P::value_type;

  public:
  static inline bool ntt(value_type* x, const value_type* wtab, const value_type* winvtab, value_type const p)
  {
    return P::core::ntt(x, wtab, winvtab, p);
  }

  static inline value_type* get_omegas(P& p) { return &p.base.omegas[0][0]; }
  static inline value_type* get_shoupomegas(P& p) { return &p.base.shoupomegas[0][0]; }
};

} // tests

} // nfl

template<size_t degree, size_t modulus, class T>
int run()
{
  std::cout << std::endl << "Polynomials of degree " << degree << " with " << modulus << " bit coefficients and " << sizeof(T) * 8 << " bit limbs" << std::endl;
  std::cout << "======================================================================" << std::endl;

  using poly_t = nfl::poly_from_modulus<T, degree, modulus>;
  using poly_proxy_t = nfl::tests::poly_tests_proxy<poly_t>;

  auto start = std::chrono::steady_clock::now();
  poly_t *resa = alloc_aligned<poly_t, 32>(REPETITIONS);
  std::fill(resa, resa + REPETITIONS, 0);

  std::fill(resa, resa + REPETITIONS, nfl::uniform());
  start = std::chrono::steady_clock::now();
  for (unsigned i = 0; i < REPETITIONS ; i++)
  {
    poly_t& p = resa[i];
    ntt_new(&p(0, 0), poly_proxy_t::get_omegas(p), poly_proxy_t::get_shoupomegas(p), degree, p.get_modulus(0));
  }
  auto end = std::chrono::steady_clock::now();
  std::cout << "Time per NTT (org): " << get_time_us(start, end, REPETITIONS) << " us" << std::endl;

  std::fill(resa, resa + REPETITIONS, nfl::uniform());
  start = std::chrono::steady_clock::now();
  for (unsigned i = 0; i < REPETITIONS ; i++)
  {
    poly_t& p = resa[i];
    poly_proxy_t::ntt(&p(0, 0), poly_proxy_t::get_omegas(p), poly_proxy_t::get_shoupomegas(p), p.get_modulus(0));
  }
  end = std::chrono::steady_clock::now();
  std::cout << "Time per NTT (lib): " << get_time_us(start, end, REPETITIONS) << " us" << std::endl;

  return 0;
}

int main()
{
  return run<1024, 124, uint64_t>();
}
