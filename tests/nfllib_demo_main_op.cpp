#include "nfllib_demo_main.hpp"
#include <chrono>
#include <iostream>
#include "tools.h"

#define REPETITIONS 10
#define NOISE_UB 4
#define SIGMA 4
#define LOOP 1

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

template <class P>
__attribute__((noinline)) static void encrypt(P& resa, P& resb, P const & pka, P const & pkb, P const & pkaprime, P const & pkbprime, nfl::FastGaussianNoise<uint8_t, typename P::value_type, 2> *g_prng)
{
  // u
  P &tmpu = *alloc_aligned<P, 32>(1, nfl::gaussian<uint8_t, typename P::value_type, 2>(g_prng));
  tmpu.ntt_pow_phi();
  
  // 2*e_1
  P &tmpe1 = *alloc_aligned<P, 32>(1, nfl::gaussian<uint8_t, typename P::value_type, 2>(g_prng, 2));
  tmpe1.ntt_pow_phi();

  // 2*e_2
  P &tmpe2 = *alloc_aligned<P, 32>(1, nfl::gaussian<uint8_t, typename P::value_type, 2>(g_prng, 2));
  tmpe2.ntt_pow_phi();

  // resa = pka * u + 2*e_2
  resa = tmpu * pka + tmpe1;

  // resb = pkb * u + 2*e_2
  resb = tmpu * pkb + tmpe2;
}

template <class P>
__attribute__((noinline)) static void decrypt(P& tmp, P const & resa, P const& resb, P const& s, P const& sprime, typename P::value_type const modulus)
{
  tmp = resb - resa * s;
  
  tmp.invntt_pow_invphi();
  for(auto & v : tmp)
  {
    v= (v<modulus/2) ? v%2 : 1-v%2;
  }
}

template <class P>
bool test_mulmod_shoup(P* resa, P* resb, P* resc, P* resd)
{
  typedef P poly_t;
  std::fill(resa, resa + REPETITIONS, nfl::uniform());
  std::fill(resb, resb + REPETITIONS, nfl::uniform());
  for (unsigned i = 0; i < REPETITIONS ; i++)
  {
    resd[i] = nfl::compute_shoup(resb[i]);
  }
  auto start = std::chrono::steady_clock::now();
  for (unsigned i = 0; i < REPETITIONS ; i++)
  {
    resc[i] = nfl::shoup(resa[i] * resb[i], resd[i]);
  }
  auto end = std::chrono::steady_clock::now();
  poly_t& ptmp = *alloc_aligned<poly_t, 32>(1);
  for (unsigned i = 0; i < REPETITIONS ; i++)
  {
    ptmp = nfl::ops::make_op<nfl::ops::mulmod_shoup<typename P::value_type, nfl::simd::serial>>(resa[i], resb[i], resd[i]);
    if (resc[i] != ptmp) {
      std::cerr << "error with vectorized mul_shoup" << std::endl;
	  return false;
    }
  }
  std::cout << "Time per polynomial multiplication (NTT form with Shoup) c=a*b: " << get_time_us(start, end, REPETITIONS) << " us" << std::endl;
  return true;
}

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

#ifdef TEST_SETCOEFFS
  end = std::chrono::steady_clock::now();
  std::cout << "Time per polynomial setcoeffs: " << get_time_us(start, end, 4*REPETITIONS) << " us" << std::endl;
#endif

#ifdef TEST_UNIFORM_GENERATION
  start = std::chrono::steady_clock::now();
  std::fill(resa, resa + REPETITIONS, nfl::uniform());
  std::fill(resb, resb + REPETITIONS, nfl::uniform());
  end = std::chrono::steady_clock::now();
  std::cout << "Time per polynomial generation (uniform): " << get_time_us(start, end, 2*REPETITIONS) << " us" << std::endl;
#endif

#ifdef TEST_BOUNDED_GENERATION
  start = std::chrono::steady_clock::now();
  std::fill(resa, resa + REPETITIONS, nfl::non_uniform(NOISE_UB));
  std::fill(resb, resb + REPETITIONS, nfl::non_uniform(NOISE_UB));
  end = std::chrono::steady_clock::now();
  std::cout << "Time per polynomial generation (bounded 6 bits): " << get_time_us(start, end, 2*REPETITIONS) << " us" << std::endl;
#endif

#ifdef TEST_GAUSSIAN_GENERATION
  nfl::FastGaussianNoise<uint8_t, T, 2> fg_prng(20, 128, 1<<14);
  start = std::chrono::steady_clock::now();
  std::fill(resa, resa + REPETITIONS, nfl::gaussian<uint8_t, T, 2>(&fg_prng));
  std::fill(resb, resb + REPETITIONS, nfl::gaussian<uint8_t, T, 2>(&fg_prng));
  end = std::chrono::steady_clock::now();
  std::cout << "Time per polynomial generation (gaussian sigma=20 k=128): " << get_time_us(start, end, 2*REPETITIONS) << " us" << std::endl;
#endif

#ifdef TEST_NTT
  std::fill(resa, resa + REPETITIONS, nfl::uniform());
  std::fill(resb, resb + REPETITIONS, nfl::uniform());
  start = std::chrono::steady_clock::now();
  for (unsigned j = 0 ; j < LOOP; j++) for (unsigned i = 0; i < REPETITIONS ; i++)
  {
    resa[i].ntt_pow_phi();
  }
  end = std::chrono::steady_clock::now();
  std::cout << "Time per polynomial NTT: " << get_time_us(start, end, LOOP*REPETITIONS) << " us" << std::endl;
#endif

#ifdef TEST_INVNTT
  std::fill(resa, resa + REPETITIONS, nfl::uniform());
  std::fill(resb, resb + REPETITIONS, nfl::uniform());
  start = std::chrono::steady_clock::now();
  for (unsigned i = 0; i < REPETITIONS ; i++)
  {
    resa[i].invntt_pow_invphi();
  }
  end = std::chrono::steady_clock::now();
  std::cout << "Time per polynomial inverse NTT: " << get_time_us(start, end, REPETITIONS) << " us" << std::endl;
#endif

#ifdef TEST_ADDITIONS_INPLACE
  std::fill(resa, resa + REPETITIONS, nfl::uniform());
  std::fill(resb, resb + REPETITIONS, nfl::uniform());
  start = std::chrono::steady_clock::now();
  for (unsigned i = 0; i < REPETITIONS ; i++)
  {
    resa[i] = resa[i] + resb[i];
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
    resc[i] = resa[i] + resb[i];
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
    resc[i] = resa[i] - resb[i];
  }
  end = std::chrono::steady_clock::now();
  std::cout << "Time per polynomial subtraction c=a-b: " << get_time_us(start, end, REPETITIONS) << " us" << std::endl;
#endif

#ifdef TEST_MULTIPLICATIONS
  std::fill(resa, resa + REPETITIONS, nfl::uniform());
  std::fill(resb, resb + REPETITIONS, nfl::uniform());
  for (unsigned i = 0; i < REPETITIONS ; i++)
  {
    resc[i] = resa[i] * resb[i];
  }
  start = std::chrono::steady_clock::now();
  for (unsigned i = 0; i < REPETITIONS ; i++)
  {
    resc[i] = resa[i] * resb[i];
  }
  end = std::chrono::steady_clock::now();
  std::cout << "Time per polynomial multiplication (NTT form) c=a*b: " << get_time_us(start, end, REPETITIONS) << " us" << std::endl;
#endif

#ifdef TEST_MULTIPLICATIONS_SHOUP
  if (!test_mulmod_shoup(resa, resb, resc, resd)) {
	  return 1;
  }
#endif


#ifdef TEST_FMAS
  std::fill(resa, resa + REPETITIONS, nfl::uniform());
  std::fill(resb, resb + REPETITIONS, nfl::uniform());
  start = std::chrono::steady_clock::now();
  for (unsigned i = 0; i < REPETITIONS ; i++)
  {
    resc[i] = resa[i] * resb[i] + resc[i];
  }
  end = std::chrono::steady_clock::now();
  std::cout << "Time per polynomial FMA (NTT form) c+=a*b: " << get_time_us(start, end, REPETITIONS) << " us" << std::endl;
#endif

#ifdef TEST_FMAS_SHOUP
  std::fill(resa, resa + REPETITIONS, nfl::uniform());
  std::fill(resb, resb + REPETITIONS, nfl::uniform());
  for (unsigned i = 0; i < REPETITIONS ; i++)
  {
    resd[i] = nfl::compute_shoup(resb[i]);
  }
  start = std::chrono::steady_clock::now();
  for (unsigned i = 0; i < REPETITIONS ; i++)
  {
    resc[i] = resc[i] + nfl::shoup(resa[i] * resb[i], resd[i]);
  }
  end = std::chrono::steady_clock::now();
  std::cout << "Time per polynomial FMA (NTT form with Shoup) c+=a*b: " << get_time_us(start, end, REPETITIONS) << " us" << std::endl;
#endif

#ifdef TEST_LWE_SYMMETRIC
  // Get some internal values
  unsigned int polyDegree = degree;
  unsigned int nbModuli = poly_t::nmoduli;

  // Some needed polynomials
  poly_t &tmp = *alloc_aligned<poly_t, 32>(1),
	  &pka = *alloc_aligned<poly_t, 32>(1),
	  &pkaprime = *alloc_aligned<poly_t, 32>(1),
	  &pkb = *alloc_aligned<poly_t, 32>(1),
	  &pkbprime = *alloc_aligned<poly_t, 32>(1);
  
  // Gaussian sampler
  nfl::FastGaussianNoise<uint8_t, T, 2> g_prng(SIGMA, 128, 1<<10);
  
  // This step generates a secret key
  poly_t &s = *alloc_aligned<poly_t, 32>(1, nfl::gaussian<uint8_t, T, 2>(&g_prng));
  poly_t &sprime = *alloc_aligned<poly_t, 32>(1);
  s.ntt_pow_phi();
  sprime = nfl::compute_shoup(s);

  // This step generates a public key
  pka = nfl::uniform();
  pkb = nfl::gaussian<uint8_t, T, 2>(&g_prng, 2);
  pkb.ntt_pow_phi();

  // pkb = pkb + pka * s;
  pkb = pkb + nfl::shoup(pka * s, sprime);
  pkaprime = nfl::compute_shoup(pka);
  pkbprime = nfl::compute_shoup(pkb);
  

  start = std::chrono::steady_clock::now();
  // Generate REPETITIONS encryption of 0 formed of polynomials (resa[i],resb[i])
  for (unsigned j = 0; j < LOOP; j++)
  for (unsigned i = 0; i < REPETITIONS ; i++)
  {
    encrypt(resa[i], resb[i], pka, pkb, pkaprime, pkbprime, &g_prng);
  }
  end = std::chrono::steady_clock::now();
  std::cout << "Time per LWE-like symmetric encryption: " << get_time_us(start, end, LOOP*REPETITIONS) << " us" << std::endl;

  // Polynomial to store the sum of the plaintexts, should be 0 at the end
  start = std::chrono::steady_clock::now();
  // Decrypt benchmark
  for (unsigned j = 0; j < LOOP; j++)
  for (unsigned i = 0; i < REPETITIONS; i++)
  {
    decrypt(tmp, resa[i], resb[i], s, sprime, tmp.get_modulus(0));
  }
  end = std::chrono::steady_clock::now();
  std::cout << "Time per LWE-like symmetric decryption: " << get_time_us(start, end, LOOP*REPETITIONS) << " us" << std::endl;

  poly_t &finalres = *alloc_aligned<poly_t, 32>(1);
  // Now we decrypt the ciphertexts and add the results
  // As we have encryptions of 0 the sum should be 0
  for (unsigned i = 0; i < REPETITIONS; i++)
  {
    decrypt(tmp, resa[i], resb[i], s, sprime, tmp.get_modulus(0));
    // Add up the results of decryption
    finalres = finalres + tmp;
  }

  // Test correctness
  for (auto val: finalres)
  {
    if (val != 0)
    {
      std::cout << "ERROR finalres is " << val << std::endl;
      return 1;
    }
  }
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


