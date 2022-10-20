#ifndef NFL_CORE_HPP
#define NFL_CORE_HPP

#include <type_traits>
#include <vector>
#include <numeric>
#include <algorithm>

#include "nfl/poly.hpp"
#include "nfl/ops.hpp"
#include "nfl/algos.hpp"
#include "nfl/permut.hpp"
#include <type_traits>

namespace nfl {

/* poly methods implementation
 */
template<class T, size_t Degree, size_t NbModuli>
template<class Op, class... Args> inline poly<T, Degree, NbModuli>::poly(ops::expr<Op, Args...> const& expr) {
  *this = expr;
}

template<class T, size_t Degree, size_t NbModuli>
template<class Op, class... Args> inline poly<T, Degree, NbModuli>& poly<T, Degree, NbModuli>::operator=(ops::expr<Op, Args...> const& expr) {
  using E = ops::expr<Op, Args...>;
  constexpr size_t vector_size = E::simd_mode::template elt_count<T>::value;
  constexpr size_t vector_bound = degree / vector_size * vector_size;
  static_assert(vector_bound == degree, "no need for a footer");

  for(size_t cm = 0; cm < nmoduli; ++cm) {
    for(size_t j = 0; j < vector_bound; j+= vector_size) {
      E::simd_mode::store(&(*this)(cm,j), expr.template load<typename E::simd_mode>(cm, j));
    }
  }
  return *this;
}

template<class T, size_t Degree, size_t NbModuli>
inline poly<T, Degree, NbModuli>::operator bool() const {
  return std::find_if(begin(), end(),
                      [](value_type v) { return v != 0; }) != end();
}

template<class T, size_t Degree, size_t NbModuli>
typename poly<T, Degree, NbModuli>::core poly<T, Degree, NbModuli>::base;

// *********************************************************
//    Constructor
// *********************************************************

template<class T, size_t Degree, size_t NbModuli> poly<T, Degree, NbModuli>::core::core()
{

  static_assert(aggregated_modulus_bit_size % params<T>::kModulusBitsize == 0,
      "core: CRITICAL, aggregated_modulus_bit_size is not a multiple of kModulusBitsize (see params.hpp)");
  static_assert(aggregated_modulus_bit_size / params<T>::kModulusBitsize <= params<T>::kMaxNbModuli,
      "core: CRITICAL. Not enough moduli of this size to reach the requested aggregated_modulus_bit_size (see params.hpp)");
  static_assert(degree <= params<T>::kMaxPolyDegree,
      "core: CRITICAL, degree is not lower or equal than kMaxPolyDegree (see params.hpp)");
  initialize();
}

template<class T, size_t Degree, size_t NbModuli>
poly<T, Degree, NbModuli>::poly() : poly(std::integral_constant<T, 0>::value) {
}

template<class T, size_t Degree, size_t NbModuli>
poly<T, Degree, NbModuli>::poly(value_type v, bool reduce_coeffs) {
  set(v, reduce_coeffs);
}

template<class T, size_t Degree, size_t NbModuli>
poly<T, Degree, NbModuli>::poly(std::initializer_list<value_type> values, bool reduce_coeffs) {
  set(values, reduce_coeffs);
}

template<class T, size_t Degree, size_t NbModuli>
template<class It>
poly<T, Degree, NbModuli>::poly(It first, It last, bool reduce_coeffs) {
  set(first, last, reduce_coeffs);
}

template<class T, size_t Degree, size_t NbModuli>
void poly<T, Degree, NbModuli>::set(value_type v, bool reduce_coeffs) {
  if (v == 0) {
    // CRITICAL: the object must be 32-bytes aligned to avoid vectorization issues
    assert((unsigned long)(this->_data) % 32 == 0);
    std::fill(begin(), end(), 0);
  }
  else {
    set({v}, reduce_coeffs);
  }
}

template<class T, size_t Degree, size_t NbModuli>
void poly<T, Degree, NbModuli>::set(std::initializer_list<value_type> values, bool reduce_coeffs) {
  set(values.begin(), values.end(), reduce_coeffs);
}

template <class T, size_t Degree, size_t NbModuli>
template <class It>
void poly<T, Degree, NbModuli>::set(It first, It last, bool reduce_coeffs) {
  // CRITICAL: the object must be 32-bytes aligned to avoid vectorization issues
  assert((unsigned long)(this->_data) % 32 == 0);

  // The initializer needs to have either less values than the polynomial degree
  // (and the remaining coefficients are set to 0), or be fully defined (i.e.
  // the degree*nmoduli coefficients needs to be provided)
  size_t size = std::distance(first, last);
  if (size > degree && size != degree * nmoduli) {
    throw std::runtime_error(
        "core: CRITICAL, initializer of size above degree but not equal "
        "to nmoduli*degree");
  }

  auto* iter = begin();
  auto viter = first;

  for (size_t cm = 0; cm < nmoduli; cm++) {
    auto const p = get_modulus(cm);

    if (size != degree * nmoduli) viter = first;

    // set the coefficients
    size_t i = 0;
    for (; i < degree && viter < last; ++i, ++viter, ++iter) {
      *iter = reduce_coeffs ? (*viter) % p : *viter;
    }

    // pad with zeroes if needed
    for (; i < degree; ++i, ++iter) {
      *iter = 0;
    }
  }
}

// ****************************************************
// Random polynomial generation functions
// ****************************************************

// Sets a pre-allocated random polynomial in FFT form
// uniformly random, else the coefficients are uniform below the bound

template<class T, size_t Degree, size_t NbModuli>
poly<T, Degree, NbModuli>::poly(uniform const& u) {
  set(u);
}

template<class T, size_t Degree, size_t NbModuli>
void poly<T, Degree, NbModuli>::set(uniform const &) {
  // CRITICAL: the object must be 32-bytes aligned to avoid vectorization issues
  assert((unsigned long)(this->_data) % 32 == 0);

  // In uniform mode we need randomness for all the polynomials in the CRT
  fastrandombytes((unsigned char *)data(), sizeof(poly));

  for (unsigned int cm = 0; cm < nmoduli; cm++) {
    // In the uniform case, instead of getting a big random (within the general
    // moduli), We rather prefer, for performance issues, to get smaller
    // randoms for each module The mask should be the same for all moduli
    // (because they are the same size) But for generality we prefer to compute
    // it for each moduli so that we could have moduli of different bitsize

    value_type mask =
        (1ULL << (int)(floor(log2(get_modulus(cm))) + 1)) - 1;

    for (size_t i = 0; i < degree; i++) {
      // First remove the heavy weight bits we dont need
      value_type tmp = _data[i + degree * cm] & mask;

      // When the random is still too large, reduce it
      if (tmp >= get_modulus(cm)) {
        tmp -= get_modulus(cm);
      }
      _data[i + degree * cm] = tmp;
    }
  }
#ifdef CHECK_STRICTMOD
  for (size_t cm = 0; cm < nmoduli; cm++) {
    for (size_t i = 0; i < degree; i++) {
      assert(_data[i + degree * cm] < get_modulus(cm));
    }
  }
#endif
}

template<class T, size_t Degree, size_t NbModuli>
poly<T, Degree, NbModuli>::poly(non_uniform const& mode) {
  set(mode);
}

template<class T, size_t Degree, size_t NbModuli>
void poly<T, Degree, NbModuli>::set(non_uniform const& mode) {
  // CRITICAL: the object must be 32-bytes aligned to avoid vectorization issues
  assert((unsigned long)(this->_data) % 32 == 0);

  uint64_t const upper_bound = mode.upper_bound;
  uint64_t const amplifier = mode.amplifier;
  // In bounded mode upper_bound must be below the smaller of the moduli
  for (unsigned int cm = 0; cm < nmoduli; cm++) {
    if (upper_bound >= get_modulus(cm)) {
      throw std::runtime_error(
          "core: upper_bound is larger than the modulus");
    }
  }

  // We play with the rnd pointer (in the uniform case), and thus
  // we need to remember the allocated pointer to free it at the end
  value_type rnd[degree];

  // Get some randomness from the PRNG
  fastrandombytes((unsigned char *)rnd, sizeof(rnd));

  // upper_bound is below the moduli so we create the same mask for all the
  // moduli
  value_type mask =
      (1ULL << (int)(floor(log2(2*upper_bound-1)) + 1)) - 1;

  if(amplifier == 1){
    for (unsigned int i = 0; i < degree; i++) {

      // First remove the heavy weight bits we dont need
      value_type tmp = rnd[i] & mask;

      // When the random is still too large, reduce it
      // In order to follow strictly a uniform distribution we should
      // get another rnd but in order to follow the proofs of security
      // strictly we should also take noise from a gaussian ...
      if (tmp >= 2*upper_bound-1) {
        tmp -= 2*upper_bound-1;
      }
      // Center the noise
      if (tmp >= upper_bound) {
        for (unsigned int cm = 0; cm < nmoduli; cm++) {
          _data[degree * cm + i] = get_modulus(cm) + tmp - (2*upper_bound-1);
        }
      }
      else
      {
        for (unsigned int cm = 0; cm < nmoduli; cm++) {
          _data[degree * cm + i] = tmp;
        }
      }
    }
  }
  else
  {
    for (unsigned int i = 0; i < degree; i++) {

      // First remove the heavy weight bits we dont need
      value_type tmp = rnd[i] & mask;

      // When the random is still too large, reduce it
      // In order to follow strictly a uniform distribution we should
      // get another rnd but in order to follow the proofs of security
      // strictly we should also take noise from a gaussian ...
      if (tmp >= 2*upper_bound-1) {
        tmp -= 2*upper_bound-1;
      }
      // Center the noise
      if (tmp >= upper_bound) {
        for (unsigned int cm = 0; cm < nmoduli; cm++) {
          _data[degree * cm + i] = get_modulus(cm) + tmp*amplifier - (2*upper_bound-1)*amplifier;
        }
      }
      else
      {
        for (unsigned int cm = 0; cm < nmoduli; cm++) {
          _data[degree * cm + i] = tmp*amplifier;
        }
      }
    }
  }
#ifdef CHECK_STRICTMOD
  for (size_t cm = 0; cm < nmoduli; cm++) {
    for (size_t i = 0; i < degree; i++) {
      assert(_data[i + degree * cm] < get_modulus(cm));
    }
  }
#endif
}

template<class T, size_t Degree, size_t NbModuli>
template<class in_class, unsigned _lu_depth>
poly<T, Degree, NbModuli>::poly(gaussian<in_class, T, _lu_depth> const& mode) {
  set(mode);
}

template<class T, size_t Degree, size_t NbModuli>
template<class in_class, unsigned _lu_depth>
void poly<T, Degree, NbModuli>::set(gaussian<in_class, T, _lu_depth> const& mode) {
  // CRITICAL: the object must be 32-bytes aligned to avoid vectorization issues
  assert((unsigned long)(this->_data) % 32 == 0);

  uint64_t const amplifier = mode.amplifier;

  // We play with the rnd pointer (in the uniform case), and thus
  // we need to remember the allocated pointer to free it at the end
  signed_value_type rnd[degree];

  // Get some randomness from the PRNG
  mode.fg_prng->getNoise((value_type *)rnd, degree);

  if (amplifier != 1) for (unsigned int i = 0; i < degree; i++) rnd[i]*= amplifier;
  for (size_t cm = 0; cm < nmoduli; cm++)
  {
    for (size_t i = 0 ; i < degree; i++)
    {
      if(rnd[i]<0)
        _data[degree*cm+i] = get_modulus(cm) + rnd[i];
      else
        _data[degree*cm+i] = rnd[i];
    }
    //memcpy(_data+degree*cm, rnd, degree*sizeof(value_type));
  }


#ifdef CHECK_STRICTMOD
  for (size_t cm = 0; cm < nmoduli; cm++) {
    for (size_t i = 0; i < degree; i++) {
      assert(_data[i + degree * cm] < get_modulus(cm));
    }
  }
#endif
}

template<class T, size_t Degree, size_t NbModuli>
poly<T, Degree, NbModuli>::poly(ZO_dist const& mode) {
  set(mode);
}

template<class T, size_t Degree, size_t NbModuli>
void poly<T, Degree, NbModuli>::set(ZO_dist const& mode) {
  uint8_t rnd[Degree];
  fastrandombytes(rnd, sizeof(rnd));
  value_type *ptr = &_data[0];
  for (size_t cm = 0; cm < NbModuli; ++cm) {
    const T pm = params<T>::P[cm] - 1u;
    /* sample {-1, 0, 1} */
    for (size_t i = 0; i < Degree; ++i, ptr++) {
      *ptr = rnd[i] <= mode.rho ? (pm + (rnd[i] & 2)) : 0u;
      if (*ptr > params<T>::P[cm]) {
        *ptr -= params<T>::P[cm];
      }
    }
  }
}

template<class T, size_t Degree, size_t NbModuli>
poly<T, Degree, NbModuli>::poly(hwt_dist const& mode) {
  set(mode);
}

template<class T, size_t Degree, size_t NbModuli>
void poly<T, Degree, NbModuli>::set(hwt_dist const& mode) {
  assert(mode.hwt > 0 && mode.hwt <= Degree);
  std::vector<size_t> hitted(mode.hwt);
  std::iota(hitted.begin(), hitted.end(), 0U); // select the first hwt positions.
  std::vector<size_t> rnd(hitted.size());
  auto rnd_end = rnd.end();
  auto rnd_ptr = rnd_end;
  /* Reservoir Sampling: uniformly select hwt coefficients. */
  for (size_t k = mode.hwt; k < degree; ++k)
  {
    size_t pos = 0;
    size_t reject_sample = std::numeric_limits<size_t>::max() / k;
    /* sample uniformly from [0, k) using reject sampling. */
    for (;;) {
      if (rnd_ptr == rnd_end)
      {
        fastrandombytes((unsigned char *)rnd.data(), rnd.size() * sizeof(size_t));
        rnd_ptr = rnd.begin();
      }
      pos = *rnd_ptr++;
      if (pos <= reject_sample * k) {
        pos %= k;
        break;
      }
    }
    if (pos < mode.hwt)
      hitted[pos] = k;
  }

  std::sort(hitted.begin(), hitted.end()); // for better locality ?
  std::memset(_data, 0x0, N * sizeof(value_type)); // clear up all
  fastrandombytes((unsigned char *)rnd.data(), rnd.size() * sizeof(size_t));
  for (size_t cm = 0, offset = 0; cm < NbModuli; ++cm, offset += degree) {
    const T pm = params<T>::P[cm] - 1u;
    rnd_ptr = rnd.begin();
    for (size_t pos : hitted) {
        _data[pos + offset] = ((*rnd_ptr++) & 2U); // {-1, 1}
        _data[pos + offset] = (_data[pos + offset] > 0 ? 1 : pm);
    }
  }
  std::memset(hitted.data(), 0x0, hitted.size() * sizeof(size_t)); // erase from memory
}

// *********************************************************
// Helper functions
// *********************************************************


template<class T, size_t Degree, size_t NbModuli>
std::ostream& operator<<(std::ostream& outs, poly<T, Degree, NbModuli> const& p)
{
  bool first = true;
  std::string term;
  if (typeid(T) == typeid(uint64_t)) term = "ULL";
  else if (typeid(T) == typeid(uint32_t)) term = "UL";
  else term = "U";

  outs << "{ ";
  for(auto v : p)
  {
    if (first)
    {
      first = false;
      outs << v ;
    }
    else
    {
      outs << term << ", " << v;
    }
  }
  return outs << term << " }";
}


// *********************************************************
//     Operations over polynomials
// *********************************************************

// Apply mulmodShoup to all the coefficients of the polynomials (much faster)
// This is a polynomial multiplication mod X**degree+1 if polynomials
// in a coefficient representation have been processed through nttAndPowPhi
// op2prime must have been obtained from op2 with allocandcomputeShouppoly
//

// *********************************************************
//     Number-Theoretic Functions
// *********************************************************

// *****************************************************************
// NTT functions from David Harvey
// From : http://web.maths.unsw.edu.au/~davidharvey/papers/fastntt/
// Most of the functions have been modified for our needs
// *****************************************************************

// Number-Theoretic Transform: replaces a coefficient representation of a
// polynomial by the values of the polynomial on a set of points. This
// allows to do coefficient wise multiplication of polynomials instead
// of the usual convolution product
// - x points to the polynomial to which we will apply the NTT
// - wtab is a pre-computed table with the powers of a root of unity
// - winvtab is another pre-computed table with the powers of the inverse
// of a root of unity
// - K is the log_2 of the degree of the polynomial
// - p is the modulus coefficient operations are done with
// The NTT is computed in-place on x, so x is the output
template<class T, size_t Degree, size_t NbModuli> bool poly<T, Degree, NbModuli>::core::ntt(value_type* x, const value_type* wtab, const value_type* winvtab, value_type const p)
{
#ifdef CHECK_STRICTMOD
  for (size_t i = 0 ; i < degree ; i++)
  {
    ASSERT_STRICTMOD(x[i] < p);
  }
#endif

#ifdef NTT_STRICTMOD
  value_type* x_orig = x;
#endif

  if (degree == 1)
    return true;

  // special case
  if (degree == 2)
  {
    value_type u0 = x[0];
    value_type u1 = x[1];
    value_type t0 = u0 + u1;
    value_type t1 = u0 - u1;
    t0 -= (t0 >= 2*p) ? (2*p) : 0;
    t1 += ((typename std::make_signed<value_type>::type) t1 < 0) ? (2*p) : 0;
    x[0] = t0;
    x[1] = t1;
    return true;
  }
  const size_t M = ops::ntt_loop<CC_SIMD, poly>::run(x, wtab, winvtab, p);

  typedef typename std::make_signed<value_type>::type signed_value_type;
  // last two layers
  for (size_t r = 0; r < M; r++, x += 4)
  {
    value_type u0 = x[0];
    value_type u1 = x[1];
    value_type u2 = x[2];
    value_type u3 = x[3];

    value_type v0 = u0 + u2;
    v0 -= (v0 >= 2*p) ? (2*p) : 0;
    value_type v2 = u0 - u2;
    v2 += ((signed_value_type) v2 < 0) ? (2*p) : 0;

    value_type v1 = u1 + u3;
    v1 -= (v1 >= 2*p) ? (2*p) : 0;
    value_type t = u1 - u3 + 2*p;

    value_type q = ((greater_value_type) t * winvtab[1]) >> params<T>::kModulusRepresentationBitsize;
    value_type v3 = t * wtab[1] - q * p;

    value_type z0 = v0 + v1;
    z0 -= (z0 >= 2*p) ? (2*p) : 0;
    value_type z1 = v0 - v1;
    z1 += ((signed_value_type) z1 < 0) ? (2*p) : 0;

    value_type z2 = v2 + v3;
    z2 -= (z2 >= 2*p) ? (2*p) : 0;
    value_type z3 = v2 - v3;
    z3 += ((signed_value_type) z3 < 0) ? (2*p) : 0;

    x[0] = z0;
    x[1] = z1;
    x[2] = z2;
    x[3] = z3;
  }

#ifdef NTT_STRICTMOD
  for (size_t i = 0; i < degree; i++)
  {
    x_orig[i]-= ((x_orig[i]>=p)? p : 0);
    ASSERT_STRICTMOD(x_orig[i] < p);
  }
#endif

  return true;
}


// Inverse NTT: replaces NTT values representation by the classic
// coefficient representation, parameters are the same as for ntt() with
// inv_wtab, inv_winvtab tables for the inverse operation and invK the
// inverse of the polynomialDegree
template<class T, size_t Degree, size_t NbModuli> inline bool poly<T, Degree, NbModuli>::core::inv_ntt(value_type * x, const value_type* const inv_wtab, const value_type* const inv_winvtab,
    const value_type invK, value_type const p)
{
  alignas(32) value_type y[degree+1];

  if (degree == 1)
    return true;

  // bit-reverse
  //r_set<0, 1, degree>{}(y, x);
  permut<degree>::compute(y, x);

  ntt(y, inv_wtab, inv_winvtab, p);

  // bit-reverse again
  permut<degree>::compute(x, y);
  //r_set<0, 1, degree>{}(x, y);
  return true;
}


// This function is private but it is kept here to group all the functions
// by David Harvey
// This function just computes the powers of the root of unity and its inverse
// It is used by initialize() and results are used by ntt and inv_ntt
template<class T, size_t Degree, size_t NbModuli> inline void poly<T, Degree, NbModuli>::core::prep_wtab(value_type* wtab, value_type* wtabshoup, value_type w, size_t cm)
{
  auto const p = get_modulus(cm);
  unsigned K = degree;

  while (K >= 2)
  {
    greater_value_type wi = 1;     // w^i
    for (size_t i = 0; i < K/2; i++)
    {
      *wtab++ = wi;
      *wtabshoup++ = ((greater_value_type) wi << params<T>::kModulusRepresentationBitsize) / p;
      wi = ops::mulmod<T, simd::serial>{}(static_cast<T>(wi), w, cm);
    }
    w = ops::mulmod<T, simd::serial>{}(w, w, cm);
    K /= 2;
  }
}
// *****************************************************************
// END OF NTT functions from David Harvey
// From : http://web.maths.unsw.edu.au/~davidharvey/papers/fastntt/
// Most of the functions have been modified for our needs
// *****************************************************************

// Number-Theoretic Transform: replaces a coefficient representation of a
// polynomial by the values of the polynomial on a set of points. This
// allows to do coefficient wise multiplication of polynomials instead
// of the usual convolution product. In order to have polynomial
// operations mod X**n + 1 we must multiply the polynomial coefficients
// by phi powers before doing the NTT
template<class T, size_t Degree, size_t NbModuli> inline void poly<T, Degree, NbModuli>::core::ntt_pow_phi(poly& op)
{
  op = nfl::shoup(op * reinterpret_cast<poly const &>(phis), reinterpret_cast<poly const &>(shoupphis));
  for(size_t currentModulus = 0; currentModulus < nmoduli; ++currentModulus) {
    poly::core::ntt(&op(currentModulus, 0), omegas[currentModulus], shoupomegas[currentModulus], get_modulus(currentModulus));
  }
}


// Inverse NTT: replaces NTT values representation by the classic
// coefficient representation and multiplies the coefficients with the
// inverse powers of phi
// In order to have polynomial operations mod X**n + 1 as we want
// we must multiply the polynomial coefficients by invphi powers before doing the inverse-NTT
template<class T, size_t Degree, size_t NbModuli> inline void poly<T, Degree, NbModuli>::core::invntt_pow_invphi(poly &op)
{
  for(size_t currentModulus = 0; currentModulus < nmoduli; ++currentModulus) {
    poly::core::inv_ntt(&op(currentModulus, 0), invomegas[currentModulus], shoupinvomegas[currentModulus], invpolyDegree[currentModulus], get_modulus(currentModulus));
  }
  op = nfl::shoup(op * reinterpret_cast<poly const &>(invpoly_times_invphis), reinterpret_cast<poly const &>(shoupinvpoly_times_invphis));
}

// *********************************************************
//    Initializing function
// *********************************************************

// Initializing function, called from the constructor. Computes all
// data needed to do NFL (NTT, inverse-NTT, lifting) operations
// Depends on the template parameters.
// TODO: All that is done in this function can be generated at compile
// time through meta-programming and statically defined.
template<class T, size_t Degree, size_t NbModuli> void  poly<T, Degree, NbModuli>::core::initialize()
{
  for(size_t currentModulus = 0; currentModulus < nmoduli; ++currentModulus)
  {
    value_type phi, invphi, omega, invomega, temp;
    shoupomegas[currentModulus] = omegas[currentModulus] + degree;
    shoupinvomegas[currentModulus] = invomegas[currentModulus] + degree;

    // We start by computing phi
    // The roots in the array are primitve
    // X=2*params<T>::kMaxPolyDegree-th roots
    // Squared log2(X)-log2(degree) times they become
    // degree-th roots as required by the NTT
    // But first we get phi = sqrt(omega) squaring them
    // log2(X/2)-log2(degree) times
    phi = params<T>::primitive_roots[currentModulus];
    for (unsigned int i = 0 ;
        i < static_log2<params<T>::kMaxPolyDegree>::value - static_log2<degree>::value ; i++)
    {
      phi = ops::mulmod<T, simd::serial>{}(phi, phi, currentModulus);
    }

    // Now that temp = phi we initialize the array of phi**i values
    // Initialized to phi**0
    temp = 1;
    for (unsigned int i = 0 ; i < degree ; i++)
    {
      phis[currentModulus][i] = temp;
      shoupphis[currentModulus][i] = ((greater_value_type) temp << params<T>::kModulusRepresentationBitsize) / get_modulus(currentModulus);
      // phi**(i+1)
      temp = ops::mulmod<T, simd::serial>{}(temp, phi, currentModulus);
    }
    // At the end of the loop temp = phi**degree

    // Computation of invphi
    // phi**(2*polydegree)=1 -> temp*phi**(degree-1) = phi**(-1)
    invphi = ops::mulmod<T, simd::serial>{}(temp, phis[currentModulus][degree-1], currentModulus);

    // Computation of the inverse of degree using the inverse of kMaxPolyDegree
    invpolyDegree[currentModulus] = ops::mulmod<T, simd::serial>{}(params<T>::invkMaxPolyDegree[currentModulus],
        static_cast<T>(params<T>::kMaxPolyDegree/degree), currentModulus);

    // Now we can compute the table invpoly_times_invphis
    temp = invpolyDegree[currentModulus];
    for (unsigned int i = 0 ; i < degree ; i++)
    {
      invpoly_times_invphis[currentModulus][i] = temp;
      shoupinvpoly_times_invphis[currentModulus][i] = ((greater_value_type) temp << params<T>::kModulusRepresentationBitsize)
        / get_modulus(currentModulus);
      // This is invpolyDegree*invphi**(i+1)
      temp = ops::mulmod<T, simd::serial>{}(temp, invphi, currentModulus);
    }

    // For the omegas it is easy, we just use the function of David Harvey modified for our needs
    omega = ops::mulmod<T, simd::serial>{}(phi, phi, currentModulus);
    prep_wtab(omegas[currentModulus], shoupomegas[currentModulus], omega, currentModulus);

    // And again for the invomegas
    invomega = ops::mulmod<T, simd::serial>{}(invphi, invphi, currentModulus);
    prep_wtab(invomegas[currentModulus], shoupinvomegas[currentModulus], invomega, currentModulus);
  }
}

}

#endif
