namespace nfl {

namespace ops {

// Fused Multiplications-Additions
// FMA with division for modular reduction (expensive)
template<class type, class tag> struct muladd;

template<class T>
struct muladd<T, simd::serial> {
  using simd_mode = simd::serial;
  T operator()(T rop, T x, T y, size_t cm) const
  {
    auto const p = params<T>::P[cm];
    ASSERT_STRICTMOD(x<p && y<p);
    using greater_value_type = typename params<T>::greater_value_type;
    greater_value_type res = (greater_value_type)(x) * y;
    ASSERT_STRICTMOD(res < (std::numeric_limits<greater_value_type>::max()-rop));
    res += rop;
    res %= p;

    ASSERT_STRICTMOD(res == ((((greater_value_type)(x) * y) % p) + rop) % p);
    return res;
  }
};

template<>
struct muladd<uint64_t, simd::serial> {
  using simd_mode = simd::serial;
  using T = uint64_t;
  T operator()(T rop, T x, T y, size_t cm) const
  {
    auto const p = params<T>::P[cm];
    auto const shift = params<T>::kModulusRepresentationBitsize;
    ASSERT_STRICTMOD((x<p) && (y<p));
    using greater_value_type = typename params<T>::greater_value_type;
    using value_type = typename params<T>::value_type;
    greater_value_type res = (greater_value_type)x * y;
    ASSERT_STRICTMOD(res < (std::numeric_limits<greater_value_type>::max()-rop));
    greater_value_type q = ((greater_value_type)params<T>::Pn[cm] * (res >> shift)) + (res<<2) ;
    value_type r  = res - (q>>shift) * p;
    if (r >= p) { r -= p; }
    r += rop;
    if (r >= p) { r -= p; }
    ASSERT_STRICTMOD(r == ((((greater_value_type)(x) * y) % p) + rop) % p);
    return r;
  }
};

// FMA using Shoup's precomputed approach again
// Works only if yshoup = ((uint128_t) y << 64) / p (for T=uint64_t)
// has been pre-computed AND p is small enough (2 bits less than limb).
// OUTPUT: x * y mod p
template<class T, class tag> struct muladd_shoup;

template<class T>
struct muladd_shoup<T, simd::serial> {
  using simd_mode = simd::serial;
T operator()(T rop, T x, T y, T yprime, size_t cm) const
{
  auto const p = params<T>::P[cm];
  using greater_value_type = typename params<T>::greater_value_type;
  ASSERT_STRICTMOD(rop<p && x<p && y<p);
  T q =
    ((greater_value_type) x * yprime) >> params<T>::kModulusRepresentationBitsize;
#ifdef ENFORCE_STRICTMOD
  using value_type = typename params<T>::value_type;
  value_type res;
  res = x * y - q * p;
  res -= ((res>=p) ? p : 0);
  rop += res;
#else
  rop += x * y - q * p;
#endif
  ASSERT_STRICTMOD(rop - ((rop>=p) ? p : 0)<p);
  return rop - ((rop>=p) ? p : 0);
}
};



} // ops

} // nfl
