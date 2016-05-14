#ifndef  NFL_POLY_P_HPP
#define NFL_POLY_P_HPP

#include <memory>
#include "nfl/poly.hpp"
#include "nfl/aligned_allocator.hpp"

namespace nfl {

template<class T, size_t Degree, size_t NbModuli>
class poly_p {
  typedef poly<T, Degree, NbModuli> poly_type;
  typedef std::shared_ptr<poly_type> ptr_type;
 
public:
  using value_type = typename poly_type::value_type;
  using greater_value_type = typename poly_type::greater_value_type;
  static constexpr size_t nmoduli = poly_type::nmoduli;
  static constexpr size_t degree = poly_type::degree;
  static constexpr size_t nbits = poly_type::nbits;
  static constexpr size_t aggregated_modulus_bit_size = poly_type::aggregated_modulus_bit_size;

public:
  poly_p(poly_p const& o):
    _p(o._p)
  { }

  // AG: we need to declare this one or it will be catched by the generic one
  // below!
  poly_p(poly_p& o):
    _p(const_cast<poly_p const&>(o)._p)
  { }

  poly_p(poly_p&& o):
    _p(std::move(o._p))
  { }

  template <class... Args>
  poly_p(Args&& ... args):
    _p(make_pointer(std::forward<Args>(args)...))
  { }

  poly_p(poly_type const&) = delete;
  poly_p(poly_type&&) = delete;

public:
  poly_type& poly_obj()
  {
    detach();
    return *_p;
  }

  poly_type const& poly_obj() const { return *_p; }

public:
 template <class O>
 poly_p& operator=(O&& o) {
   poly_obj() = std::forward<O>(o);
   return *this;
 }

 poly_p& operator=(std::initializer_list<T> values) {
   poly_obj() = std::forward<std::initializer_list<T>>(values);
   return *this;
 }

  poly_p& operator=(poly_p const& o) 
  {
    if (this != &o) {
      _p = o._p;
    }
    return *this;
  }

  poly_p& operator=(poly_p&& o) 
  {
    if (this != &o) {
      _p = std::move(o._p);
    }
    return *this;
  }

  auto operator+(poly_p const& o) const -> decltype(poly_type() + poly_type())
  {
    return poly_obj() + o.poly_obj();
  }
  auto operator-(poly_p const& o) const -> decltype(poly_type() - poly_type())
  {
    return poly_obj() - o.poly_obj();
  }
  auto operator*(poly_p const& o) const -> decltype(poly_type() * poly_type())
  {
    return poly_obj() * o.poly_obj();
  }
  auto operator+(poly_type const& o) const -> decltype(poly_type() + poly_type())
  {
    return poly_obj() + o;
  }
  auto operator-(poly_type const& o) const -> decltype(poly_type() - poly_type())
  {
    return poly_obj() - o;
  }
  auto operator*(poly_type const& o) const -> decltype(poly_type() * poly_type())
  {
    return poly_obj() * o;
  }
  bool operator==(poly_p const& o) const
  {
    if (_p.get() == o._p.get()) {
      return true;
    }
    return *this == o.poly_obj();
  }
  bool operator!=(poly_p const& o) const
  {
    if (_p.get() == o._p.get()) {
      return false;
    }
    return *this != o.poly_obj();
  }

  template <class O>
  bool operator==(O const& o) const
  {
    return poly_obj() == o;
  }
  template <class O>
  bool operator!=(O const& o) const
  {
    return poly_obj() != o;
  }

  value_type& operator()(size_t cm, size_t i) { return poly_obj()(cm, i); }
  value_type const& operator()(size_t cm, size_t i) const { return poly_obj()(cm, i); }
  template<class M> auto load(size_t cm, size_t i) const -> decltype(M::load(&(this->operator()(cm, i)))) { return M::load(&(*this)(cm, i)); }

  static constexpr value_type get_modulus(size_t n) { return poly_type::get_modulus(n); }

  /* ntt stuff - public API
   */
  void ntt_pow_phi() { poly_obj().ntt_pow_phi();}
  void invntt_pow_invphi() { poly_obj().invntt_pow_invphi(); }

  /* manual serializers
  */
  void serialize_manually(std::ostream& outputstream) {
    poly_obj().serialize_manually(outputstream);
  }
  void deserialize_manually(std::istream& inputstream) {
    poly_obj().deserialize_manually(inputstream);
  }

  /* serializer (cereal)
  */
  template<class Archive> void serialize(Archive & archive) { 
    archive( poly_obj() );
  }

  /* set */
  void set(value_type v, bool reduce_coeffs = true) { poly_obj().set(v, reduce_coeffs); };
  void set(uniform const& mode) { poly_obj().set(mode); };
  void set(non_uniform const& mode) { poly_obj().set(mode); };
  template <class in_class, unsigned _lu_depth> void set(gaussian<in_class, T, _lu_depth> const& mode) { poly_obj().set(mode); };
  void set(std::initializer_list<value_type> values, bool reduce_coeffs = true) { poly_obj().set(values, reduce_coeffs); };
  void set(std::array<value_type, Degree> values, bool reduce_coeffs = true) { poly_obj().set(values, reduce_coeffs); };
  template <class It> void set(It first, It last, bool reduce_coeffs = true) { poly_obj().set(first, last, reduce_coeffs); };

private:
  template <class... Args>
  static ptr_type make_pointer(Args&& ... args)
  {
    aligned_allocator<poly_type, 32> alloc; 
    return std::allocate_shared<poly_type>(alloc, std::forward<Args>(args)...);
  }

  void detach()
  {
    // TODO: AG: I think there are issues w/ multithreading
    if (!_p.unique()) {
      // Copy the underlying object
      _p = make_pointer(*_p);
    }
  }

public:
  void set_mpz(mpz_t const& v) { poly_obj().set_mpz(v); };
  void set_mpz(std::array<mpz_t, Degree> const& values) { poly_obj().set_mpz(values); };
  void set_mpz(mpz_class const& v) { poly_obj().set_mpz(v); };
  void set_mpz(std::array<mpz_class, Degree> const& values) { poly_obj().set_mpz(values); };
  void set_mpz(std::initializer_list<mpz_class> const& values) { poly_obj().set_mpz(values); };
  template<class It> void set_mpz(It first, It last) { poly_obj().set_mpz(first, last); };

  inline std::array<mpz_t, Degree> poly2mpz() { return poly_obj().poly2mpz(); };
  inline void poly2mpz(std::array<mpz_t, Degree> & array) { poly_obj().poly2mpz(array); };
  inline void mpz2poly(std::array<mpz_t, Degree> const& array) { poly_obj().mpz2poly(array); };

  inline static constexpr size_t bits_in_moduli_product() { return poly_type::bits_in_moduli_product(); };
  inline static constexpr mpz_t& moduli_product() { return poly_type::moduli_product(); };
  inline static constexpr mpz_t& modulus_shoup() { return poly_type::modulus_shoup(); };
  inline static constexpr std::array<mpz_t, nmoduli> lifting_integers() { return poly_type::lifting_integers(); };

private:
  ptr_type _p;
};

/* misc type adaptor
 */
template<class T, size_t Degree, size_t AggregatedModulusBitSize>
using poly_p_from_modulus = poly_p<T, Degree, AggregatedModulusBitSize / params<T>::kModulusBitsize>;

/* stream operator
 */
template<class T, size_t Degree, size_t NbModuli>
std::ostream& operator<<(std::ostream& os, nfl::poly_p<T, Degree, NbModuli> const& p)
{
  return os << p.poly_obj();
}

/* unary operators
*/
template<class E, class T, size_t Degree, size_t NbModuli>
auto shoup(E const& e, poly_p<T, Degree, NbModuli> const& m) -> decltype(shoup(e, m.poly_obj()))
{
  return shoup(e, m.poly_obj());
}

template<class T, size_t Degree, size_t NbModuli>
auto compute_shoup(poly_p<T, Degree, NbModuli> const& p) -> decltype(compute_shoup(p.poly_obj()))
{
  return compute_shoup(p.poly_obj());
}

/* operator overloads - includes expression templates
 */
#define DECLARE_BINARY_OPERATOR_P(SYM, NAME)\
template<class T, size_t Degree, size_t NbModuli>\
auto SYM(nfl::poly_p<T, Degree, NbModuli> const& op0, nfl::poly_p<T, Degree, NbModuli> const& op1) -> decltype(ops::make_op<ops::NAME<T, CC_SIMD>>(op0, op1)) {\
  return ops::make_op<ops::NAME<T, CC_SIMD>>(op0, op1);\
}\
template<class T, size_t Degree, size_t NbModuli, class Op, class...Args>\
auto SYM(nfl::poly_p<T, Degree, NbModuli> const& op0, nfl::ops::expr<Op, Args...> const& op1) -> decltype(ops::make_op<ops::NAME<typename nfl::ops::expr<Op, Args...>::value_type, typename nfl::ops::expr<Op, Args...>::simd_mode>>(op0, op1)) {\
  return ops::make_op<ops::NAME<typename nfl::ops::expr<Op, Args...>::value_type, typename nfl::ops::expr<Op, Args...>::simd_mode>>(op0, op1);\
}\
template<class T, size_t Degree, size_t NbModuli, class Op, class...Args>\
auto SYM(nfl::ops::expr<Op, Args...> const& op0, nfl::poly_p<T, Degree, NbModuli> const& op1) -> decltype(ops::make_op<ops::NAME<typename nfl::ops::expr<Op, Args...>::value_type, typename nfl::ops::expr<Op, Args...>::simd_mode>>(op0, op1)){\
  return ops::make_op<ops::NAME<typename nfl::ops::expr<Op, Args...>::value_type, typename nfl::ops::expr<Op, Args...>::simd_mode>>(op0, op1);\
}

DECLARE_BINARY_OPERATOR_P(operator-, submod)
DECLARE_BINARY_OPERATOR_P(operator+, addmod)
DECLARE_BINARY_OPERATOR_P(operator*, mulmod)

} // nfl

#endif
