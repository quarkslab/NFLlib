#ifndef  NFL_POLY_P_HPP
#define NFL_POLY_P_HPP

#include <memory>
#include "nfl/poly.hpp"

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
  poly_p& operator=(O&& o)
  {
    poly_obj() = std::forward<T>(o);
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

  static constexpr value_type get_modulus(size_t n) { return poly_type::get_modulus(n); }

private:
  template <class... Args>
  static ptr_type make_pointer(Args&& ... args)
  {
    return std::make_shared<poly_type>(std::forward<Args>(args)...);
  }

  void detach()
  {
    // TODO: AG: I think there are issues w/ multithreading
    if (!_p.unique()) {
      // Copy the underlying object
      _p = make_pointer(*_p);
    }
  }

private:
  ptr_type _p;
};

template<class T, size_t Degree, size_t NbModuli>
std::ostream& operator<<(std::ostream& os, nfl::poly_p<T, Degree, NbModuli> const& p)
{
  return os << p.poly_obj();
}

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

} // nfl

#endif
