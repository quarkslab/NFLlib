#include <cereal/archives/binary.hpp>
#include <chrono>
#include <iostream>
#include <iterator>
#include <nfl.hpp>
#include <sstream>
#include "tools.h"

#define REPETITIONS 500

template <class T>
double get_time_us(T const& start, T const& end, uint32_t N) {
  auto diff = end - start;
  return (long double)(std::chrono::duration_cast<std::chrono::microseconds>(
                           diff)
                           .count()) /
         N;
}

/// Temporary class to compare manual serialization and serialization using
/// cereal
template <class T, size_t degree, size_t modulus>
class poly_temp {
  using poly_p = nfl::poly_p_from_modulus<T, degree, modulus>;
  static constexpr size_t N = poly_p::degree * poly_p::nmoduli;

 public:
  T _data[N] __attribute__((aligned(32))) = {0};
  poly_temp() {}
  void set(poly_p const& p) {
    for (size_t cm = 0; cm < poly_p::nmoduli; cm++) {
      for (size_t i = 0; i < poly_p::degree; i++) {
        _data[cm * poly_p::degree + i] = p(cm, i);
      }
    }
  }
  T& operator()(size_t cm, size_t i) { return _data[cm * degree + i]; }

  /// cereal serializer
  template <class Archive>
  void serialize(Archive& archive) {
    archive(_data);  // serialize coefficients by passing them to the archive
  }

  /// manual serializers
  inline void oserialize_manually(std::ostream& outputstream) {
    outputstream.write(reinterpret_cast<char*>(_data), N * sizeof(T));
  }
  inline void iserialize_manually(std::istream& inputstream) {
    inputstream.read(reinterpret_cast<char*>(_data), N * sizeof(T));
  }
} __attribute__((aligned(32)));

template <size_t degree, size_t modulus, class T>
bool run() {
  using poly_t = nfl::poly_from_modulus<T, degree, modulus>;

  // Stream
  std::stringstream ss(std::stringstream::in | std::stringstream::out |
                       std::stringstream::binary);

  // Archives
  cereal::BinaryOutputArchive oarchive(ss);  // output archive
  cereal::BinaryInputArchive iarchive(ss);   // input archive

  // define a random polynomial
  poly_t& p0 = *alloc_aligned<poly_t, 32>(1, nfl::uniform());

  // serialize the polynomial
  oarchive(p0);

  // copy p0 into p1
  poly_t& p1 = *alloc_aligned<poly_t, 32>(1, p0);

  // free p0
  free_aligned(1, &p0);

  // define a new polynomial p2
  poly_t& p2 = *alloc_aligned<poly_t, 32>(1);

  // deserialize the stream into p2
  iarchive(p2);

  // Verify that p1 == p2
  bool ret_value = (p1 == p2);

  // Cleaning
  free_aligned(1, &p1);
  free_aligned(1, &p2);

  return ret_value;
}

template <size_t degree, size_t modulus, class T>
bool run_p() {
  using poly_p = nfl::poly_p_from_modulus<T, degree, modulus>;

  // Stream
  std::stringstream ss(std::stringstream::in | std::stringstream::out |
                       std::stringstream::binary);

  // Archives
  cereal::BinaryOutputArchive oarchive(ss);  // output archive
  cereal::BinaryInputArchive iarchive(ss);   // input archive

  // define a random polynomial
  poly_p* p0 = new poly_p(nfl::uniform());

  // serialize the polynomial
  oarchive(*p0);

  // copy p0 into p1
  poly_p p1 = *p0;

  // delete p0
  delete p0;

  // define a new polynomial p2
  poly_p p2;

  // deserialize the stream into p2
  iarchive(p2);

  // Verify that p1 == p2
  return (p1 == p2);
}

template <size_t degree, size_t modulus, class T>
bool run_perf() {
  using poly_p = nfl::poly_p_from_modulus<T, degree, modulus>;
  using poly_tmp = poly_temp<T, degree, modulus>;

  auto start = std::chrono::steady_clock::now();
  auto end = std::chrono::steady_clock::now();

  // Stream
  std::stringstream ss(std::stringstream::in | std::stringstream::out |
                       std::stringstream::binary);

  // Archives
  cereal::BinaryOutputArchive oarchive(ss);  // output archive
  cereal::BinaryInputArchive iarchive(ss);   // input archive

  // array of random polynomials
  std::array<poly_p, REPETITIONS> poly_p_array;
  std::fill(poly_p_array.begin(), poly_p_array.end(), nfl::uniform());

  // array of poly_tmp
  poly_tmp* poly_tmp_array = alloc_aligned<poly_tmp, 32>(REPETITIONS);
  for (size_t i = 0; i < REPETITIONS; i++) {
    poly_tmp_array[i].set(poly_p_array[i]);
  }

  // Benchmark cereal serialization
  start = std::chrono::steady_clock::now();
  for (size_t i = 0; i < REPETITIONS; i++) {
    oarchive(poly_tmp_array[i]);
  }
  end = std::chrono::steady_clock::now();
  std::cout << "Cereal output serialization: "
            << get_time_us(start, end, REPETITIONS) << " us" << std::endl;

  start = std::chrono::steady_clock::now();
  for (size_t i = 0; i < REPETITIONS; i++) {
    iarchive(poly_tmp_array[i]);
  }
  end = std::chrono::steady_clock::now();
  std::cout << "Cereal input serialization: "
            << get_time_us(start, end, REPETITIONS) << " us" << std::endl;

  // reset stringstream
  ss.clear();
  ss.str("");

  // Benchmark cereal serialization
  start = std::chrono::steady_clock::now();
  for (size_t i = 0; i < REPETITIONS; i++) {
    poly_tmp_array[i].oserialize_manually(ss);
  }
  end = std::chrono::steady_clock::now();
  std::cout << "Manual output serialization: "
            << get_time_us(start, end, REPETITIONS) << " us" << std::endl;

  start = std::chrono::steady_clock::now();
  for (size_t i = 0; i < REPETITIONS; i++) {
    poly_tmp_array[i].iserialize_manually(ss);
  }
  end = std::chrono::steady_clock::now();
  std::cout << "Manual input serialization: "
            << get_time_us(start, end, REPETITIONS) << " us" << std::endl;

  // Verify serialization
  bool ret_value = true;
  for (size_t i = 0; i < REPETITIONS; i++) {
    // reset stringstream
    ss.clear();
    ss.str("");

    // poly_p serialized to ss2, and read as a poly_tmp
    oarchive(poly_p_array[i]);
    poly_tmp* p_tmp = alloc_aligned<poly_tmp, 32>(1);
    p_tmp->iserialize_manually(ss);
    for (size_t cm = 0; cm < poly_p::nmoduli; cm++) {
      for (size_t j = 0; j < poly_p::degree; j++) {
        ret_value &= ((*p_tmp)(cm, j) == poly_p_array[i](cm, j));
      }
    }

    // reset stringstream
    ss.clear();
    ss.str("");

    // poly_tmp serialized to ss3, and read as a poly_p
    p_tmp->oserialize_manually(ss);
    poly_p p;
    iarchive(p);
    for (size_t cm = 0; cm < poly_p::nmoduli; cm++) {
      for (size_t j = 0; j < poly_p::degree; j++) {
        ret_value &= ((*p_tmp)(cm, j) == p(cm, j));
      }
    }

    free_aligned(1, p_tmp);
  }

  free_aligned(REPETITIONS, poly_tmp_array);

  return ret_value;
}

int main(int argc, char const* argv[]) {
  return not(run<CONFIG>() and run_p<CONFIG>() and run_perf<CONFIG>());
}
