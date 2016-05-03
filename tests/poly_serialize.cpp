#include <sstream>
#include <cereal/archives/binary.hpp>
#include <nfl.hpp>
#include "tools.h"


template <size_t degree, size_t modulus, class T>
bool run() {
  using poly_t = nfl::poly_from_modulus<T, degree, modulus>;

  // Stream
  std::stringstream ss;

  // Archives
  cereal::BinaryOutputArchive oarchive(ss); // output archive
  cereal::BinaryInputArchive iarchive(ss); // input archive

  // define a random polynomial
  poly_t& p0 = *alloc_aligned<poly_t, 32>(1, nfl::uniform());

  // serialize the polynomial
  oarchive(p0);

  // define a new polynomial
  poly_t& p1 = *alloc_aligned<poly_t, 32>(1);

  // deserialize the stream
  iarchive(p1);

  // Verify that p0 == p1
  bool ret_value = (p0 == p1);
  
  // Cleaning
  free_aligned(1, &p0);
  free_aligned(1, &p1);

  return ret_value;
}

template <size_t degree, size_t modulus, class T>
bool run_p() {
  using poly_p = nfl::poly_p_from_modulus<T, degree, modulus>;

  // Stream
  std::stringstream ss;

  // Archives
  cereal::BinaryOutputArchive oarchive(ss); // output archive
  cereal::BinaryInputArchive iarchive(ss); // input archive

  // define a random polynomial
  poly_p p0{nfl::uniform()};

  // serialize the polynomial
  oarchive(p0);

  // define a new polynomial
  poly_p p1;

  // deserialize the stream
  iarchive(p1);

  // Verify that p0 == p1
  return (p0 == p1);
}

int main(int argc, char const* argv[]) {
  return not(run<CONFIG>() and run_p<CONFIG>());
}
