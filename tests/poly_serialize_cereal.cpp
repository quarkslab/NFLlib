#include <cereal/archives/binary.hpp>
#include <iostream>
#include <nfl.hpp>
#include <sstream>
#include "tools.h"

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

int main(int argc, char const* argv[]) {
  return not(run<CONFIG>() and run_p<CONFIG>());
}
