/*
 * File:   lattisigns512-20130329/fastrandombytes.c
 * Author: Gim Güneysu, Tobias Oder, Thomas Pöppelmann, Peter Schwabe
 * Public Domain
 */

#include <inttypes.h>
#include <iostream>
#include "nfl/prng/crypto_stream_salsa20.h"
#include "nfl/prng/randombytes.h"

namespace nfl {

static int init = 0;
static unsigned char key[32];
static unsigned char nonce[8] = {0};

void fastrandombytes(unsigned char *r, unsigned long long rlen) {
  unsigned long long n = 0;
  int i;
  if (!init) {
    randombytes(key, 32);
    init = 1;
  }
  nfl_crypto_stream_salsa20_amd64_xmm6(r, rlen, nonce, key);

  // Increase 64-bit counter (nonce)
  for (i = 0; i < 8; i++) n ^= ((unsigned long long)nonce[i]) << 8 * i;
  n++;
  for (i = 0; i < 8; i++) nonce[i] = (n >> 8 * i) & 0xff;
}
}