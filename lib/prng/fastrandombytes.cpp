/*
 * File:   lattisigns512-20130329/fastrandombytes.c
 * Author: Gim Güneysu, Tobias Oder, Thomas Pöppelmann, Peter Schwabe
 * Public Domain
 */

#include <inttypes.h>
#include <iostream>
#include "nfl/prng/crypto_stream_salsa20.h"
#include "nfl/prng/randombytes.h"

#define nfl_crypto_stream_salsa20_KEYBYTES 32
#define nfl_crypto_stream_salsa20_NONCEBYTES 8

namespace nfl {

static int init = 0;
static unsigned char key[nfl_crypto_stream_salsa20_KEYBYTES];
static unsigned char nonce[nfl_crypto_stream_salsa20_NONCEBYTES] = {0};

void fastrandombytes(unsigned char *r, unsigned long long rlen) {
  unsigned long long n = 0;
  int i;
  if (!init) {
    randombytes(key, nfl_crypto_stream_salsa20_KEYBYTES);
    init = 1;
  }
  nfl_crypto_stream_salsa20_amd64_xmm6(r, rlen, nonce, key);

  // Increase 64-bit counter (nonce)
  for (i = 0; i < nfl_crypto_stream_salsa20_NONCEBYTES; i++) n ^= ((unsigned long long)nonce[i]) << 8 * i;
  n++;
  for (i = 0; i < nfl_crypto_stream_salsa20_NONCEBYTES; i++) nonce[i] = (n >> 8 * i) & 0xff;
}
}