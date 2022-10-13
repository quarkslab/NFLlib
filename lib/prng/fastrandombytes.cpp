/*
 * File:   lattisigns512-20130329/fastrandombytes.c
 * Author: Gim Güneysu, Tobias Oder, Thomas Pöppelmann, Peter Schwabe
 * Public Domain
 */

#include <inttypes.h>
#include <string.h>
#include <iostream>
#include "nfl/prng/crypto_stream_salsa20.h"
#include "nfl/prng/randombytes.h"

namespace nfl {

static size_t constexpr crypto_stream_salsa20_KEYBYTES = 32;
static size_t constexpr crypto_stream_salsa20_NONCEBYTES = 8;

static int init = 0;
static unsigned char key[crypto_stream_salsa20_KEYBYTES];
static unsigned char nonce[crypto_stream_salsa20_NONCEBYTES] = {0};

void fastrandombytes_seed(unsigned char *s, unsigned long long slen) {
    memcpy(key, s, slen);
    memset(key + slen, 0, crypto_stream_salsa20_KEYBYTES - slen);
    memset(nonce, 0, crypto_stream_salsa20_NONCEBYTES);
    init = -1;
}

void fastrandombytes_reseed() {
    init = 0;
}

void fastrandombytes(unsigned char *r, unsigned long long rlen) {
  unsigned long long n = 0;
  int i;
  if (!init) {
    randombytes(key, crypto_stream_salsa20_KEYBYTES);
    init = 1;
  }
  nfl_crypto_stream_salsa20_amd64_xmm6(r, rlen, nonce, key);

  // Increase 64-bit counter (nonce)
  for (i = 0; i < crypto_stream_salsa20_NONCEBYTES; i++) n ^= ((unsigned long long)nonce[i]) << 8 * i;
  n++;
  for (i = 0; i < crypto_stream_salsa20_NONCEBYTES; i++) nonce[i] = (n >> 8 * i) & 0xff;
}
}
