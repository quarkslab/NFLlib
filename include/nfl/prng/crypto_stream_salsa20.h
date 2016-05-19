#ifndef CRYPTO_STREAM_H
#define CRYPTO_STREAM_H

#include <string>

extern "C" {
int nfl_crypto_stream_salsa20_amd64_xmm6(unsigned char *c, unsigned long long clen,
                                     const unsigned char *n,
                                     const unsigned char *k);
}

#endif
