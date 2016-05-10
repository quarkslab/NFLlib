#ifndef CRYPTO_STREAM_H
#define CRYPTO_STREAM_H

#ifdef __cplusplus
#include <string>
#endif

namespace nfl
{
#define crypto_stream_salsa20_KEYBYTES 32
#define crypto_stream_salsa20_NONCEBYTES 8

#ifdef __cplusplus
extern "C" {
#endif
extern
int crypto_stream_salsa20_amd64_xmm6(
        unsigned char *c,unsigned long long clen,
  const unsigned char *n,
  const unsigned char *k
);
    
#ifdef __cplusplus
}
#endif
#define crypto_stream_salsa20 crypto_stream_salsa20_amd64_xmm6

}

#endif
