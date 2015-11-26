#ifndef NFL_ARCH_HPP
#define NFL_ARCH_HPP

#include "nfl/arch/common.hpp"

#ifdef NFL_OPTIMIZED

#if defined __AVX2__ && defined NTT_AVX2
#define CC_SIMD nfl::simd::avx2
#elif defined __SSE4_2__ && defined NTT_SSE
#define CC_SIMD nfl::simd::sse
#endif

#endif

#ifndef CC_SIMD
#define CC_SIMD nfl::simd::serial
#endif

#endif
