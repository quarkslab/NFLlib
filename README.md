# NFLlib
An NTT-based Fast Lattice library

## Goal

NFLlib is an efficient and open-source C++ library dedicated to ideal lattice cryptography. It is specialized in polynomial rings quotiented by a cyclotomic polynomial whose degree is a power of two. The library combines algorithmic optimizations (Chinese Remainder Theorem, optimized Number Theoretic Transform) together with programming optimization techniques (SSE and AVX2 specializations, C++ expression templates, etc.).

## License

MIT

# Install Steps

You need cmake, GMP and Mpfr, as well as a C++11 compiler to build NFLLib.

To build, test and install a production version of nfllib, run the following:

```
$> mkdir _build
$> cd _build
$> cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/nfllib
$> make
$> make test
$> make install
```

The following CMake options are relevant:

Option                              | Description
------------------------------------|---------------------------------
`-DCMAKE_INSTALL_PREFIX=<value>`    | Where the library is installed
`-DCMAKE_BUILD_TYPE=Release\|Debug` | The basic compiler configuration
`-DNFL_OPTIMIZED=ON`                | Enable SSE/AVX-based optimization

# SSE/AVX optimizations

To use SSE-based optimizations, compile the code with the flags `-DNFL_OPTIMIZED=ON -DNTT_SSE`.

To use AVX-based optimizations, compile the code with the flags `-DNFL_OPTIMIZED=ON -DNTT_AVX2`.

# Getting started

In order to get the documentation you need Sphinx and the Alabaster theme. Installation procedures are described at (we recommend using pip for both installations) :

  * http://www.sphinx-doc.org/en/stable/install.html
  * https://anaconda.org/pypi/alabaster

After getting both Sphinx and the Alabaster theme build the documentation and show it with:

```
$> sphinx-build . build
$> your_favorite_browser build/nfl.html
```

If you have issues building the documentation please contact the developper team.

# Contributors

This library is an extension/evolution of the NTTTools module from [XPIR](https://github.com/XPIR-team/XPIR) done by members of [CryptoExperts](https://www.cryptoexperts.com), [INP ENSEEIHT](http://www.enseeiht.com), [Quarkslab](http://www.quarkslab.com) (in alphabetical order).
