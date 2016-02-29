.. nfllib documentation master file, created by
   sphinx-quickstart on Thu Mar 19 09:32:40 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to nfllib's documentation!
==================================

Contents:

.. toctree::
   :maxdepth: 2


Install
-------

You need cmake, GMP and Mpfr, as well as a C++11 compiler to build NFLLib.

To build, test and install a production version of nfllib, run the following:

.. code:: sh

    $> mkdir _build
    $> cd _build
    $> cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/nfllib
    $> make
    $> make test
    $> make install

Polynomials
-----------

Polynomials are the core of NFLlib. In NFLlib, a polynomial is described as a triplet:

:type:

    The integer type used to actually store the polynomial. Must be uint16_t, uint32_t or uint64_t.

:degree:

    The degree of the polynomial. Must be a power of two.

:modulus_bitsize:

    The number of bits of its coefficients. Must be equal to the bit size of the type used to store the polynomial minus two, times an integer (e.g. if type is uint32_t, modulus_bitsize can be 30, 60, 90, etc.)   

Programatically, you can retreive a given polynomial type using the type generator ``nfl::poly_from_modulus<type, degree, modulus_bitsize>``, for instance using a type alias:

.. code-block:: c++

    #include "nfl.hpp"
    using poly_type = nfl::poly_from_modulus<uint32_t, 512, 30>

which defines a polynomial type of degree 512, with 30 bits for each
coefficient, stored on unsigned 32 bits integers. Some parameters combination
are invalid, in that case you should receive a compile-time error!

Initialization
==============

.. highlight:: c++

Polynomials can be initialized in multiple ways:

:null polynomial:

    The default constructor builds a polynomial with all coefficients set to zero:: 

        poly_type P; // the null polynomial P(X) = 0

:constant polynomial:

    One can build a constant polynomial by passing an integral value to the constructor::

        poly_type P(1); // constant polynomial P(X) = 1

:any polynomial:

    One can specify the value of each coefficient, higher degrees default to zero::

        poly_type P{1,2,3}; // P(X) = 1 + 2X + 3X^2

:uniform polynomial:

    A polynomial initialized with an uniform distribution::

        poly_type P(nfl::uniform());

:non-uniform polynomial:

    A polynomial initialized with a bounded uniform distribution::

        uint32_t constexpr bound = 8;
        poly_type P(nfl::non_uniform(8));

:gaussian polynomial:

    A polynomial initialized with a gaussian distribution::

        // Define the gaussian generator object (parameters are explained below)
        nfl::FastGaussianNoise<uint8_t, typename poly_type::value_type, 2> g_prng(8, 128, 1<<15);
        // Define a struct to access the generator via a constructor (parameters are explained below)
        nfl::gaussian<typename poly_type::value_type> gaussian_struct(&g_prng, 2);
        // This initializes P as a polynomial with coefficients following closely a discrete gaussian in 
        // the sense of [1] with standard deviation of 8. The output distribution is within a statistical 
        // distance of the target distribution of 2^{-128} as long as the polynomial degree is below 2^{15}.
        // As usual in cryptography, the coefficients of this noise polynomial have been multiplied by 2 as
        // requested by the second parameter of the struct constructor.
        // [1]: Sampling from discrete Gaussians for lattice-based cryptography on a constrained device
        // Nagarjun C. Dwarakanath, Steven Galbraith, AAECC 2014 (25)

        poly_type P(gaussian_struct);

Gaussian generators
===================

Polynomials with coefficients following a discrete gaussian are needed in most cryptographic applications. These are obtained using a gaussian generator which is defined by two triplets::
        
        nfl::FastGaussianNoise<internal_array_index_type,  output_type, recursion_level> gaussian_prng_instance(sigma, security, max_degree);


The first triplet ``<internal_array_type, output_type, recursion_level>`` defines the type through a template and exists only to improve performance using some parameters at compilation time. The second triplet ``(sigma, security, max_degree)`` is used by the constructor to define the generator. 

 

The different parameters are used as follows:

:internal_array_index_type:

    May be uint8_t or uint16_t, the latter type can improve computational performance but increases memory usage as the internal lookup tables used are larger.

:output_type:

    The type used to store the output noise. Note that if the noise generated is signed and if the output type is unsigned the sign is just removed. If the application context requires working modulo a number it is thus best to retrieve a signed result and add the modulus when needed to get a positive result.

:recursion_level:

    May be 1 or 2. Recursion level 2 may improve computational performance but increases memory usage. 

:sigma:

    Standard deviation (in the sense of [1]) of the discrete gaussian from which are sampled the polynomial coefficients.

    [1]: Sampling from discrete Gaussians for lattice-based cryptography on a constrained device Nagarjun C. Dwarakanath, Steven Galbraith, AAECC 2014 (25)
 
:security:

    The statistical distance between the output distribution and the target distribution is below 2^{-security}  as long as the degree of the polynomial is not to large.

:degree:
  
    Maximum degree for a polynomial to respect the statistical distance requested by the parameter ``security``.


Operations on Polynomials
=========================

Once a polynomial is created, you can do basic arithmetic on them, as in::

    poly_type P = P0 + P1 * P2;

The library optimizes these operations to avoid the creation of temporary polynomials and vectorizes the operations.

.. Concepts
.. --------
.. 
.. :math:`a^2 + b^2 = c^2`
.. 
.. .. math::
.. 
..    (a + b)^2 = a^2 + 2ab + b^2
.. 
..    (a - b)^2 = a^2 - 2ab + b^2
.. 
.. API
.. ---
.. 
.. Add
.. ***
.. 
.. .. doxygenfunction:: nfl::add
.. 
.. 
.. Sub
.. ***
.. 
.. .. doxygenfunction:: nfl::sub
.. 
.. Operators
.. *********
.. 
.. Operator versions of ``add``, ``sub`` and the likes using expression templates to avoid the creation of useless intermediate polynomials. It takes a polynomial or a polynomial expression as arguments and returns a polynomial expression that gets evaluated when assigned to a ``poly<T,D,N>``.
.. 
.. Both version internally use the ``addmod`` functor to perform the computation.
.. 
.. 
.. 
.. 
.. 
.. 
.. 
.. License
.. -------
.. 
.. 
.. 
.. 
.. 
.. Indices and tables
.. ==================
.. 
.. * :ref:`genindex`
.. * :ref:`search`
.. 
