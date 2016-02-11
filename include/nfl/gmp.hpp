#ifndef NFL_GMP_HPP
#define NFL_GMP_HPP

#include "nfl/poly.hpp"

namespace nfl
{

template<class T, size_t Degree, size_t NbModuli>
typename poly<T, Degree, NbModuli>::GMP poly<T, Degree, NbModuli>::gmp;

template<class T, size_t Degree, size_t NbModuli>
poly<T, Degree, NbModuli>::poly(mpz_t v) {
  set(v);
}

template<class T, size_t Degree, size_t NbModuli>
poly<T, Degree, NbModuli>::poly(mpz_class v) {
  set(v);
}

template<class T, size_t Degree, size_t NbModuli>
poly<T, Degree, NbModuli>::poly(mpz_t* values) {
  set(values);
}

template<class T, size_t Degree, size_t NbModuli>
poly<T, Degree, NbModuli>::poly(mpz_class* values) {
  set(values);
}

template<class T, size_t Degree, size_t NbModuli>
void poly<T, Degree, NbModuli>::set(mpz_t v)
{
  auto* iter = begin();
  
  for(size_t cm = 0; cm < nmoduli; ++cm) 
  {
    *iter++ = mpz_fdiv_ui(v, params<T>::P[cm]);
    for(size_t i = 1; i < degree; i++)
    {  
      *iter++ = 0;
    }
  }
}

template<class T, size_t Degree, size_t NbModuli>
void poly<T, Degree, NbModuli>::set(mpz_class v)
{
  set(v.get_mpz_t());
}

template<class T, size_t Degree, size_t NbModuli>
void poly<T, Degree, NbModuli>::set(mpz_t* values)
{
	*this = gmp.mpz2poly(values);
}

template<class T, size_t Degree, size_t NbModuli>
void poly<T, Degree, NbModuli>::set(mpz_class* values)
{
	mpz_t* values_mpz_t = new mpz_t[nmoduli*degree];
	for (size_t i=0; i<nmoduli*degree; i++) {
		mpz_init_set(values_mpz_t[i], values[i].get_mpz_t());
	}
	*this = gmp.mpz2poly(values_mpz_t);
	for (size_t i=0; i<nmoduli*degree; i++) {
		mpz_clear(values_mpz_t[i]);
	}
}

template<class T, size_t Degree, size_t NbModuli>
poly<T, Degree, NbModuli>::GMP::GMP()
{
	// Compute product of moduli
	mpz_init_set_ui(moduli_product, 1);
	for (unsigned i = 0; i < nmoduli; i++)
	{
		mpz_mul_ui(moduli_product, moduli_product, params<T>::P[i]);
	}

	bits_in_moduli_product = mpz_sizeinbase(moduli_product, 2);

	// Compute Shoup value for optimized reduction modulo "moduli_product"
	mpz_t one;
	mpz_init_set_ui(one, 1);
	mpz_init(modulus_shoup);
	mpz_mul_2exp(modulus_shoup, one, bits_in_moduli_product + sizeof(T)*CHAR_BIT + params<T>::logkMaxNbModuli + 1);
	mpz_tdiv_q(modulus_shoup, modulus_shoup, moduli_product);

	// Compute the lifting coefficients
	mpz_t quotient;
	mpz_init(quotient);
	lifting_integers = (mpz_t*) malloc(nmoduli * sizeof(mpz_t));

	for (unsigned i = 0; i < nmoduli; i++)
	{
		// Current modulus
		mpz_t current_modulus;
		mpz_init_set_ui(current_modulus, params<T>::P[i]);

		// compute the product of primes except the current one
		mpz_divexact(quotient, moduli_product, current_modulus);

		// Compute the inverse of the product
		mpz_init(lifting_integers[i]);
		mpz_invert(lifting_integers[i], quotient, current_modulus);

		// Multiply by the quotient
		mpz_mul(lifting_integers[i], lifting_integers[i], quotient);
		mpz_clear(current_modulus);
	}

	// Clear
	mpz_clear(one);
	mpz_clear(quotient);
}

template<class T, size_t Degree, size_t NbModuli>
poly<T, Degree, NbModuli>::GMP::~GMP()
{
	for (unsigned i = 0; i < nmoduli; i++)
	{
		mpz_clear(lifting_integers[i]);
	}
	free(lifting_integers);
	mpz_clear(modulus_shoup);
	mpz_clear(moduli_product);
}

template<class T, size_t Degree, size_t NbModuli>
mpz_t* poly<T, Degree, NbModuli>::GMP::poly2mpz(const poly& op)
{
	mpz_t* poly_mpz;
	poly_mpz = new mpz_t[degree];
	for (unsigned i = 0; i < degree; i++)
	{
		mpz_init(poly_mpz[i]);
	}

	poly2mpz(poly_mpz, op);
	return poly_mpz;
}

template<class T, size_t Degree, size_t NbModuli>
void poly<T, Degree, NbModuli>::GMP::poly2mpz(mpz_t* rop, const poly& op)
{
	mpz_t tmp;
	mpz_init2(tmp, (bits_in_moduli_product<<2));
	for (unsigned i = 0; i < degree; i++)
	{
		mpz_set_ui(rop[i], 0);
		for (unsigned cm = 0; cm < nmoduli; cm++)
		{
			if (op(cm, i) != 0)
			{
				mpz_mul_ui(tmp, lifting_integers[cm], op(cm, i));
			}
			mpz_add(rop[i], rop[i], tmp);
		}
		// Modular reduction using Shoup
		mpz_mul(tmp, rop[i], modulus_shoup);
		mpz_tdiv_q_2exp(tmp, tmp, bits_in_moduli_product + sizeof(T)*CHAR_BIT + params<T>::logkMaxNbModuli + 1);
		mpz_mul(tmp, tmp, moduli_product);
		mpz_sub(rop[i], rop[i], tmp);
		if (mpz_cmp(rop[i], moduli_product) >= 0)
		{
			mpz_sub(rop[i], rop[i], moduli_product);
		}
	}
	mpz_clear(tmp);
}

template<class T, size_t Degree, size_t NbModuli>
poly<T, Degree, NbModuli> poly<T, Degree, NbModuli>::GMP::mpz2poly(const mpz_t* poly_mpz)
{
	poly result;
	mpz2poly(result, poly_mpz);
	return result;
}

template<class T, size_t Degree, size_t NbModuli>
void poly<T, Degree, NbModuli>::GMP::mpz2poly(poly<T, Degree, NbModuli>& rop, const mpz_t* poly_mpz)
{
	for (unsigned cm = 0; cm < nmoduli; cm++)
	{
		for (unsigned i = 0; i < degree; i++)
		{
			rop(cm, i) = mpz_fdiv_ui(poly_mpz[i], params<T>::P[cm]);
		}
	}
}

}

#endif