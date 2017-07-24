#ifndef fast_gauss_noise_h
#define fast_gauss_noise_h

#include <cstdint>
#include <cstddef>
#include <list>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <climits>
#include <cstring>
#include <tuple>
#include <typeinfo>
#include <gmp.h>
#include <mpfr.h>
#include "fastrandombytes.h"

#ifdef BOOST_RAPHSON
#include <boost/math/tools/roots.hpp>
#endif


//#define OUTPUT_BARRIERS
//#define OUTPUT_LUT_FLAGS
// Use it to test rare values do not happen too often
//#define UNITTEST_ONEMILLION

namespace nfl {

template<typename T0, typename T1>
constexpr inline auto tstbit(T0 x, T1 n) -> decltype((x << (63 - n)) >> 63 ) { return ((x << (63 - n)) >> 63 ); }

template<class in_class, class out_class>
struct output {
  out_class val;
  bool flag;
  std::list<in_class*> l_b_ptr;
};

template<class in_class, class out_class, unsigned _lu_depth>
class FastGaussianNoise {
    typedef output<in_class, out_class> output_t;
  private:
    unsigned int _bit_precision;
    unsigned int _word_precision;
    unsigned int _number_of_barriers;
    unsigned int _lu_size;
    unsigned int _flag_ctr1;
    unsigned int _flag_ctr2;
    double _sigma;
    mpfr_t _const_sigma;
    mpfr_t _center;
    int rounded_center;
    unsigned int _security;
    unsigned int _samples;  
    double _tail_bound;
    bool _verbose;

    in_class **barriers;
    output_t *lu_table;
    output_t **lu_table2;

    void check_template_params();
    void init();
    void precomputeBarrierValues(); 
    void buildLookupTables();
    void nn_gaussian_law(mpfr_t rop, const mpfr_t x_fr);
    int cmp(in_class *op1, in_class *op2); 

  public:
    static const unsigned int default_k;
    FastGaussianNoise(double sigma, unsigned int security, unsigned int samples, double center_d = 0, bool verbose = false);
    FastGaussianNoise(double sigma, unsigned int security, unsigned int samples, mpfr_t center, bool verbose = false);
    ~FastGaussianNoise();
    void getNoise(out_class * const rand_data2out, uint64_t rlen); 
};


/* NOTATIONS (from paper Sampling from discrete Gaussians 
 * for lattice-based cryptography on a constrainded device
 * Nagarjun C. Dwarakanath, Steven Galbraith
 * Journal : AAECC 2014 (25)
 *
 * m : lattice dimension, we set it to 1 in here (better results for same security) !!!!
 *
 * security : output distribution should be within 2^{-security} 
 * statistical distance from the Gaussian distribution
 *
 * sigma : parameter of the Gaussian distribution (close 
 * to the standard deviation but not equal)
 * SUM = \sum_{k=-\inf}^{\inf} exp(-(k-c)^2/(2*sigma^2))
 *
 * D_{sigma} : distribution such that for x\in \mathcal{Z}
 * Pr(x) = ro_{sigma}(x) = 1/SUM exp(-(x-c)^2/(2*sigma^2))
 *
 * Warning : can be replaced in some works by s = sqrt(2*pi)*sigma
 * if Pr(x) = 1/SUM exp(-pi*(x-c)^2/s^2)
 *
 * tail_bound : tail bound we only output points with 
 * norm(point-c) < tail_bound*sigma
 *
 * bit_precision : bit precision when computing the probabilities 
 *
 * delta_tailbound : statistical distance introduced 
 * by the tail bound
 *
 * delta_epsi : statistical distance introduced by 
 * probability approximation
 */



// HELPER FUNCTIONS

/* /!\ Warning only for x64 */
extern __inline__ uint64_t rdtsc(void) 
{
  uint64_t a, d;
  __asm__ volatile ("rdtsc" : "=a" (a), "=d" (d));
  return (d<<32) | a;
}

// Function to be used in a Newton-Raphson solver
struct funct
{
	funct(double const& target) : k(target){}
	std::tuple<double, double> operator()(double const& x)
	{
		return std::make_tuple(x*x - 2*log(x) - 1 - 2*k*log(2), 2*x-2/x); 
	}
	private:
	double k;
};

inline
double newton_raphson(double k, double max_guess, int digits)
{
  unsigned max_counter = 1U<<15;
  std::tuple<double, double> values;
  double delta;
  double guess = max_guess;
  for (unsigned counter = 0 ; counter < max_counter; counter++)
  {
    values = funct(k)(guess);
    delta = std::get<0>(values)/std::get<1>(values);
    guess -= delta;
    if ( fabs(delta)/fabs(guess) < pow(10.0,-digits) ) break;
  }
  // In case there is a flat zone in the function
  while(0.95*guess*0.95*guess - 2*log(0.95*guess) - 1 - 2*k*log(2)>=0) guess*=0.95;
  // Test result
  if(guess*guess - 2*log(guess) - 1 - 2*k*log(2)<0)
  {
    std::cout << "FastGaussianNoise: WARNING Newton-Raphson failed the generator is NOT secure" << std::endl;
  }
  return guess;
}


// Constructors


template<class in_class, class out_class, unsigned _lu_depth>
FastGaussianNoise<in_class, out_class, _lu_depth>::FastGaussianNoise( double sigma, 
		unsigned int security, unsigned int samples, double center_d /*=0*/, bool verbose /*=true*/):
	_sigma(sigma),
	_security(security),
  _samples(samples),
  _verbose(verbose)
{
  // Check template parameters
  check_template_params();
  
  //Center cannot be initialized before the constructor
  mpfr_init_set_d(_center, center_d, MPFR_RNDN);
  rounded_center = round(center_d);

  //Initialization functions
	init(); 
	precomputeBarrierValues();
	buildLookupTables();
}

template<class in_class, class out_class, unsigned _lu_depth>
FastGaussianNoise<in_class, out_class, _lu_depth>::FastGaussianNoise( double sigma, 
		unsigned int security, unsigned int samples, mpfr_t center /*=0*/, bool verbose /*=true*/):
	_sigma(sigma),
	_security(security),
  _samples(samples),
  _verbose(verbose)
{
  // Check template parameters
  check_template_params();

  //Center cannot be initialized before the constructor
  mpfr_init_set(_center, center, MPFR_RNDN);
  rounded_center = mpfr_get_d(_center, MPFR_RNDN);

  //Initialization functions
	init(); 
	precomputeBarrierValues();
	buildLookupTables();
}


// Check template parameters
template<class in_class, class out_class, unsigned _lu_depth>
void FastGaussianNoise<in_class, out_class, _lu_depth>::check_template_params() 
{
  // Lookup tables can only have depth 1 or 2
  if (_lu_depth != 1 && _lu_depth != 2)
  {
    std::cout << "FastGaussianNoise: CRITICAL _lu_depth must be 1 or 2" << std::endl;
    exit(1);
  }

  // Lookup tables can only have uint8_t or uint16_t indexes
  if (typeid(in_class) == typeid(uint8_t))
    _lu_size = 1<<8;
  else if(typeid(in_class) == typeid(uint16_t))
    _lu_size = 1<<16;
  else
  {
    std::cout << "FastGaussianNoise: CRITICAL in_class must be uint8_t or uint16_t" << std::endl;
    exit(2);
  }
}


// Compute some values (precision, number of barriers and outputs)
template<class in_class, class out_class, unsigned _lu_depth>
void FastGaussianNoise<in_class, out_class, _lu_depth>::init() 
{
	double epsi, k;

	/* Lemma 1: 
	 * tail_bound >= sqrt(1 + 2*log(tail_bound) + 2*k*log(2))  
	 * IMPLIES delta_tailbound < 2**(-k)
   *
   * -epsi <= -k - log2(2 * tail_bound * sigma)) 
	 * IMPLIES 2 * tail_bound * sigma * 2**(-epsi) <= 2**(-k)
	 * 
   * Lemma 2 (and Lemma 1): 
	 * 2 * tail_bound * sigma * 2**(-epsi) <= 2**(-k)
	 * IMPLIES  delta_epsilon < 2**(-k)
   */

	// THUS setting 
	k = _security + 1 + ceil(log(_samples)/log(2));
	// IMPLIES (delta_tailbound + delta_epsilon)*_samples < 2**(-_security)
	// We can thus generate vectors of _samples samples securely

  // To compute the tail bound we use Newton-Raphson with three digits of precision
	// WE HAVE tail_bound >= sqrt(1 + 2*log(tail_bound) + 2*k*log(2)) 
  // IS EQUIV TO tail_bound**2 -2*log(tail_bound) - 1- 2*k*log(2) >= 0
	int digits = 3;
#ifdef BOOST_RAPHSON
  double max_guess = 1+2*k*log(2);
  double min_guess = sqrt(1 + 2*k*log(2)), guess = min_guess;
	_tail_bound = boost::math::tools::newton_raphson_iterate(funct(k),guess, min_guess, max_guess, digits);
#else
  double min_guess = sqrt(1 + 2*k*log(2));
  _tail_bound = newton_raphson(k, min_guess, digits);
#endif
	// We now can compute the precision needed 
	epsi = k + log2(2 * _tail_bound * _sigma);
	_bit_precision = ceil(epsi);
	_word_precision = ceil(_bit_precision/(8.0*sizeof(in_class)));
  // For the same cost we get a simpler situation and more precision
	_bit_precision = _word_precision*8*sizeof(in_class);
  // And the number of probabilities that must be computed
  // Only half of them are computed (others are symmetric)
  // but all values are stored to speedup noise generation
  _number_of_barriers = 1+2*ceil(_tail_bound*_sigma);
  if ( ((uint64_t)_number_of_barriers >> (sizeof(out_class) * 8 - 1)) != 0)
    std::cout << "FastGaussianNoise: WARNING out_class too small to contain some of the (signed) results" << std::endl;
  if ( ((uint64_t)_number_of_barriers >> (sizeof(int) * 8 - 1)) != 0)
    std::cout << "FastGaussianNoise: WARNING outputs are above 2**31, unexpected results" << std::endl;

  // Finally we precompute 1/(2*sigma**2) to accelerate things
	mpfr_inits2(_bit_precision, _const_sigma, nullptr);
  mpfr_set_d(_const_sigma, _sigma, MPFR_RNDN);
	mpfr_sqr(_const_sigma, _const_sigma, MPFR_RNDN);
	mpfr_mul_ui(_const_sigma, _const_sigma, 2, MPFR_RNDN); 
	mpfr_ui_div(_const_sigma, 1, _const_sigma, MPFR_RNDN);

  // Give some feedback
  if (_verbose) std::cout << "FastGaussianNoise: " << _number_of_barriers << 
    " barriers with " << _bit_precision << 
    " bits of precision will be computed" << std::endl;
}


template<class in_class, class out_class, unsigned _lu_depth>
void FastGaussianNoise<in_class, out_class, _lu_depth>::precomputeBarrierValues() 
{
  // Declare and init mpfr vars
  mpfr_t sum, tmp, tmp2;
  mpfr_t *mp_barriers;
  mpfr_inits2(_bit_precision, sum, tmp, tmp2, nullptr);

  // This var is used to export mpfr values
  mpz_t int_value;
  mpz_init2(int_value, _bit_precision);

  // Init SUM = \sum_{k=-ceil(tail_bound*sigma)}^{ceil(tail_bound*sigma)} exp(-(k+round(c)-c)^2/(2*sigma^2))
  // and compute on the loop with the barriers
  mpfr_set_ui(sum, 0, MPFR_RNDN);

  // Allocate memory for the barrier pointers
  barriers = (in_class **) malloc(_number_of_barriers*sizeof(in_class *)); 
  mp_barriers = (mpfr_t *) malloc(_number_of_barriers*sizeof(mpfr_t)); 

  // Now loop over the barriers
  for (int i = 0; i < (int)_number_of_barriers; i++)
  {
    // Init mpfr var
    mpfr_init2(mp_barriers[i], _bit_precision);

    // Compute the barrier value (without normalization)
    mpfr_set_si(tmp2, rounded_center+i-((int)_number_of_barriers-1)/2, MPFR_RNDN);
    nn_gaussian_law(tmp, tmp2);
    //mpfr_out_str(stdout, 10, 0, tmp2, MPFR_RNDN);
    //std::cout << std::endl;
    if (i==0) mpfr_set(mp_barriers[0], tmp, MPFR_RNDN);
    else 
    {
      mpfr_add(mp_barriers[i], mp_barriers[i-1], tmp, MPFR_RNDN);
    }

    // Add the probability to the sum
    mpfr_add(sum, sum, tmp, MPFR_RNDN);
  }	
  
  // Invert the sum and scale it 
  mpfr_ui_div(sum, 1, sum, MPFR_RNDN);
  mpfr_set_ui(tmp, 2, MPFR_RNDN);
  mpfr_pow_ui(tmp, tmp, _bit_precision, MPFR_RNDN);
  mpfr_sub_ui(tmp, tmp, 1, MPFR_RNDN);
  mpfr_mul(sum, sum, tmp, MPFR_RNDN);

  // Now that we got the inverted sum normalize and export
  for (unsigned i = 0; i < _number_of_barriers; i++)
  {
    // Allocate space
    barriers[i] = (in_class *) calloc(_word_precision,sizeof(in_class)); 
	  
    mpfr_mul(mp_barriers[i], mp_barriers[i], sum, MPFR_RNDN);  
    mpfr_get_z(int_value, mp_barriers[i], MPFR_RNDN);
    mpz_export((void *) (barriers[i] + ((int)_word_precision - 
           (int)ceil( (float)mpz_sizeinbase(int_value, 256)/sizeof(in_class) ))), nullptr, 1, sizeof(in_class), 0, 0, int_value);
#ifdef OUTPUT_BARRIERS
    mpz_out_str(stdout, 10, int_value);
    std::cout << " = Barriers[" << i << "] = " << std::endl;
    if (sizeof(in_class) == 1) for (unsigned j = 0 ; j < _word_precision; j++)
      printf("%.2x", barriers[i][j]);
    if (sizeof(in_class) == 2) for (unsigned j = 0 ; j < _word_precision; j++)
      printf("%.4x", barriers[i][j]);
    std::cout <<  std::endl;
#endif 
    mpfr_clear(mp_barriers[i]);
  }
	mpfr_clears(sum, tmp, tmp2, nullptr);
  mpz_clear(int_value);
	mpfr_free_cache();
  free(mp_barriers);
}



//Build lookup tables used during noise generation
template<class in_class, class out_class, unsigned _lu_depth>
void FastGaussianNoise<in_class, out_class, _lu_depth>::buildLookupTables() 
{
	unsigned lu_index1 = 0, lu_index2 = 0;
  _flag_ctr1 = _flag_ctr2 = 0;

  // Allocate space for the lookup tables
  lu_table = new output_t[_lu_size](); 
  if(_lu_depth == 2) 
    lu_table2 = (output_t **) calloc(_lu_size,sizeof(output_t *)); 

	// We start building the first dimension of the lookup table
  // corresponding to the first in_class word of the barriers
	for (int64_t val = -((int)_number_of_barriers-1)/2 + rounded_center, b_index = 0; val <= ((int)_number_of_barriers-1)/2 + rounded_center && lu_index1 < _lu_size;) 
  {

		while (lu_index1 < barriers[b_index][0] && lu_index1 < _lu_size) 
    {
      lu_table[lu_index1].val = val;
			lu_index1++;
		}

		// Flag the entry
		lu_table[lu_index1].val = val;
		lu_table[lu_index1].flag = true;
    _flag_ctr1++;
		// If _lu_depth == 1 we have to list the barriers here
    if (_lu_depth == 1)
    {
#ifdef OUTPUT_LUT_FLAGS
      std::cout << "FastGaussianNoise: flagged lu_table[" 
       << lu_index1 << "] for barriers " << val;
#endif
      
      // Prepare the first element of the chained list of barrier
			lu_table[lu_index1].l_b_ptr.push_back(barriers[b_index++]);
      val++;
			// If more that one barrier is present, we add them to the chained list 
			while ( (b_index<_number_of_barriers) && (lu_index1 == barriers[b_index][0])) 
      {
			  lu_table[lu_index1].l_b_ptr.push_back(barriers[b_index]);
#ifdef OUTPUT_LUT_FLAGS
      std::cout << "FastGaussianNoise: flagged lu_table[" << lu_index1 << "] for barriers " << val;
#endif
				b_index++;
				val++;
			} // while
    } // if
	
    if (_lu_depth == 2)
    {
      // When we meet a barrier in an entry of the lu_table, 
      // we build another lu_table inside that entry
		  // corresponding to the next in_class word of the barriers 
		  lu_index2 = 0;
      lu_table2[lu_index1] = new output_t[_lu_size](); 
		  while (lu_index2 < _lu_size) 
      {
			  if(lu_index1 < barriers[b_index][0] || lu_index2 < barriers[b_index][1])
        {
          lu_table2[lu_index1][lu_index2].val = val;
        }
        else
        {
		      // If we are on a barrier
		      if (lu_index1 == barriers[b_index][0] && lu_index2 == barriers[b_index][1]) 
          {
			      // Flag the entry
			      lu_table2[lu_index1][lu_index2].val = val;
			      lu_table2[lu_index1][lu_index2].flag = true;
#ifdef OUTPUT_LUT_FLAGS
            std::cout << "FastGaussianNoise: flagged lu_table2[" << lu_index1 << "][" << lu_index2 << "] for barriers " << val;
#endif
            _flag_ctr2++;
			      // And prepare the first element of the chained list of barrier
			      lu_table2[lu_index1][lu_index2].l_b_ptr.push_back(barriers[b_index++]);
            val++;
			      // If more that one barrier is present, we add them to the chained list 
			      while ( (b_index<_number_of_barriers) && 
                (lu_index1 == barriers[b_index][0]) && 
                (lu_index2 == barriers[b_index][1]) ) 
            {
				      lu_table2[lu_index1][lu_index2].l_b_ptr.push_back(barriers[b_index]);
#ifdef OUTPUT_LUT_FLAGS
            std::cout << " " << val;
#endif
				      b_index++;
				      val++;
			      } // while
#ifdef OUTPUT_LUT_FLAGS
            std::cout << std::endl;
#endif
          } // if
        } // else
        lu_index2++;
      } // while
		} // if

		lu_index1++;
	}
  // Give some feedback
  if (_verbose) std::cout << "FastGaussianNoise: Lookup tables built" << std::endl;
}

template<class in_class, class out_class, unsigned _lu_depth>
void FastGaussianNoise<in_class, out_class, _lu_depth>::getNoise(out_class* const rand_outdata, uint64_t rlen) 
{
	uint64_t computed_outputs, innoise_bytesize, innoise_words, used_words;
	int64_t output;
  bool flagged;
	in_class *noise, *noise_init_ptr, input1, input2;
  float innoise_multiplier;
 
  // Expected number of input bytes per output byte. Lowering this 
  // (e.g. by a factor .5) works but could lead to segfaults.
  if (_lu_depth == 1)
  {
    innoise_multiplier = 1.05 * ((float)(_lu_size - _flag_ctr1)/(float)_lu_size) + _word_precision * ((float)_flag_ctr1/_lu_size);
  } 
  else // _lu_depth == 2
  {
    innoise_multiplier = 1.05 * ((float)(_lu_size - _flag_ctr1)/(float)_lu_size) + 2.0 * ((float)_flag_ctr1/(float)_lu_size) + _word_precision * ((float)_flag_ctr2/((float)_lu_size*_lu_size));
  }
  innoise_words = rlen * innoise_multiplier;
  innoise_bytesize = sizeof(in_class) * innoise_words;
	noise = noise_init_ptr = new in_class[innoise_words];
  used_words = 0;
  if (_verbose) std::cout << "FastGaussianNoise: Using " << " " <<innoise_multiplier*sizeof(in_class)*8/std::min(log2(_number_of_barriers),(double)sizeof(out_class)*8) << " input bits per output bit" << std::endl;
 
  // Count time for uniform noise generation 
	uint64_t start = rdtsc();
	fastrandombytes((uint8_t*)noise, innoise_bytesize);
	uint64_t stop = rdtsc();

  // Give some feedback
	if (_verbose) printf("FastGaussianNoise: Uniform noise  cycles = %.2e bits = %.2e cycles/bit = %.2e\n", (float) stop - start, (float) innoise_bytesize * 8, (float)(stop-start)/(innoise_bytesize*8));

  // Loop until all the outputs have been generated
  computed_outputs = 0;
	while (computed_outputs < rlen ) 
  {
		input1 = *noise;
		flagged = lu_table[input1].flag;

    // If flagged we have to look at the next in_class word
    if (flagged) 
    {
		  if (_lu_depth == 1)
      {
				output = lu_table[input1].val;
				for(in_class* b_ptr : lu_table[input1].l_b_ptr) 
        {
					// If the barrier value is greater than the noise
					if(cmp(b_ptr, noise)==1) break;
					output++;
				}
        // We shift the noise pointer of word_precision minus 1
        // As there another byte shift later
				noise += _word_precision - 1;
				used_words += _word_precision - 1;

      }
      else // _lu_depth == 2
      {
			  input2 = *(noise+1);
			  flagged = lu_table2[input1][input2].flag;
        // If flagged again we compare using full precision the random value 
        // with all barriers in the linked list of the lookup table
			  if (flagged) 
        {
			    output = lu_table2[input1][input2].val;
				  for(in_class* b_ptr : lu_table2[input1][input2].l_b_ptr) 
          {
					  // If the barrier value is greater than the noise
					  if(cmp(b_ptr, noise)==1) break;
					  output++;
				  }
          // We shift the noise pointer of word_precision minus 2
          // As there are two other one byte shifts later
				  noise += _word_precision - 2;
				  used_words += _word_precision - 2;
			  } // if                                     		
        else
        {
			    output = lu_table2[input1][input2].val;
        }
        noise++;
        used_words++;
      } // else
    }
    else
    { 
		  output = lu_table[input1].val;
    }
		noise++;
    used_words++;
		// Add the obtained result to the list of outputs
		rand_outdata[computed_outputs++] = (out_class) output; 

#ifdef UNITTEST_ONEMILLION
    if ( (output > _sigma * 6) || (output < -_sigma*6) )
    {
      std::cout << output << "FastGaussianNoise: Unit test failed, this should happen once in a million. Uniform input leading to this is  ";
      for (unsigned i = 0; i < _word_precision ; i++)
        printf("%.2x", *(noise-_word_precision+i));
      std::cout << std::endl;
    }
#endif

#if 1
    // If too much noise has been used regenerate it
    if ((used_words+_word_precision) >= innoise_words)
    {
      noise = noise_init_ptr;
      used_words = 0; 
      if (_verbose) std::cout << "FastGaussianNoise: All the input bits have been used, regenerating them ..." << std::endl;
 
	    fastrandombytes((uint8_t*)noise, innoise_bytesize);
    }
#endif
	}
	delete[] noise_init_ptr;
}

/* Compare two arrays word by word.
 * return 1 if op1 > op2, 0 if equals and -1 if op1 < op2 */
template<class in_class, class out_class, unsigned _lu_depth>
inline int FastGaussianNoise<in_class, out_class, _lu_depth>::cmp(in_class *op1, in_class *op2) 
{

	for (int i = 0; i < (int)_word_precision; i++) 
  {

		if (op1[i] > op2[i]) return 1;

		else if (op1[i] < op2[i]) return -1;
	}
	return 0;
}


// Compute exp(-(x-center)^2/(2*sigma^2)) this is not normalized ! (hence the nn)
template<class in_class, class out_class, unsigned _lu_depth>
void  inline FastGaussianNoise<in_class, out_class, _lu_depth>::nn_gaussian_law(mpfr_t rop, const mpfr_t x) 
{
	mpfr_sub(rop, x, _center, MPFR_RNDN);
	mpfr_sqr(rop, rop, MPFR_RNDN);
	mpfr_neg(rop, rop, MPFR_RNDN);
	mpfr_mul(rop, rop, _const_sigma, MPFR_RNDN);
	mpfr_exp(rop, rop, MPFR_RNDN);
}


template<class in_class, class out_class, unsigned _lu_depth>
FastGaussianNoise<in_class, out_class, _lu_depth>::~FastGaussianNoise() 
{
  // Freed allocated memory for the barriers
  for (unsigned ctr = 0; ctr < _number_of_barriers; ctr++)
  {
    if (barriers[ctr] != nullptr) free(barriers[ctr]); 
    barriers[ctr] = nullptr;
  }
  if (barriers != nullptr) free(barriers);
  barriers=nullptr;

  // Free other variables
  mpfr_clear(_const_sigma);
  mpfr_clear(_center);
  delete[](lu_table);
  if(_lu_depth == 2) 
  {
    for (unsigned ctr = 0 ; ctr < _lu_size; ctr++)
    {
      if (lu_table2[ctr]!=nullptr) delete[](lu_table2[ctr]);
    }
    free(lu_table2);
  }
}

}  // namespace nfl

#endif
