#include "nfl/prng/FastGaussianNoise.hpp"

#define OUTPUT_NOISE_SIZE 50000000LLU
#define REPETITIONS 10

int main(int argc, const char *argv[])
{
	typedef int32_t T;
	//FastGaussianNoise<uint16_t, T, 1> rng(300, 128, 1UL<<10, true);
	nfl::FastGaussianNoise<uint8_t, T, 2> rng(3.19, 128, 1UL << 19, 5.25, true);

	T *noise = new T[OUTPUT_NOISE_SIZE]();
	bzero(noise, OUTPUT_NOISE_SIZE);

	for (unsigned i = 0; i < REPETITIONS ; i++)
	{
		uint64_t start = nfl::rdtsc();
		rng.getNoise(noise, OUTPUT_NOISE_SIZE);
		uint64_t stop = nfl::rdtsc();

	  printf("FastGaussianNoise: Gaussian noise cycles = %.2e bits = %.2e cycles/bit = %.2e\n", (float) stop - start, (float) (OUTPUT_NOISE_SIZE*CHAR_BIT), (float)(stop - start) / (OUTPUT_NOISE_SIZE*CHAR_BIT));
  } 
#if(1)
  std::cout << "FastGaussianNoise: Noise generated below" << std::endl;
	std::cout << "[";
	for (unsigned int i = 0; i < (1 << 19); i++)
	{
		printf("%i,", noise[i]);
	}
	std::cout << "]" << std::endl;
#endif
	delete[] noise;

	return 0;
}






