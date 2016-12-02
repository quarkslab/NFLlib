#ifndef NFL_TESTS_TOOLS_H
#define NFL_TESTS_TOOLS_H

#include <chrono>

template <class T, size_t Align, class... Args>
T* alloc_aligned(size_t n, Args&& ... args)
{
	T* ret;
	if (posix_memalign((void**) &ret, Align, sizeof(T)*n) != 0) {
		throw std::bad_alloc();
	}
	for (size_t i = 0; i < n; i++) {
		new (&ret[i]) T(std::forward<Args>(args)...);
	}
	return ret;
}

template <class T>
void free_aligned(size_t n, T* p)
{
	for (size_t i = 0; i < n; i++) {
		p[i].~T();
	}
	free(p);
}

template <class T>
double get_time_us(T const& start, T const& end, uint32_t N)
{
  auto diff = end-start;
  return (long double)(std::chrono::duration_cast<std::chrono::microseconds>(diff).count())/N;
}

#endif
