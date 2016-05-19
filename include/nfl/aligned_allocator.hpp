
/*
 * Aligned memory allocator
 *
 * Adapted from
 * http://jmabille.github.io/blog/2014/12/06/aligned-memory-allocator/
 */

namespace nfl {

namespace detail {
inline void* _aligned_malloc(size_t size, size_t alignment) {
  void* res = 0;
  void* ptr = malloc(size + alignment);
  if (ptr != 0) {
    res = reinterpret_cast<void*>(
        (reinterpret_cast<size_t>(ptr) & ~(size_t(alignment - 1))) + alignment);
    *(reinterpret_cast<void**>(res) - 1) = ptr;
  }
  return res;
}
}

inline void* aligned_malloc(size_t size, size_t alignment) {
#if MALLOC_ALREADY_ALIGNED
  return malloc(size);
#elif HAS_MM_MALLOC
  return _mm_malloc(size, alignment);
#elif HAS_POSIX_MEMALIGN
  void* res;
  const int failed = posix_memalign(&res, size, alignment);
  if (failed) res = 0;
  return res;
#elif (defined _MSC_VER)
  return _aligned_malloc(size, alignment);
#else
  return detail::_aligned_malloc(size, alignment);
#endif
}

namespace detail {
inline void _aligned_free(void* ptr) {
  if (ptr != 0) free(*(reinterpret_cast<void**>(ptr) - 1));
}
}

inline void aligned_free(void* ptr) {
#if MALLOC_ALREADY_ALIGNED
  free(ptr);
#elif HAS_MM_MALLOC
  _mm_free(ptr);
#elif HAS_POSIX_MEMALIGN
  free(ptr);
#elif defined(_MSC_VER)
  _aligned_free(ptr);
#else
  detail::_aligned_free(ptr);
#endif
}

template <class T, int N>
class aligned_allocator {
 public:
  typedef T value_type;
  typedef T& reference;
  typedef const T& const_reference;
  typedef T* pointer;
  typedef const T* const_pointer;
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;

  template <class U>
  struct rebind {
    typedef aligned_allocator<U, N> other;
  };

  inline aligned_allocator() throw() {}
  inline aligned_allocator(const aligned_allocator&) throw() {}

  template <class U>
  inline aligned_allocator(const aligned_allocator<U, N>&) throw() {}

  inline ~aligned_allocator() throw() {}

  inline pointer address(reference r) { return &r; }
  inline const_pointer address(const_reference r) const { return &r; }

  pointer allocate(size_type n,
                   typename std::allocator<void>::const_pointer hint = 0);
  inline void deallocate(pointer p, size_type);

  template <typename... Args>
  inline void construct(pointer p, Args&&... args) {
    new (p) value_type( std::forward<Args>(args)... );
  }
  inline void destroy(pointer p) { p->~value_type(); }

  inline size_type max_size() const throw() {
    return size_type(-1) / sizeof(T);
  }

  inline bool operator==(const aligned_allocator&) { return true; }
  inline bool operator!=(const aligned_allocator& rhs) {
    return !operator==(rhs);
  }
};

template <class T, int N>
typename aligned_allocator<T, N>::pointer aligned_allocator<T, N>::allocate(
    size_type n, typename std::allocator<void>::const_pointer hint) {
  pointer res = reinterpret_cast<pointer>(aligned_malloc(sizeof(T) * n, N));
  if (res == 0) throw std::bad_alloc();
  return res;
}

template <class T, int N>
inline void aligned_allocator<T, N>::deallocate(pointer p, size_type) {
  aligned_free(p);
}

}
