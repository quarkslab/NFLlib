cmake_minimum_required(VERSION 2.8)
project(nfl CXX ASM-ATT)
enable_testing()
if(APPLE)
    cmake_policy(SET CMP0042 NEW)
endif()

option(NFL_OPTIMIZED "Build an optimized version of NFL" OFF)
if(NFL_OPTIMIZED)
    message(STATUS "Tring to use optimized version of NFL")
    add_definitions(-DNFL_OPTIMIZED)

    include(CheckCXXSourceCompiles)
    set(CMAKE_CXX_FLAGS_OLD ${CMAKE_CXX_FLAGS})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

    check_cxx_source_compiles("#include <immintrin.h>\nint main() { return 0;}" NFLLIB_USE_AVX)
    check_cxx_source_compiles("#include <pmmintrin.h>\nint main() { return 0;}" NFLLIB_USE_SSE)
    if(NFLLIB_USE_AVX)
        message(STATUS "Using AVX vector engine")
        add_definitions(-DNTT_AVX2)
    elseif(NFLLIB_USE_SSE)
        message(STATUS "Using SSE vector engine")
        add_definitions(-DNTT_SSE)
    else()
        message(ERROR "Asking for optimized mode but neither SSE not AVX support detected")
    endif()
endif()

include(CMakeCompilers.txt)


set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)

include_directories(include include/nfl include/nfl/arch include/nfl/opt include/nfl/opt/arch include/nfl/prng)


file(GLOB NFLLIB_ASMS lib/prng/*.s)
file(GLOB NFLLIB_SRCS lib/prng/*.cpp lib/params/*.cpp)
file(GLOB NFLLIB_HDRS_0 include/nfl.hpp)
file(GLOB NFLLIB_HDRS_1 include/nfl/*.hpp)
file(GLOB NFLLIB_HDRS_2 include/nfl/arch/*.hpp)
file(GLOB NFLLIB_HDRS_3 include/nfl/opt/*.hpp)
file(GLOB NFLLIB_HDRS_4 include/nfl/opt/arch/*.hpp)
file(GLOB NFLLIB_PRNG_HDRS include/nfl/prng/*.hpp include/nfl/prng/*.h)

# shared library
add_library(nfllib SHARED ${NFLLIB_ASMS} ${NFLLIB_SRCS})
target_link_libraries(nfllib ${GMP_LIBRARY} ${MPFR_LIBRARY})
install(TARGETS nfllib LIBRARY DESTINATION lib)

# static library
add_library(nfllib_static STATIC ${NFLLIB_ASMS} ${NFLLIB_SRCS})
target_link_libraries(nfllib_static ${GMP_LIBRARY} ${MPFR_LIBRARY})
install(TARGETS nfllib_static ARCHIVE DESTINATION lib)

install(FILES ${NFLLIB_} DESTINATION include)
install(FILES ${NFLLIB_HDRS_0} DESTINATION include)
install(FILES ${NFLLIB_HDRS_1} DESTINATION include/nfl)
install(FILES ${NFLLIB_HDRS_2} DESTINATION include/nfl/arch)
install(FILES ${NFLLIB_HDRS_3} DESTINATION include/nfl/opt)
install(FILES ${NFLLIB_HDRS_4} DESTINATION include/nfl/opt/arch)
install(FILES ${NFLLIB_PRNG_HDRS} DESTINATION include/nfl/prng)

add_subdirectory(tests)
