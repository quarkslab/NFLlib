set(CONFIGS
    "8, 60, uint32_t"
    "128, 14, uint16_t"
    "1024, 60, uint32_t"
    "8192, 124, uint64_t"
    "32768, 124, uint64_t"
)

#always double check during tests
add_definitions(-DCHECK_STRICTMOD)

#if cereal exists
find_package(CEREAL)
if(CEREAL_FOUND)
MESSAGE(STATUS "Found CEREAL: ${CEREAL_INCLUDE_DIR}")
include_directories(${CEREAL_INCLUDE_DIR})
endif()

foreach(CONFIG ${CONFIGS})
    string(REGEX REPLACE ", " "_" CONFIG_ID ${CONFIG})

    foreach(VERSION op func)

        file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/nfllib_demo${CONFIG_ID}_${VERSION}.cpp"
            "#define CONFIG ${CONFIG}\n#include \"${CMAKE_CURRENT_SOURCE_DIR}/nfllib_demo_main_${VERSION}.cpp\"\n")
        add_executable(nfllib_demo${CONFIG_ID}_${VERSION}
                       EXCLUDE_FROM_ALL
                       ${CMAKE_CURRENT_BINARY_DIR}/nfllib_demo${CONFIG_ID}_${VERSION}.cpp)
        target_link_libraries(nfllib_demo${CONFIG_ID}_${VERSION} nfllib)
        add_test(build_nfllib_demo${CONFIG_ID}_${VERSION} "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target  nfllib_demo${CONFIG_ID}_${VERSION})
        add_test(run_nfllib_demo${CONFIG_ID}_${VERSION} nfllib_demo${CONFIG_ID}_${VERSION})
        set_tests_properties(run_nfllib_demo${CONFIG_ID}_${VERSION} PROPERTIES DEPENDS nfllib_demo${CONFIG_ID}_${VERSION})

    endforeach()

    foreach(OP nfl_add nfl_sub nfl_mul nfl_eq nfl_neq nfl_stream poly_p poly_set poly_mpz poly_serialize_manually)

        file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${OP}${CONFIG_ID}.cpp"
            "#define CONFIG ${CONFIG}\n#include \"${CMAKE_CURRENT_SOURCE_DIR}/${OP}.cpp\"\n")
        add_executable(${OP}${CONFIG_ID}
                       EXCLUDE_FROM_ALL
                       ${CMAKE_CURRENT_BINARY_DIR}/${OP}${CONFIG_ID}.cpp)
        target_link_libraries(${OP}${CONFIG_ID} nfllib)
        add_test(build_${OP}${CONFIG_ID} "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${OP}${CONFIG_ID})
        add_test(run_${OP}${CONFIG_ID} ${OP}${CONFIG_ID})
        set_tests_properties(run_${OP}${CONFIG_ID} PROPERTIES DEPENDS build_${OP}${CONFIG_ID})

    endforeach()

    if(CEREAL_FOUND)

        file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/poly_serialize_cereal${CONFIG_ID}.cpp"
            "#define CONFIG ${CONFIG}\n#include \"${CMAKE_CURRENT_SOURCE_DIR}/poly_serialize_cereal.cpp\"\n")
        add_executable(poly_serialize_cereal${CONFIG_ID}
                       EXCLUDE_FROM_ALL
                       ${CMAKE_CURRENT_BINARY_DIR}/poly_serialize_cereal${CONFIG_ID}.cpp)
        target_link_libraries(poly_serialize_cereal${CONFIG_ID} nfllib)
        add_test(build_poly_serialize_cereal${CONFIG_ID} "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target  poly_serialize_cereal${CONFIG_ID})
        add_test(run_poly_serialize_cereal${CONFIG_ID} poly_serialize_cereal${CONFIG_ID})
        set_tests_properties(run_poly_serialize_cereal${CONFIG_ID} PROPERTIES DEPENDS poly_serialize_cereal${CONFIG_ID})

    endif()

endforeach()

add_executable(prng_demo
               EXCLUDE_FROM_ALL
               prng_demo_main.cpp)
target_link_libraries(prng_demo nfllib)
add_test(build_prng_demo "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target prng_demo)
add_test(run_prng_demo prng_demo)
set_tests_properties(run_prng_demo PROPERTIES DEPENDS build_prng_demo)

add_executable(ntt_perfs EXCLUDE_FROM_ALL ntt_perfs.cpp)
target_link_libraries(ntt_perfs nfllib)

add_executable(ntt_multi multi0.cpp multi1.cpp)
